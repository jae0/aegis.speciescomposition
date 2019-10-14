
speciescomposition_carstm = function( p=NULL, DS=NULL, sppoly=NULL, redo=FALSE, ... ) {

  #\\ Note inverted convention: depths are positive valued
  #\\ i.e., negative valued for above sea level and positive valued for below sea level

  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  # -----------------------


  if ( DS=="carstm_inputs") {

    fn = file.path( p$modeldir, paste( "speciescomposition", "carstm_inputs", p$auid,
      p$inputdata_spatial_discretization_planar_km,
      round(p$inputdata_temporal_discretization_yr, 6),
      "rdata", sep=".") )

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found

    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))

    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = speciescomposition.db( p=p, DS="speciescomposition" )

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area

    M$StrataID = as.character( M$StrataID )  # match each datum to an area
    M$tiyr = M$yr + M$dyear
    M[,p$variabletomodel] = M$speciescomposition.mean
    M$$tag = "observations"

    APS = as.data.frame(sppoly)
    APS$StrataID = as.character( APS$StrataID )
    APS$tag ="predictions"
    APS$[,p$variabletomodel] = NA


    pb = aegis.bathymetry::bathymetry_parameters( p=p, project_class =="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
    BI = bathymetry_carstm ( p=pb, DS="carstm_modelled" )  # unmodeled!
    jj = match( as.character( APS$StrataID), as.character( BI$StrataID) )
    APS$z = BI$z.predicted[jj]
    jj =NULL
    BI = NULL

    ps = aegis.substrate::substrate_parameters( p=p, project_class =="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
    SI = substrate_carstm ( p=ps, DS="carstm_modelled" )  # unmodeled!
    jj = match( as.character( APS$StrataID), as.character( SI$StrataID) )
    APS$substrate.grainsize = SI$substrate.grainsize.predicted[jj]
    jj =NULL
    SI = NULL

    pt = aegis.temperature::temperature_parameters( p=p, project_class =="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
    TI = temperature_carstm ( p=pt, DS="carstm_modelled" )  # unmodeled!
    jj = match( as.character( APS$StrataID), as.character( TI$StrataID) )
    APS$temperature = TI$temperature.predicted[jj]
    jj =NULL
    TI = NULL

    vn = c("t", "tag", "StrataID", "z", "substrate.grainsize", "temperature")
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")

    M$tiyr = M$yr + M$dyear
    M = rbind( M[, names(APS)], APS )
    APS = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    sppoly = NULL

    save( M, file=fn, compress=TRUE )
    return( M )
  }


  # ------------



  if ( DS=="carstm_modelled") {

    auids = paste(  p$auid, p$inputdata_spatial_discretization_planar_km,
      round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )

    fn = file.path( p$modeldir, paste("speciescomposition", "carstm_modelled", p$carstm_modelengine, auids, "rdata", sep="." ) )
    fn_fit = file.path( p$modeldir, paste( "speciescomposition", "carstm_modelled_fit", p$carstm_modelengine, auids,  "rdata", sep=".") )

    if (!redo)  {
         if (file.exists(fn)) {
        load( fn)
        return( res )
      }
      if (DS=="carstm_modelled_fit") {
        if (file.exists(fn_fit)) {
          load( fn_fit )
          return( fit )
        }
      }
    }

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    res = list(StrataID = sppoly[["StrataID"]])  # init results list
    res$strata = as.numeric(res$StrataID)

    M = speciescomposition_carstm ( p=p, DS="carstm_inputs" )  # 16 GB in RAM just to store!
    M$strata  = as.numeric( M$StrataID)

    fit  = NULL

    if ( grepl("glm", p$carstm_modelengine) ) {

      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      # s = summary(fit)
      # AIC(fit)  # 104487274
      # reformat predictions into matrix form
      ii = which( M$tag=="predictions" & M$strata %in% M[ which(M$tag=="observations"), "strata"] )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2

      res[,paste(p$variabletomodel, "predicted", sep=".")] =  preds$fit
      res[,paste(p$variabletomodel, "predicted_se", sep=".")] =  preds$se.fit
      res[,paste(p$variabletomodel, "predicted_lb", sep=".")] =  preds$fit - preds$se.fit
      res[,paste(p$variabletomodel, "predicted_ub", sep=".")] =  preds$fit + preds$se.fit
      save( res, file=fn, compress=TRUE )
    }

    if ( grepl("gam", p$carstm_modelengine) ) {
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      s = summary(fit)
      AIC(fit)  # 104487274
      # reformat predictions into matrix form
      ii = which( M$tag=="predictions" & M$strata %in% M[ which(M$tag=="observations"), "strata"] )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      res[,paste(p$variabletomodel, "predicted", sep=".")] =  preds$fit
      res[,paste(p$variabletomodel, "predicted_se", sep=".")] =  preds$se.fit
      res[,paste(p$variabletomodel, "predicted_lb", sep=".")] =  preds$fit - preds$se.fit
      res[,paste(p$variabletomodel, "predicted_ub", sep=".")] =  preds$fit + preds$se.fit
      save( res, file=fn, compress=TRUE )
    }





    if ( grepl("inla", p$carstm_modelengine) ) {

      H = carstm_hyperparameters( sd(M[,p$variabletomodel], na.rm=TRUE), alpha=0.5, median( M[,p$variabletomodel], na.rm=TRUE) )

      M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

      M$zi = discretize_data( M$z, p$discretization$z )
      M$year = floor(M$tiyr)
      M$dyear  =  M$tiyr - M$year
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      # reformat predictions into matrix form
      ii = which(
        M$tag=="predictions" &
        M$strata %in% res$strata &
        M$year %in% p$yrs
      )  # filter by strata and years in case additional data in other areas and times are used in the input data

      matchfrom = list( strata=M$strata[ii], year=as.character(M$year[ii]), dyear=M$dyear[ii] )
      matchto   = list( strata=res$strata, year=as.character(p$yrs), dyear=factor(p$dyears) )

      vn = paste( p$variabletomodel, "predicted", sep=".")
      res[vn] = reformat_to_array( input=fit$summary.fitted.values[ ii, "mean" ], matchfrom=matchfrom, matchto=matchto )

      vn = paste( p$variabletomodel, "predicted_lb", sep=".")
      res[vn] = reformat_to_array( input=fit$summary.fitted.values[ ii, "0.025quant" ], matchfrom=matchfrom, matchto=matchto )

      vn = paste( p$variabletomodel, "predicted_ub", sep=".")
      res[vn] = reformat_to_array( input=fit$summary.fitted.values[ ii, "0.975quant" ], matchfrom=matchfrom, matchto=matchto )

      # random effects results ..
      if (exists("summary.random", fit)) {

        nstrata = length(res$StrataID)

        # reformat predictions into matrix form
        ii = which(
          M$tag=="predictions" &
          M$strata %in% res$strata &
          M$year %in% p$yrs
        )  # filter by strata and years in case additional data in other areas and times are used in the input data


        if (exists("iid_error", fit$summary.random)) {
          # IID random effects
          vn = paste( p$variabletomodel, "random_sample_iid", sep=".")
          matchfrom = list( strata=M$strata[ii], year=M$year[ii], dyear=M$dyear[ii] )
          matchto   = list( strata=res$strata, year=p$yrs, dyear=factor(p$dyears) )
          res[vn] = reformat_to_array( input=fit$summary.random$iid_error[ ii, "mean" ], matchfrom=matchfrom, matchto=matchto )
          # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="1950", dyear="0") )
        }

        if (exists("strata", fit$summary.random)) {

          if (nrow(fit$summary.random$strata) == nstrata*2) {
            # CAR random effects can be of variable length depending upon model construct:

            # a single spatial and nonspatial effect (no grouping across time)
            jj = 1:nstrata
            matchfrom = list( strata=fit$summary.random$strata$ID[jj]  )
            matchto   = list( strata=res$strata  )
            vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
            res[vn] = reformat_to_array( fit$summary.random$strata[ jj, "mean" ], matchfrom=matchfrom, matchto=matchto )

            vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
            res[vn] = reformat_to_array( fit$summary.random$strata[ jj+max(jj), "mean" ], matchfrom=matchfrom, matchto=matchto )
            # carstm_plot( p=p, res=res, vn=vn  )

          } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$ny ) {
            # spatial and nonspatial effects grouped by year

            matchfrom = list( strata=M$strata[ii], year=M$year[ii] )
            matchto   = list( strata=res$strata, year=p$yrs )

            vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
            res[vn] = reformat_to_array( fit$summary.random$strata[ ii, "mean" ], matchfrom=matchfrom, matchto=matchto )

            vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
            res[vn] = reformat_to_array( fit$summary.random$strata[ ii+max(ii), "mean" ], matchfrom=matchfrom, matchto=matchto )
            # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000" ) )

          } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$nt ) {

            # need to test/fix ...
            matchfrom = list( StrataID=M$StrataID[ii], year=M$year[ii], dyear=M$dyear[ii] ),
            matchto   = list( StrataID=res$StrataID, year=p$yrs, dyear=factor(p$dyears) )

            vn = paste( p$variabletomodel, "random_strata_nonspatial", sep=".")
            res[vn] = reformat_to_array( fit$summary.random$strata[ ii, "mean" ], matchfrom=matchfrom, matchto=matchto )

            vn = paste( p$variabletomodel, "random_strata_spatial", sep=".")
            res[vn] = reformat_to_array( fit$summary.random$strata[ ii+max(ii), "mean" ], matchfrom=matchfrom, matchto=matchto )
            # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8") )

          }
        }
      save( res, file=fn, compress=TRUE )
    }

    return( res )
  }

}
