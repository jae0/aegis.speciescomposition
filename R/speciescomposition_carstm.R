
speciescomposition_carstm = function( p=NULL, DS="parameters", redo=FALSE, ... ) {

  if (DS=="parameters") {

    if (is.null(p)) {
      p = aegis.speciescomposition::speciescomposition_parameters(...)
    } else {

      # if p is passed, assume it is a secondary call ... overwrite nonrelevent params to force use of project defaults
      p$data_root = NULL
      p$datadir  = NULL
      p$data_transformation= NULL
      p$carstm_modelcall = NULL  # defaults to generic
      p$carstm_model_tag = NULL

      p = aegis.speciescomposition::speciescomposition_parameters(p=p, ...)
    }

    p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.substrate", "aegis.coastline", "aegis.polygons", "aegis.temperature", "aegis.speciescomposition", "aegis.survey", "carstm" ) ) )


    if ( !exists("areal_units_source", p)) p$areal_units_source = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( p$spatial_domain == "SSE" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "groundfish_strata" #.. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      if ( !exists("inputdata_temporal_discretization_yr", p)) p$inputdata_temporal_discretization_yr = 1/12  #  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( p$spatial_domain == "snowcrab" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      if ( !exists("inputdata_spatial_discretization_planar_km", p))p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      if ( !exists("inputdata_temporal_discretization_yr", p)) p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }
    }

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )

        if ( !exists("carstm_model_label", p)) p$carstm_model_label = "production"

        p$carstm_modelcall = paste(
          'inla( formula = ', p$variabletomodel,
          ' ~ 1
            + f( dyri, model="ar1", hyper=H$ar1 )
            + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
            family = "normal",
            data= M,
            control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, config=TRUE),
            control.results = list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor = list(compute=FALSE, link=1 ),
            # control.fixed = list(prec.intercept = 0.1),
            control.fixed = H$fixed,  # priors for fixed effects, generic is ok
            # control.inla = list(cmin = 0 ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.inla = list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            control.inla = list( cmin = 0, h=1e-4, tolerance=1e-9, strategy="adaptive", optimise.strategy="smart"), # restart=3), # restart a few times in case posteriors are poorly defined
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            verbose=TRUE
          )'
        )
      }
        #    + f(tiyr, model="ar1", hyper=H$ar1 )
        # + f(year,  model="ar1", hyper=H$ar1 )

      if ( grepl("glm", p$carstm_modelengine) ) {
        if ( !exists("carstm_model_label", p)) p$carstm_model_label = "default_glm"
        p$carstm_modelcall = paste(
          'glm( formula =',  p$variabletomodel,
          ' ~ 1 + AUID + t + z + substrate.grainsize +tiyr,
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }

      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv" ) ) )
        if ( !exists("carstm_model_label", p)) p$carstm_model_label = "default_gam"
         p$carstm_modelcall = paste(
          'gam( formula =',  p$variabletomodel,
          ' ~ 1 + AUID + s(t) + s(z) + s(substrate.grainsize) + s(year) + s(dyear),
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }
    }

    p = carstm_parameters( p=p )  #generics

    p$carstm_inputs_aggregated = FALSE

    return(p)
  }


  # ------------------------

  if ( DS=="carstm_inputs") {

    # prediction surface
    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))
    sppoly = areal_units( p=p )  # will redo if not found
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    fn = carstm_filenames( p=p, projectname="speciescomposition", projecttype="carstm_inputs", areal_units_fn=areal_units_fn )

    fn = file.path( p$modeldir, fn)

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")


    # do this immediately to reduce storage for sppoly (before adding other variables)

    M = speciescomposition_db( p=p, DS="speciescomposition"  )

    # globally remove all unrealistic data
      # p$quantile_bounds_data = c(0.0005, 0.9995)
    if (exists("quantile_bounds_data", p)) {
      TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds_data, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
      keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
      if (length(keep) > 0 ) M = M[ keep, ]
        # this was -1.7, 21.8 in 2015
    }

    M = planar2lonlat(M, proj.type=p$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]


    M$AUID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
    M = M[ which(!is.na(M$AUID)),]

    names(M)[which(names(M)=="yr") ] = "year"
    M = M[ which(M$year %in% p$yrs), ]
    M$tiyr = lubridate::decimal_date ( M$timestamp )
    M$dyear = M$tiyr - M$year


    # reduce size
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")
    # M$plon = floor(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
    # M$plat = floor(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km


    pB = bathymetry_carstm( p=p, DS="parameters", variabletomodel="z" )
    pS = substrate_carstm( p=p, DS="parameters", variabletomodel="substrate.grainsize" )
    pT = temperature_carstm( p=p, DS="parameters", variabletomodel="t" )

    if (!(exists(pB$variabletomodel, M ))) M[,pB$variabletomodel] = NA
    if (!(exists(pS$variabletomodel, M ))) M[,pS$variabletomodel] = NA
    if (!(exists(pT$variabletomodel, M ))) M[,pT$variabletomodel] = NA

    kk = which(!is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0 ) M[kk, pB$variabletomodel] = bathymetry_lookup( p=p, locs=M[kk, c("lon", "lat")], source_data_class="aggregated_rawdata" )

    kk = which(!is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0 ) M[kk, pS$variabletomodel] = substrate_lookup(  p=p, locs=M[kk, c("lon", "lat")], source_data_class="aggregated_rawdata" )

    kk = which(!is.finite(M[, pT$variabletomodel]))
    if (length(kk) > 0 ) M[kk, pT$variabletomodel] = temperature_lookup(  p=p, locs=M[kk, c("lon", "lat")], timestamp=M$timestamp, source_data_class="aggregated_rawdata" )


    # if any still missing then use a mean depth by AUID
    kk =  which( !is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0) {
      AD = bathymetry_db ( p=pB, DS="aggregated_data"   )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( AD[, paste(pB$variabletomodel, "mean", sep="." )], AD$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      M[kk, pB$variabletomodel] = oo[jj ]
    }

    # if any still missing then use a mean substrate by AUID
    kk =  which( !is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0) {
      AD = substrate_db ( p=pS, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( AD[, paste(pS$variabletomodel, "mean", sep="." )], AD$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      M[kk, pS$variabletomodel] = oo[jj ]
    }


    # substrate coverage poor .. add from modelled results
    kk =  which( !is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0) {
      SI = carstm_summary( p=pS )
      jj = match( as.character( M$AUID[kk]), as.character( SI$AUID) )
      M[kk, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    }


    # if any still missing then use a mean temp by AUID
    kk =  which( !is.finite(M[, pT$variabletomodel]))
    if (length(kk) > 0) {
      AD = temperature_db ( p=pT, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")

      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      AD$uid = paste(AD$AUID, AD$year, AD$dyear, sep=".")

      M_dyear_discret = discretize_data( M$dyear, p$discretization$dyear )  # AD$dyear is discretized. . match discretization
      M$uid =  paste(M$AUID, M$year, M_dyear_discret, sep=".")

      oo = tapply( AD[, paste(pT$variabletomodel, "mean", sep="." )], AD$uid, FUN=median, na.rm=TRUE )

      jj = match( as.character( M$uid[kk]), as.character( names(oo )) )
      M[kk, pT$variabletomodel] = oo[jj ]
    }


    if( exists("spatial_domain", p)) M = geo_subset( spatial_domain=p$spatial_domain, Z=M ) # need to be careful with extrapolation ...  filter depths


    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL

    M = M[ which(is.finite(M[, pB$variabletomodel] )),]
    M = M[ which(is.finite(M[, pS$variabletomodel] )),]
    M = M[ which(is.finite(M[, pT$variabletomodel] )),]

    M$tag = "observations"


    APS = as.data.frame(sppoly)
    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA

    BI = carstm_summary( p=pB )
    jj = match( as.character( APS$AUID), as.character( BI$AUID) )
    APS[, pB$variabletomodel] = BI[[ paste(pB$variabletomodel,"predicted",sep="." ) ]] [jj]
    jj =NULL
    BI = NULL

    SI = carstm_summary( p=pS )
    jj = match( as.character( APS$AUID), as.character( SI$AUID) )
    APS[, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    jj =NULL
    SI = NULL

    # to this point APS is static, now add time dynamics (teperature)
    # ---------------------

    vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "AUID" )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")
    APS$year = floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$year


    TI = carstm_summary( p=pT )
    TI = TI[[ paste(pT$variabletomodel,"predicted",sep="." )]]

    auid_map = match( APS$AUID, dimnames(TI)$AUID )
    year_map = match( as.character(APS$year), dimnames(TI)$year )

    dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
    dyear_map = as.numeric( cut( APS$dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )

    dindex = cbind(auid_map, year_map, dyear_map )

    APS[, pT$variabletomodel] = TI[ dindex]
    jj =NULL
    TI = NULL

    M = rbind( M[, names(APS)], APS )
    APS = NULL

    M$auid  = as.numeric( factor(M$AUID) )

    M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
    M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

    M$tiyr  = floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

    M$year = floor( M$tiyr)
    M$year_factor = as.numeric( factor( M$year, levels=p$yrs))
    M = M[ is.finite(M$year_factor), ]
    M$dyear =  M$tiyr - M$year   # revert dyear to non-discretized form

    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)

    save( M, file=fn, compress=TRUE )
    return( M )
  }


} # end function



