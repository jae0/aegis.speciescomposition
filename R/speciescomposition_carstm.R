
  speciescomposition_carstm = function( DS="", p=NULL ) {

    ddir = project.datadirectory( "aegis", "speciescomposition" )
    dir.create( ddir, showWarnings=FALSE, recursive=TRUE )

    infix = paste( p$spatial_domain,  p$taxa, sep=".")


  # ---------------------

  if (DS=="carstm_auid") {
    # translate param values from one project to a unified representation
    # must be first to catch p
    P = speciescomposition_carstm(
      DS="parameters",
      project_class = "carstm", # defines which parameter class / set to load
      project_name = "speciescomposition",
      yrs = p$yrs,
      spatial_domain = p$spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_overlay = p$areal_units_overlay, # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = p$areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = p$areal_units_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = p$inputdata_temporal_discretization_yr,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      auid = p$auid
    )

    return(P)
  }


    # -----------------------




  if (DS=="parameters") {

    if ( !exists("project_name", p)) p$project_name = "speciescomposition"

    p = carstm_parameters( p=p, DS="basic" )  #generics

    if ( !exists("areal_units_strata_type", p)) p$areal_units_strata_type = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( p$spatial_domain == "SSE" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "groundfish_strata" #.. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      p$inputdata_temporal_discretization_yr = 1/12  #  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
    }

    if ( p$spatial_domain == "snowcrab" ) {
      if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
      if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
      # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }
    }

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla.default"  # {model engine}.{label to use to store}


    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )
        p$carstm_modelcall = paste(
          'inla( formula = ', p$variabletomodel,
          ' ~ 1
            + f(tiyr2, model="seasonal", season.length=10 )
            + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(gsi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
            + f(iid_error, model="iid", hyper=H$iid),
            family = "normal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            num.threads=4,
            #blas.num.threads=4,
            verbose=TRUE
          )'
        )
      }
        #    + f(tiyr, model="ar1", hyper=H$ar1 )
        # + f(year,  model="ar1", hyper=H$ar1 )

      if ( grepl("glm", p$carstm_modelengine) ) {
        p$carstm_modelcall = paste(
          'glm( formula =',  p$variabletomodel,
          ' ~ 1 + StrataID + t + z + substrate.grainsize +tiyr,
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }

      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv" ) ) )
        p$carstm_modelcall = paste(
          'gam( formula =',  p$variabletomodel,
          ' ~ 1 + StrataID + s(t) + s(z) + s(substrate.grainsize) + s(yr) + s(dyear),
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }
    }
    return(p)
  }


  # ------------------------


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

      M = M[ which(M$yr %in% p$yrs), ]
      M$tiyr = lubridate::decimal_date ( M$date )

      # globally remove all unrealistic data
  #    keep = which( M[,p$variabletomodel] >= -3 & M[,p$variabletomodel] <= 25 ) # hard limits
  #    if (length(keep) > 0 ) M = M[ keep, ]

      # p$quantile_bounds_data = c(0.0005, 0.9995)
      if (exists("quantile_bounds_data", p)) {
        TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
        keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
      }

      M$dyear = M$tiyr - M$yr
      M$dyear = discretize_data( M$dyear, seq(0, 1, by=p$inputdata_temporal_discretization_yr), digits=6 )

      # reduce size
      M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
      # levelplot(z.mean~plon+plat, data=M, aspect="iso")

      M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      M$lon = NULL
      M$lat = NULL
      M = M[ which(is.finite(M$StrataID)),]
      M$StrataID = as.character( M$StrataID )  # match each datum to an area

      names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel

      M$tag = "observations"

      APS = as.data.frame(sppoly)
      APS$StrataID = as.character( APS$StrataID )
      APS$tag ="predictions"
      APS[,p$variabletomodel] = NA


      pB = bathymetry_carstm( p=p, DS="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
      pS = substrate_carstm( p=p, DS ="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
      pT = temperature_carstm( p=p, DS ="carstm_auid" ) # transcribes relevant parts of p to load



      M$plon = round(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      M$plat = round(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km

      BI = carstm_model ( p=pB, DS="aggregated_data" )  # unmodeled!
      jj = match( as.character( M$StrataID), as.character( BI$StrataID) )
      M[, pB$variabletomodel] = BI[jj, paste(pB$variabletomodel,"predicted",sep="." )]
      jj =NULL
      BI = NULL

      SI = carstm_model ( p=pS, DS="aggregated_data" )  # unmodeled!
      jj = match( as.character( M$StrataID), as.character( SI$StrataID) )
      M[, pS$variabletomodel] = SI[jj, paste(pS$variabletomodel,"predicted",sep="." )]
      jj =NULL
      SI = NULL

      TI = carstm_model ( p=pT, DS="aggregated_data" )  # unmodeled!
      jj = match( as.character( M$StrataID), as.character( TI$StrataID) )
      M[, pT$variabletomodel] = TI[jj, paste(pT$variabletomodel,"predicted",sep="." )]
      jj =NULL
      TI = NULL



      BI = carstm_model ( p=pB, DS="carstm_modelled" )  # unmodeled!
      jj = match( as.character( APS$StrataID), as.character( BI$StrataID) )
      APS[, pB$variabletomodel] = BI[jj, paste(pB$variabletomodel,"predicted",sep="." )]
      jj =NULL
      BI = NULL

      SI = carstm_model ( p=pS, DS="carstm_modelled" )  # unmodeled!
      jj = match( as.character( APS$StrataID), as.character( SI$StrataID) )
      APS[, pS$variabletomodel] = SI[jj, paste(pS$variabletomodel,"predicted",sep="." )]
      jj =NULL
      SI = NULL

      TI = carstm_model ( p=pT, DS="carstm_modelled" )  # unmodeled!
      jj = match( as.character( APS$StrataID), as.character( TI$StrataID) )
      APS[, pT$variabletomodel] = TI[jj, paste(pT$variabletomodel,"predicted",sep="." )]
      jj =NULL
      TI = NULL

      vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel,  pT$variabletomodel, "tag", "StrataID" )
      APS = APS[, vn]

      # expand APS to all time slices
      n_aps = nrow(APS)
      APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
      names(APS) = c(vn, "tiyr")

      M$tiyr = M$yr + M$dyear
      M = rbind( M[, names(APS)], APS )
      APS = NULL

      M$strata  = as.numeric( M$StrataID)
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
      M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
      M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

      M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
      M$tiyr2 = M$tiyr  # a copy

      M$year = floor(M$tiyr)
      M$dyear  =  factor( as.character( trunc(  (M$tiyr - M$year )/ p$tres )*p$tres), levels=p$dyears)

      save( M, file=fn, compress=TRUE )
      return( M )
    }


  } # end function



