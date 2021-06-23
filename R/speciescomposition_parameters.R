

speciescomposition_parameters = function( p=list(), project_name="speciescomposition", project_class="core",  ... ) {


  # ---------------------
  # deal with additional passed parameters
  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args


  # ---------------------
  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "lubridate",  "lattice",
    "parallel", "sf", "GADMTools", "INLA" , "data.table") ) )

  p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.coastline",
    "aegis.polygons", "aegis.substrate", "aegis.temperature", "aegis.survey", "aegis.speciescomposition" ) ) )


  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )


  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )


  p = parameters_add_without_overwriting( p,
    spatial_domain = "SSE",  # canada.east.highres and canada.east.superhighres result in memory overflow
    spatial_domain_subareas = c( "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    aegis_dimensionality="space-year"
  )

  p$quantile_bounds =c(0, 0.95) # trim upper bounds

  p = spatial_parameters( p=p)

  # define focal years for modelling and interpolation

  if (!exists("year.assessment", p )) {
    message("need probably want to assign current year.assessment, using current year for now")
    p$year.assessment = lubridate::year(lubridate::now())
  }

  p = parameters_add_without_overwriting( p, yrs=1999:p$year.assessment, timezone="America/Halifax" )  # default
  p = temporal_parameters(p=p)

  p = parameters_add_without_overwriting( p,
    additional.data=c("groundfish", "snowcrab", "USSurvey_NEFSC", "lobster"),
    taxa =  "maxresolved",
    varstomodel = c( "pca1", "pca2", "ca1", "ca2" ),
    inputdata_spatial_discretization_planar_km = p$pres/2, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
    inputdata_temporal_discretization_yr = 1/12  # ie., controls resolution of data prior to modelling to reduce data set and speed up modelling;; use 1/12 -- monthly or even 1/4.. if data density is low
  )



  # ---------------------

  if (project_class=="core") {
    p$project_class = "core"
    return(p)
  }

  # ---------------------

  if (project_class=="carstm") {
    # simple run of carstm. There are two types:
    #   one global, run directly from  polygons defined in aegis.bathymetry/inst/scripts/99.bathymetry.carstm.R
    #   and one that is called secondarily specific to a local project's polygons (eg. snow crab)
    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )
    p$project_class = "carstm"

    if (!exists("variabletomodel", p)) stop( "The dependent variable, p$variabletomodel needs to be defined")

    # over-rides
    p$inputdata_spatial_discretization_planar_km = 0.5  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
    p$inputdata_temporal_discretization_yr = 1/52  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling


    # defaults in case not provided ...
    p = parameters_add_without_overwriting( p,
      areal_units_xydata = "speciescomposition_db(p=p, DS='areal_units_input')",
      areal_units_type = "tesselation", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_resolution_km = 1, # default in case not provided ... 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
      areal_units_proj4string_planar_km =  p$aegis_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      areal_units_overlay = "none",
      areal_units_timeperiod = "none",
      tus="yr",
      fraction_todrop = 1/11,
      fraction_cv = 1.0,
      fraction_good_bad = 0.9,
      areal_units_constraint_nmin=3,  # best compromise
      areal_units_constraint_ntarget=20,
      nAU_min = 30,
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "default",
      carstm_inputs_prefilter = "rawdata"
    )


    if ( !exists("carstm_inputdata_model_source", p))  p$carstm_inputdata_model_source = list()
    p$carstm_inputdata_model_source = parameters_add_without_overwriting( p$carstm_inputdata_model_source,
      bathymetry = "stmv",  # "stmv", "hybrid", "carstm"
      substrate = "stmv",  # "stmv", "hybrid", "carstm"
      temperature = "carstm"  # "stmv", "hybrid", "carstm"
    )


    if ( grepl("inla", p$carstm_modelengine) ) {
      if ( !exists("carstm_model_formula", p)  ) {

        p$carstm_model_formula = as.formula( paste(
         p$variabletomodel, ' ~ 1',
            ' + f( uid, model="iid" ) ',
            ' + f( season, model="rw2", hyper=H$rw2, cyclic=TRUE ) ',
            ' + f( time, model="ar1",  hyper=H$ar1 ) ',
            ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE ) ',
            ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
            ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
#             ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
            ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group))'
          ) )
      }

      if ( !exists("carstm_model_family", p)  )  p$carstm_model_family = "gaussian"
    }

    p = carstm_parameters( p=p )  #generics

    if ( p$inputdata_spatial_discretization_planar_km >= p$areal_units_resolution_km ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$areal_units_resolution_km " )
    }

    message ("p$areal_units_resolution_km: ", p$areal_units_resolution_km)

    return(p)
  }


  # ---------------------

  if (project_class %in% c("stmv", "default")) {

    p$libs = c( p$libs, project.library ( "stmv" ) )
    p$project_class = "stmv"

    p = parameters_add_without_overwriting( p,
      stmv_model_label="default",
      stmv_variables = list(
        LOCS=c("plon", "plat"),
        TIME="tiyr"
      ),  # required as fft has no formulae
      inputdata_spatial_discretization_planar_km = p$pres / 4, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
      inputdata_temporal_discretization_yr = 1/12,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="fft",
      stmv_variogram_method = "fft",
      stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
      stmv_rsquared_threshold = 0.01, # lower threshold  .. ignore
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_force_complete_method = "linear_interp"
    )


    p = parameters_add_without_overwriting( p,
      stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1/2, 5 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
      stmv_distance_scale = p$stmv_distance_statsgrid * c( 1, 2, 3, 4, 5, 10, 15), # km ... approx guesses of 95% AC range
      stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 5, 10, 15 ),  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2, 3, 4, 8) # finalizing preds using linear interpolation
    )

    if (!exists("stmv_global_modelformula", p)) p$stmv_global_modelformula = "none"


      if (p$stmv_local_modelengine == "fft"  |  p$stmv_twostep_space == "fft"  ) {
        nu = 0.5  # exponential smoothing
        ac_local = 0.1  # ac at which to designate "effective range"
        p = parameters_add_without_overwriting( p,
          stmv_fft_filter = "matern tapered lowpass modelled fast_predictions", #  act as a low pass filter first before matern with taper
          stmv_autocorrelation_fft_taper = 0.9,  # benchmark from which to taper
          stmv_autocorrelation_localrange = ac_local,  # for output to stats
          stmv_autocorrelation_interpolation = c(0.25, 0.1, 0.05, 0.01),
          stmv_lowpass_nu = nu, # exp
          stmv_lowpass_phi = stmv::matern_distance2phi( distance=p$pres/2, nu=nu, cor=ac_local )
        )
      }


    # intervals of decimal years... fractional year breaks finer than the default 10 units (taking daily for now..)
    #.. need to close right side for "cut" .. controls resolution of data prior to modelling
    if (!exists("dyear_discretization_rawdata", p)) p$dyear_discretization_rawdata = c( {c(1:365)-1}/365, 1)


    # default to serial mode
    p = parameters_add_without_overwriting( p,
      stmv_runmode = list(
        globalmodel = FALSE,
        scale = rep("localhost", 1),
        interpolate_correlation_basis = list(
          cor_0.25 = rep("localhost", 1),
          cor_0.1  = rep("localhost", 1),
          cor_0.05 = rep("localhost", 1),
          cor_0.01 = rep("localhost", 1)
        ),
        interpolate_predictions = list(
          c1 = rep("localhost", 1),
          c2 = rep("localhost", 1),
          c3 = rep("localhost", 1),
          c4 = rep("localhost", 1),
          c5 = rep("localhost", 1),
          c6 = rep("localhost", 1),
          c7 = rep("localhost", 1)
        ),
        save_intermediate_results = TRUE,
        save_completed_data = TRUE # just a dummy variable with the correct name
      )
    )


    p = aegis_parameters(p=p, DS="stmv" ) # generics for aegis.* projects


    if ( p$inputdata_spatial_discretization_planar_km >= p$pres ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$pres " )
    }
#     message ("p$stmv_distance_statsgrid: ", p$stmv_distance_statsgrid)


    return(p)
  }


  # ---------------------

  if (project_class %in% c("hybrid")) {

    if (!exists("variabletomodel", p)) stop( "The dependent variable, p$variabletomodel needs to be defined")

    p = parameters_add_without_overwriting( p,
      stmv_model_label="default",
      stmv_variables = list(
        LOCS=c("plon", "plat"),
        TIME="tiyr"
      ),  # required as fft has no formulae
      inputdata_spatial_discretization_planar_km = p$pres / 4, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
      inputdata_temporal_discretization_yr = 1/12,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="carstm",
      stmv_local_covariates_carstm = "",  # only model covariates
      stmv_local_all_carstm = "",  # ignoring au
      stmv_local_modelcall = paste(
        'inla(
          formula =', p$variabletomodel, ' ~ 1
            + f( uid, model="iid" )
            + f(space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
          family = "gaussian",
          data= dat,
          control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE),  # config=TRUE if doing posterior simulations
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list(compute=FALSE, link=1 ),
          control.fixed=H$fixed,  # priors for fixed effects, generic is ok
          verbose=FALSE
        ) '
      ),   # NOTE:: this is a local model call
      stmv_filter_depth_m = TRUE,
      stmv_local_model_distanceweighted = TRUE,
      stmv_rsquared_threshold = 0.01, # lower threshold  .. ignore
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_distance_prediction_limits =c( 3, 25 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
      stmv_distance_interpolation = c(  2.5 , 5, 10, 15, 20, 40, 80, 150, 200 ) , # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_runmode = list(
        carstm = rep("localhost", 1),
        globalmodel = FALSE,
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )  # ncpus for each runmode
    )


    p = parameters_add_without_overwriting( p,
      stmv_distance_prediction_limits = p$stmv_distance_statsgrid * c( 1, 2 ), # range of permissible predictions km (i.e  stats grid to upper limit based upon data density)
      stmv_distance_interpolation = p$stmv_distance_statsgrid * c( 1/2, 1, 2 ),  # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_distance_interpolate_predictions = p$stmv_distance_statsgrid * c( 1/2, 1, 2) # finalizing preds using linear interpolation
    )


    p = parameters_add_without_overwriting( p,
      stmv_runmode = list(
        carstm = rep("localhost", 1),
        globalmodel = FALSE,
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )
    )

    p = aegis_parameters( p=p, DS="stmv" )  # get defaults

    # intervals of decimal years... fractional year breaks finer than the default 10 units (taking daily for now..)
    #.. need to close right side for "cut" .. controls resolution of data prior to modelling
    if (!exists("dyear_discretization_rawdata", p)) p$dyear_discretization_rawdata = c( {c(1:365)-1}/365, 1)
    if ( p$inputdata_spatial_discretization_planar_km >= p$pres ) {
      warning( "p$inputdata_spatial_discretization_planar_km >= p$pres " )
    }
    message ("p$stmv_distance_statsgrid: ", p$stmv_distance_statsgrid)

    return(p)
  }

}
