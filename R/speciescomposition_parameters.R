

speciescomposition_parameters = function( p=list(), project_name="speciescomposition", project_class="core",  ... ) {


  # ---------------------
  # deal with additional passed parameters
  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args


  # ---------------------
  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "lubridate",  "lattice",
    "parallel", "sf", "INLA" , "data.table") ) )

  p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.coastline",
    "aegis.polygons", "aegis.substrate", "aegis.temperature", "aegis.survey", "aegis.speciescomposition", "bio.taxonomy" ) ) )


  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )


  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )

  p = parameters_add_without_overwriting( p, runlabel="1999_present" )

  p = parameters_add_without_overwriting( p,
    spatial_domain = "SSE",  # canada.east.highres and canada.east.superhighres result in memory overflow
    spatial_domain_subareas = c( "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    dimensionality="space-time" # dimensionality of output data predictions (season is modelled but only a single slice kept for storage issues)
  )

  p$quantile_bounds =c(0.005, 0.995) # trim upper bounds

  p = spatial_parameters( p=p)

  # define focal years for modelling and interpolation

  if (!exists("year.assessment", p )) if (exists("yrs", p)) p$year.assessment=max(p$yrs)

  if (!exists("year.assessment", p )) {
    message("probably want to assign current year.assessment, using current year for now")
    p$year.assessment = lubridate::year(lubridate::now())
  }
  
  yrs_default = 1999:p$year.assessment
  p = parameters_add_without_overwriting( p, yrs=yrs_default, timezone="America/Halifax" )  # default
  p = temporal_parameters(p=p)


  p$discretization = discretizations(p=p$discretization)  # key for discretization levels

  p = parameters_add_without_overwriting( p,
    additional.data=c("groundfish", "snowcrab", "USSurvey_NEFSC", "lobster"),
    taxa =  "maxresolved",
    varstomodel = c( "pca1", "pca2", "pca3", "ca1", "ca2", "ca3" ),
    inputdata_spatial_discretization_planar_km = p$pres/2, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
    inputdata_temporal_discretization_yr = 1/12  # ie., controls resolution of data prior to modelling to reduce data set and speed up modelling;; use 1/12 -- monthly or even 1/4.. if data density is low
  )



    # basic selection criteria
  p = parameters_add_without_overwriting( p,
    selection = list(
      # biologicals=list(
      #   spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=groundfish_survey_species_code )
      # ),
      survey=list(
        data.source = c("groundfish", "snowcrab" ),
        yr = p$yrs,      # time frame for comparison specified above
        # months=6:8,
        # dyear = c(150,250)/365, #  summer = which( (x>150) & (x<250) ) , spring = which(  x<149 ), winter = which(  x>251 )
        # ranged_data="dyear"
        settype = c(1,2,5,8),
        # gear = c("Western IIA trawl", "Yankee #36 otter trawl"),
        # strata_toremove=c("Gulf", "Georges_Bank", "Spring", "Deep_Water"),  # <<<<< strata to remove from standard strata-based analysis
        # polygon_enforce=TRUE
        greater_than = c("sa", "sa_towdistance"),
        sa = 0.001,  # km^2 (snowcrab)
        sa_towdistance = 0.001 
      )
    )
  )

  # NOTE: groundfish settypes:
    # 1=stratified random,
    # 2=regular survey,
    # 3=unrepresentative(net damage),
    # 4=representative sp recorded(but only part of total catch),
    # 5=comparative fishing experiment,
    # 6=tagging,
    # 7=mesh/gear studies,
    # 8=explorartory fishing,
    # 9=hydrography



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

    if (!exists("carstm_model_label", p)) p$carstm_model_label = "1970_present"

    if (exists("carstm_model_label", p)) {

      if (p$carstm_model_label == "1999_present"){
          p$yrs = 1999:p$year.assessment
          p$areal_units_timeperiod = p$carstm_model_label 
      } else if (p$carstm_model_label == "1970_present"){
          p$yrs = 1970:p$year.assessment
          p$areal_units_timeperiod = p$carstm_model_label 
      }
    }
 
 
    p = temporal_parameters(p=p)  # redo in case of user-specificed params
    


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
      fraction_todrop = 0.025,
      fraction_cv = 0.95,
      fraction_good_bad = 1.0,
      areal_units_constraint_ntarget=8, 
      areal_units_constraint_nmin=1,  # granularity options for areal_units
      areal_units_constraint="none",
      # areal_units_constraint_ntarget =  floor(length(p$yrs)/2),  # n time slices req in each au
      # areal_units_constraint_nmin=5,  # best compromise
      spbuffer=5, 
      lenprob=0.95,   # these are domain boundary options for areal_units
      n_iter_drop=0, 
      sa_threshold_km2=16, 
      nAU_min = 50,
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "1999_present",
      carstm_inputs_prefilter = "rawdata"
    )

    if ( !exists("carstm_prediction_surface_parameters", p))  {
        # generics using "default" carstm models and stmv solutions for spatial effects
        p$carstm_prediction_surface_parameters = list()
        p$carstm_prediction_surface_parameters = parameters_add_without_overwriting( p$carstm_prediction_surface_parameters,
          bathymetry = aegis.bathymetry::bathymetry_parameters( project_class="stmv" ),
          substrate = aegis.substrate::substrate_parameters(   project_class="stmv" ),
          temperature = aegis.temperature::temperature_parameters( project_class="carstm", yrs=1999:p$year.assessment, carstm_model_label="1970_present" ) 
        )
    }


    if ( grepl("inla", p$carstm_modelengine) ) {
      if ( !exists("formula", p)  ) {

        p$formula = as.formula( paste(
         p$variabletomodel, ' ~ 1',
#            ' + f( cyclic, model="seasonal", scale.model=TRUE, season.length=10, hyper=H$iid  ) ',  # cannot use seasonal as missing data in some levels
            ' + f( cyclic, model="ar1",  hyper=H$ar1  ) ',
            ' + f( time, model="ar1",  hyper=H$ar1 ) ',
            ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, hyper=H$bym2 ) ',
            ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
            ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
            ' + f( inla.group( log.substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',  # causes issues due to limited spatial range ?
            ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group))'
          ) )
      }

      if ( !exists("family", p)  )  p$family = "gaussian"
    }


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

    using_fft = FALSE
    
    if (exists("stmv_local_modelengine", p) ) if (p$stmv_local_modelengine == "fft") using_fft = TRUE  
    if (exists("stmv_twostep_space", p) ) if (p$stmv_twostep_space == "fft") using_fft = TRUE  

      if (using_fft ) {
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
            + NOTE____UPDATE_MODEL
            + f(space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
          family = "gaussian",
          data= dat,
          inla.mode="compact",
          control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE), return.marginals.predictor=TRUE,  # config=TRUE if doing posterior simulations
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

    return(p)
  }

}
