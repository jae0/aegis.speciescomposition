

speciescomposition_parameters = function( p=list(), project_name="speciescomposition", project_class="core", workflow_decentralized=FALSE, ... ) {


  # ---------------------
  # deal with additional passed parameters
  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args


  if (!exists("variabletomodel", p)) warning( "The dependent variable, p$variabletomodel needs to be defined")


  # ---------------------
  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools" ) ) )
  p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.speciescomposition" ) ) )


  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( "aegis", p$project_name ) )
  p = parameters_add_without_overwriting( p, datadir  = file.path( p$data_root, "data" ) )
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )

  # for projects that require access to default data and local data, a switch is needed to force use of default data
  if ( p$workflow_decentralized )  {
    if (exists( "modeldir_override", p)) p$modeldir = p$modeldir_override  # must also specify p$workflow_decentralized =TRUE  for override to work
  }

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )


  p = parameters_add_without_overwriting( p,
    spatial_domain = "SSE",  # canada.east.highres and canada.east.superhighres result in memory overflow
    spatial_domain_subareas = c( "SSE", "SSE.mpa" , "snowcrab"),  # this is for bathymetry_db, not stmv
    aegis_dimensionality="space-year"
  )

  p = spatial_parameters( p=p)

  # define focal years for modelling and interpolation
  p = parameters_add_without_overwriting( p, yrs = 1950:lubridate::year(lubridate::now()) )  # default
  p = temporal_parameters(p=p)


  p = parameters_add_without_overwriting( p,
    additional.data=c("groundfish", "snowcrab", "USSurvey_NEFSC", "lobster"),
    taxa =  "maxresolved",
    varstomodel = c( "pca1", "pca2", "ca1", "ca2" ),
    inputdata_spatial_discretization_planar_km = 1, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
    inputdata_temporal_discretization_yr = 1/12  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling;; use 1/12 -- monthly or even 1/4.. if data density is low
  )



  # ---------------------

  if (project_class=="core")  return(p)

  # ---------------------

  if (project_class=="carstm") {
    # simple run of carstm. There are two types:
    #   one global, run directly from  polygons defined in aegis.bathymetry/inst/scripts/99.bathymetry.carstm.R
    #   and one that is called secondarily specific to a local project's polygons (eg. snow crab)
    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )


    # defaults in case not provided ...
    p = parameters_add_without_overwriting( p,
      areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_resolution_km = 5, # default in case not provided ... 25 km dim of lattice ~ 1 hr; 5km = 79hrs; 2km = ?? hrs
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
      areal_units_overlay = "none",
      carstm_modelengine = "inla",  # {model engine}.{label to use to store}
      carstm_model_label = "default",
      carstm_inputs_aggregated = FALSE
    )


    if ( !exists("carstm_model_call", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$carstm_model_call = paste(
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
            control.fixed = H$fixed,  # priors for fixed effects, generic is ok
            control.inla = list( cmin = 0, h=1e-4, tolerance=1e-9, strategy="adaptive", optimise.strategy="smart"), # restart=3), # restart a few times in case posteriors are poorly defined
            verbose=TRUE
          )'
        )
      }
        #    + f(tiyr, model="ar1", hyper=H$ar1 )
        # + f(year,  model="ar1", hyper=H$ar1 )
    }

    p = carstm_parameters( p=p )  #generics

    return(p)
  }


  # ---------------------

  if (project_class %in% c("stmv")) {

    p = parameters_add_without_overwriting( p,
      project_class="stmv",
      stmv_model_label="default",
      stmv_variables = list(
        LOCS=c("plon", "plat"),
        TIME="tiyr"
      ),  # required as fft has no formulae
      inputdata_spatial_discretization_planar_km = p$pres / 4, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
      inputdata_temporal_discretization_yr = 1/12,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      stmv_global_modelengine = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_local_modelengine="fft",
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_runmode = list(
        carstm = rep("localhost", 1),
        globalmodel = FALSE,
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )  # ncpus for each runmode

      p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )


     p = aegis_parameters(p=p, DS="stmv" ) # generics for aegis.* projects

    return(p)
  }


  # ---------------------

  if (project_class %in% c("hybrid", "default")) {

    p = parameters_add_without_overwriting( p,
      project_class="stmv",
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
          formula = z ~ 1
            + f(auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
          family = "normal",
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
      stmv_distance_basis_interpolation = c(  2.5 , 5, 10, 15, 20, 40, 80, 150, 200 ) , # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
      stmv_nmin = 90, # min number of data points req before attempting to model in a localized space
      stmv_nmax = 1000, # no real upper bound.. just speed /RAM
      stmv_runmode = list(
        carstm = rep("localhost", 1),
        globalmodel = FALSE,
        save_intermediate_results = TRUE,
        save_completed_data = TRUE
      )  # ncpus for each runmode

      p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )


    p = aegis_parameters(p=p, DS="stmv" ) # generics for aegis.* projects

    return(p)
  }

}
