

speciescomposition_parameters = function( p=list(), project_name="speciescomposition", project_class="default", reset_data_location=FALSE, ... ) {


  if (reset_data_location) {
    # reset a few project specific params, forcing the use of defaults (below)
    p$data_root = NULL
    p$datadir  = NULL
    p$carstm_modelcall = NULL  # defaults to generic
    p$carstm_model_tag = NULL
    p$variabletomodel = NULL
    p$aegis_dimensionality = NULL
    p$data_transformation = NULL
  }

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

  if (project_class=="default")  return(p)

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


    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
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

  if (project_class=="stmv") {
    p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )
    if (!exists("stmv_variables", p)) p$stmv_variables = list()
    if (!exists("LOCS", p$stmv_variables)) p$stmv_variables$LOCS=c("plon", "plat")
    if (!exists("TIME", p$stmv_variables)) p$stmv_variables$TIME="tiyr"

    p = aegis_parameters(p=p, DS="stmv" ) # generics:
    p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }

    return(p)
  }

}
