

speciescomposition_parameters = function( p=NULL, project_name=NULL, project_class="default", ... ) {

  # ---------------------
  # deal with additional passed parameters
 p = parameters_control(p, list(...), control="add") # add passed args to parameter list, priority to args



  # ---------------------
  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools" ) ) )
  p$libs = unique( c( p$libs, project.library ( "aegis", "aegis.speciescomposition" ) ) )

  p$project_name = ifelse ( !is.null(project_name), project_name, "speciescomposition" )

  if ( !exists("data_root", p) ) p$data_root = project.datadirectory( "aegis", p$project_name )
  if ( !exists("datadir", p) )   p$datadir  = file.path( p$data_root, "data" )
  if ( !exists("modeldir", p) )  p$modeldir = file.path( p$data_root, "modelled" )

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=F, recursive=T )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=F, recursive=T )

   if (!exists("variabletomodel", p)) warning( "The dependent variable, p$variabletomodel needs to be defined")
   if (!exists("varstomodel", p) ) p$varstomodel = c( "pca1", "pca2", "ca1", "ca2" )

  if (!exists("spatial_domain", p) ) p$spatial_domain = "SSE"
  if (!exists("spatial_domain_subareas", p)) p$spatial_domain_subareas = c( "snowcrab", "SSE.mpa" )


  p$aegis_dimensionality="space-year"

  p = spatial_parameters( p=p)

  # define focal years for modelling and interpolation
  if (!exists("yrs", p)) p$yrs = c(1999:lubridate::year(lubridate::now()))  # NOTE:: this is short as groundfish species id is inconsistent
  p = temporal_parameters(p=p)

  p$taxa =  "maxresolved"


  if (project_class=="default") {
    if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    if ( !exists("inputdata_temporal_discretization_yr", p)) p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }
    return(p)
  }


  if (project_class=="carstm") {
    if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    if ( !exists("inputdata_temporal_discretization_yr", p)) p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }
    return(p)
  }



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
