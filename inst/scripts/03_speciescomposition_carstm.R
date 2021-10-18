
# species composition analysis via car

year.assessment = 2021

require( aegis.speciescomposition )

# construct basic parameter list defining the main characteristics of the study
# choose one:

# default run: 1970:present
p0 = speciescomposition_parameters(
  project_class="carstm",
  data_root = project.datadirectory( "aegis", "speciescomposition" ),
  variabletomodel = "",  # will b eover-ridden .. this brings in all pca's and ca's
  carstm_model_label = "default",
  inputdata_spatial_discretization_planar_km = 0.5,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  inputdata_temporal_discretization_yr = 1/52,  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
  year.assessment = year.assessment,
  yrs = 1970:year.assessment,
  aegis_dimensionality="space-year",
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
#     areal_units_type = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
    areal_units_type = "tesselation", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not     
    areal_units_overlay = "none"
)


# for bio.snowcrab 1999:present
p0 = speciescomposition_parameters(
  project_class="carstm",
  data_root = project.datadirectory( "aegis", "speciescomposition" ),
  variabletomodel = "",  # will b eover-ridden .. this brings in all pca's and ca's
  carstm_model_label = "1999_present",
  inputdata_spatial_discretization_planar_km = 0.5,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  inputdata_temporal_discretization_yr = 1/52,  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
  year.assessment = year.assessment,
  yrs = 1970:year.assessment,
  aegis_dimensionality="space-year",
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
#     areal_units_type = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
  areal_units_type = "tesselation", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not     
  areal_units_overlay = "none",
  carstm_lookup_parameters = list( 
    bathymetry = bathymetry_parameters( project_class="stmv", spatial_domain="SSE", stmv_model_label="default"  ),
    substrate = substrate_parameters(   project_class="stmv", spatial_domain="SSE", stmv_model_label="default"  ),
    temperature = temperature_parameters( project_class="carstm", carstm_model_label="1999_present", yrs=1999:year.assessment ) 
  )
)


if (0) { 
    # p0$fraction_todrop = 1/11 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
    # p0$fraction_cv = 0.7 #sd/mean no.
    # p0$fraction_good_bad = 0.9
    # p0$areal_units_constraint_nmin =  3
    # p0$areal_units_constraint_ntarget = 15  # length(p0$yrs)

    # p$nAU_min = 100

    # # adjust based upon RAM requirements and ncores
    # require(INLA)
    # inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
    # inla.setOption(blas.num.threads= 2 )

    # to recreate the underlying data
    xydata = speciescomposition_db(p=p0, DS="areal_units_input", redo=TRUE)

    sppoly = areal_units( p=p0, redo=TRUE, hull_alpha=20, verbose=TRUE )   
  
      plot(sppoly["AUID"])

}

 
 
  


M = speciescomposition_db( p=p0, DS="carstm_inputs", redo=TRUE  )  # will redo if not found .. .
str(M); 
M= NULL; gc()

p0$formula = NULL  # reset to force a new default below 

for ( variabletomodel in c("pca1", "pca2" )) { #  } , "ca1", "ca2",  "pca3", "ca3"))  {
    
    # variabletomodel = "pca1"
    # variabletomodel = "pca2"
    # variabletomodel = "pca3"
    
    # construct basic parameter list defining the main characteristics of the study
    p = speciescomposition_parameters( p=p0, project_class="carstm", variabletomodel = variabletomodel, mc.cores=2 )  # mc.cores == no of cores to use for posterior extraction .. can be memory intensive so keep low 

    # run model and obtain predictions
    fit = carstm_model( 
      p=p, 
      data="speciescomposition_db( p=p, DS='carstm_inputs' ) ", 
      num.threads="6:2",  # adjust for your machine
      # control.inla = list( strategy='laplace' ), # "adaptive" strategy seems to run into problems with sparse data (in current year) 
      control.inla = list( strategy='adaptive' ),
      control.mode = list(
        theta = switch( variabletomodel,
          pca1 = c(5.779, 4.330, 15.261, -1.951, 4.617, 3.950, 5.652, 3.552, 3.611 ),
          pca2 = c( 6.014, 5.303, 11.065, 20.323, 9.516, 4.379, 5.920, 3.823, 3.528 ),
          pca3 = c( 6.014, 5.303, 11.065, 20.323, 9.516, 4.379, 5.920, 3.823, 3.528 ),
          ca1 = c( 2.658, 1.229, 10.311, 19.854, 8.145, 0.250, 2.298, 5.280, 3.342 ),
          ca2 = c( 6.014, 5.303, 11.065, 20.323, 9.516, 4.379, 5.920, 3.823, 3.528 ),
          ca3 = c( 6.014, 5.303, 11.065, 20.323, 9.516, 4.379, 5.920, 3.823, 3.528 )
        ), 
        restart=TRUE
      ),  
      
      redo_fit=F, # to start optim from a solution close to the final in 2021 ... 
      verbose=TRUE 
     )
    
      # extract results
    if (0) {
      # very large files .. slow 
      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }

    fit = NULL; gc()

    res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results

    map_centre = c( (p$lon0+p$lon1)/2 - 0.5, (p$lat0+p$lat1)/2 -0.8 )
    map_zoom = 6.5


    if (0) {
      # map all :
      vn=c( "random", "space", "combined" )
      vn=c( "random", "spacetime", "combined" )
      vn="predictions"
      tmatch="2015"

      carstm_map(  res=res, vn=vn, tmatch=tmatch, 
        plot_crs = "+proj=omerc +lat_0=44.5 +lonc=-63.5 +gamma=0.0 +k=1 +alpha=332 +x_0=0 +y_0=0 +ellps=WGS84 +units=km" ,
        palette="RdYlBu",
        breaks = seq(-0.3, 0.3, by=0.1),
        plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
        tmap_zoom= c(map_centre, map_zoom),
        title=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") )  
      )

    }


    outputdir = file.path( gsub( ".rdata", "", carstm_filenames(p, returntype="carstm_modelled_fit") ), "figures" )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    graphics.off()

    for (y in res$time ){
      tmatch = as.character(y) 
      fn_root = paste( "speciescomposition", variabletomodel, paste0(tmatch, collapse=" - "), sep="_" )
      fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

      vn="predictions"
    
      carstm_map(  res=res, vn=vn, tmatch=tmatch , 
        palette="RdYlBu",
        breaks = seq(-0.3, 0.3, by=0.1),
        plot_elements=c( "isobaths", "coastline", "compass", "scale_bar", "legend" ),
        map_mode="plot",
        tmap_zoom= c(map_centre, map_zoom),
        title=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") ) ,
        outfilename=fn
      )

    }

  }



# end
