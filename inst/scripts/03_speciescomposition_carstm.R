
# species composition analysis via car

year.assessment = 2021

require( aegis.speciescomposition )


for ( variabletomodel in c("pca1", "pca2", "pca3"))  {
    
    # variabletomodel = "pca1"
    # variabletomodel = "pca2"
    # variabletomodel = "pca3"
    
    # construct basic parameter list defining the main characteristics of the study
    p = speciescomposition_parameters(
      project_class="carstm",
      data_root = project.datadirectory( "aegis", "speciescomposition" ),
      variabletomodel = variabletomodel,
      carstm_model_label = "default",
      inputdata_spatial_discretization_planar_km = 0.5,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = 1/52,  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      year.assessment = year.assessment,
      yrs = 1999:year.assessment,
      aegis_dimensionality="space-year",
      spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
 #     areal_units_type = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
       areal_units_type = "tesselation", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not     
       areal_units_overlay = "none"
    )
  

    if (0) { 
        # p$fraction_todrop = 1/11 # aggressiveness of solution finding ( fraction of counts to drop each iteration)
        # p$fraction_cv = 1.0 #sd/mean no.
        # p$fraction_good_bad = 0.9
        # p$areal_units_constraint_nmin =  3
        # p$areal_units_constraint_ntarget = 15  # length(p$yrs)

        # p$nAU_min = 100

        # # adjust based upon RAM requirements and ncores
        # require(INLA)
        # inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
        # inla.setOption(blas.num.threads= 2 )

        # to recreate the underlying data
        xydata = speciescomposition_db(p=p, DS="areal_units_input", redo=TRUE)

        sppoly = areal_units( p=p, redo=TRUE, verbose=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      
        plot(sppoly["AUID"])

    }


    M = speciescomposition_db( p=p, DS="carstm_inputs", redo=TRUE  )  # will redo if not found .. .
    # to extract fits and predictions
    M= NULL
    gc()

    
    # run model and obtain predictions
    fit = carstm_model( 
      p=p, 
      data="speciescomposition_db( p=p, DS='carstm_inputs' ) ", 
      num.threads="4:2",
      # control.inla = list( strategy='adaptive', int.strategy='eb' ),
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
        tmap_zoom= c(map_centre, map_zoom),
        title=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") ) ,
        outfilename=fn
      )

    }

  }



# end
