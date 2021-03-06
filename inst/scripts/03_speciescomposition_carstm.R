
# species composition analysis via car

year.assessment = 2020

require( aegis.speciescomposition )


for ( variabletomodel in c("pca1", "pca2"))  {
    
    # variabletomodel = "pca1"
    # variabletomodel = "pca2"
    
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
    res = carstm_model( p=p, M="speciescomposition_db( p=p, DS='carstm_inputs' ) "  )
    
      # extract results
    if (0) {
      # very large files .. slow 
      fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      fit$summary$dic$dic
      fit$summary$dic$p.eff
      fit$dyear
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    }


    res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
  
    if (0) {
      # map all :

        time_match = list(year="2019" )
      
        vn = paste(p$variabletomodel, "predicted", sep=".")
    #   vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
    #   vn = paste(p$variabletomodel, "random_space", sep=".")

        carstm_map(  res=res, vn=vn, time_match=time_match , 
          plot_crs = "+proj=omerc +lat_0=44.5 +lonc=-63.5 +gamma=0.0 +k=1 +alpha=332 +x_0=0 +y_0=0 +ellps=WGS84 +units=km" ,
          main=paste("Species composition: ", variabletomodel, "  ", paste0(time_match, collapse="-") )  
        )

    }

    plot_crs = p$aegis_proj4string_planar_km

    coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
    isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400 ), project_to=plot_crs  )


    vn = paste( variabletomodel, "predicted", sep=".")
    outputdir = file.path( gsub( ".rdata", "", dirname(res$fn_res) ), "figures", vn )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    for (y in res$year ){
      time_match = list( year=as.character(y)  )
      fn_root = paste( "speciescomposition", variabletomodel, paste0(time_match, collapse=" - "), sep="_" )
      fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

        carstm_map(  res=res, vn=vn, time_match=time_match , 
          coastline=coastline,
          isobaths=isobaths,
          palette="RdYlBu",
          breaks = seq(-0.3, 0.3, by=0.1),
          main=paste("Species composition: ", variabletomodel, "  ", paste0(time_match, collapse="-") ) ,
          outfilename=fn
        )

    }

  }



# end
