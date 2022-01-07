

  # -----------------------------
  # ordination of all years 1970 to present
  yrs = 1970:2021
  runlabel="1970_present"

  require(aegis.speciescomposition)
  


  # -----------------------------
  # prepare data

  p = speciescomposition_parameters( yrs=yrs, runlabel=runlabel )

  speciescomposition_db( DS="speciescomposition.ordination.redo", p=p )  # analsysis

  speciescomposition_db( DS="speciescomposition.redo", p=p ) # compute planar coords and remove dups

  if (0) { 
    # extract summaries and plot
    pca = speciescomposition_db( DS="pca", p=p )  # analsysis
    ca  = speciescomposition_db( DS="ca", p=p )  # analsysis

    toplot = as.data.frame( pca$loadings )
    toplot$vern = taxonomy.recode( from="spec", to="taxa", tolookup=rownames( toplot ) )$vern

    plot( PC2 ~ PC1, toplot, type="n")
    text( PC2 ~ PC1, labels=vern, data=toplot )

    plot( PC3 ~ PC1, toplot, type="n")
    text( PC3 ~ PC1, labels=vern, data=toplot )

  }


  

# -----------------------------
# carstm predictions / analysis

# parameter template 
p0 = speciescomposition_parameters(
  project_class="carstm",
  data_root = project.datadirectory( "aegis", "speciescomposition" ),
  variabletomodel = "",  # will b eover-ridden .. this brings in all pca's and ca's
  runlabel = runlabel,
  carstm_model_label = runlabel,
  inputdata_spatial_discretization_planar_km = 0.5,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  inputdata_temporal_discretization_yr = 1/52,  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
  year.assessment =  max(yrs),
  yrs = yrs,
  aegis_dimensionality="space-year",
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  areal_units_type = "tesselation",    
  areal_units_overlay = "none",
  carstm_lookup_parameters = list( 
      bathymetry = aegis.bathymetry::bathymetry_parameters( project_class="stmv" ),
      substrate = aegis.substrate::substrate_parameters(   project_class="stmv" ),
      temperature = aegis.temperature::temperature_parameters( project_class="carstm", yrs=yrs, carstm_model_label=runlabel ) 
  )
  ,
  theta0 =  list(     
    pca1 = c( 5.633, 7.664, 4.900, 3.740, 7.041, 2.208, 9.606, 3.152, 5.141, 2.335, 3.868 ),  ## z
    pca2 = c( 6.207, 8.640, 6.834, 4.381, 8.431, 1.576, 9.208, 0.871, 5.781, 1.463, 3.919 ),  ## t
    pca3 = c( 5.640, 5.510, 5.372, 3.478, 6.495, 2.323, 6.246, 3.696, 5.105, 2.094, 3.855 ),   
    ca1 = c( 5.640, 5.510, 5.372, 3.478, 6.495, 2.323, 6.246, 3.696, 5.105, 2.094, 3.855 ), 
    ca2 = c( 5.640, 5.510, 5.372, 3.478, 6.495, 2.323, 6.246, 3.696, 5.105, 2.094, 3.855 ),
    ca3 = c( 5.640, 5.510, 5.372, 3.478, 6.495, 2.323, 6.246, 3.696, 5.105, 2.094, 3.855 )
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
 
# do this once for the default (1970:present) ... the shorter are subset from it 
# .. if not then this needs to be rerun for the full set if years do not span chosen subset
M = speciescomposition_db( p=p0, DS="carstm_inputs", sppoly=areal_units( p=p0), redo=TRUE  )  # will redo if not found .. .
str(M); 
M= NULL; gc()


for ( variabletomodel in c("pca1", "pca2")) { #  , "pca3" , "ca1", "ca2",   "ca3"))  {
    
    # variabletomodel = "pca1"
    # variabletomodel = "pca2"
    # variabletomodel = "pca3"
    
    # construct basic parameter list defining the main characteristics of the study
    p0$formula = NULL  # MUST reset to force a new formulae to be created on the fly below 
    p = speciescomposition_parameters( p=p0, project_class="carstm", variabletomodel = variabletomodel, yrs=p0$yrs, runlabel=runlabel,
      mc.cores=2, theta=p0$theta0[[variabletomodel]]
    )  
    
    # run model and obtain predictions
    
    fit = carstm_model( 
      p=p, 
      data="speciescomposition_db( p=p, DS='carstm_inputs' ) ", 
      num.threads="6:2",  # adjust for your machine
      control.inla = list( strategy='adaptive', int.strategy='eb' ),  # "eb" required for stabilization
      redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
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

    map_centre = c( (p$lon0+p$lon1)/2  , (p$lat0+p$lat1)/2   )
    map_zoom = 7
    background = tmap::tm_basemap(leaflet::providers$CartoDB.Positron, alpha=0.8 )

    if (0) {
      # map all :
      vn=c( "random", "space", "combined" )
      vn=c( "random", "spacetime", "combined" )
      vn="predictions"
      tmatch="2015"

      qn = quantile(  carstm_results_unpack( res, vn )[,,"mean"], probs=c(0.1, 0.9), na.rm=TRUE  )
      brks = pretty( qn ) 

      carstm_map(  res=res, vn=vn, tmatch=tmatch, 
        plot_crs = "+proj=omerc +lat_0=44.5 +lonc=-63.5 +gamma=0.0 +k=1 +alpha=332 +x_0=0 +y_0=0 +ellps=WGS84 +units=km" ,
        palette="RdYlBu",
        breaks = seq(-0.2, 0.1 , by=0.05),
        plot_elements=c( "isobaths", "compass", "scale_bar", "legend" ),
        tmap_zoom= c(map_centre, map_zoom),
        background = background,
        title=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") )  
      )

    }


    outputdir = file.path( gsub( ".rdata", "", carstm_filenames(p, returntype="carstm_modelled_fit") ), "figures" )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    vn="predictions"
    
    qn = quantile(  carstm_results_unpack( res, vn )[,,"mean"], probs=c(0.1, 0.9), na.rm=TRUE  )
    brks = pretty( qn ) 

    graphics.off()

    for (y in res$time ){
      tmatch = as.character(y) 
      fn_root = paste( "speciescomposition", variabletomodel, paste0(tmatch, collapse=" - "), sep="_" )
      fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

    
      carstm_map(  res=res, vn=vn, tmatch=tmatch , 
        palette="RdYlBu",
        breaks = seq(-0.3, 0.3, by=0.1),
        plot_elements=c( "isobaths", "compass", "scale_bar", "legend" ),
        map_mode="view",
        tmap_zoom= c(map_centre, map_zoom),
        background=background, 
        title=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") ) ,
        outfilename=fn
      )

    }

    # pure spatial effect
    vn=c( "random", "space", "combined" )
   
    fn_root = paste( "speciescomposition", variabletomodel, "spatial_effect", sep="_" )
    fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

    carstm_map(  res=res, vn=vn, tmatch=tmatch , 
        palette="RdYlBu",
        breaks = seq(-0.3, 0.3, by=0.1),
        plot_elements=c( "isobaths", "compass", "scale_bar", "legend" ),
        map_mode="view",
        tmap_zoom= c(map_centre, map_zoom),
        background=background, 
        title=paste("Species composition: ", variabletomodel, "  ", "spatial_effect" ) ,
        outfilename=fn
    )
  }



# end
