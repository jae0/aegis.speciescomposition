 
# -----------------------------
# ordination of all years 1970 to present
year.assessment = 2021

yrs = 1999:year.assessment
runlabel="1999_present"

require(aegis.speciescomposition)

p = speciescomposition_parameters( yrs=yrs, runlabel=runlabel )


# -----------------------------
# prepare data
speciescomposition_db( DS="speciescomposition.ordination.redo", p=p )  # analsysis

speciescomposition_db( DS="speciescomposition.redo", p=p  ) # compute planar coords and remove dups

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

p0 = speciescomposition_parameters(
  project_class="carstm",
  data_root = project.datadirectory( "aegis", "speciescomposition" ),
  variabletomodel = "",  # will b eover-ridden .. this brings in all pca's and ca's
  runlabel = runlabel,
  carstm_model_label = runlabel,
  inputdata_spatial_discretization_planar_km = 0.5,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
  inputdata_temporal_discretization_yr = 1/52,  # ie., every 1 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
  year.assessment = max(yrs),
  yrs = yrs,
  aegis_dimensionality="space-year",
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
  areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  areal_units_type = "tesselation",       
  areal_units_overlay = "none",
  carstm_prediction_surface_parameters = list( 
    bathymetry = aegis.bathymetry::bathymetry_parameters( project_class="stmv" ),
    substrate = aegis.substrate::substrate_parameters(   project_class="stmv" ),
    temperature = aegis.temperature::temperature_parameters( project_class="carstm", spatial_domain="canada.east", yrs=1970:year.assessment, carstm_model_label="1970_present" ) 
  ) 
  ,
  theta = list(   
    pca1 = c(  6.502, 6.274, 6.920, 3.247, 9.126, 2.932, 12.019, 7.808, 6.398, 4.631, 3.499 ),  # good
    pca2 = c(  6.603, 5.534, 7.152, 3.170, 9.120, 2.728, 9.643, 6.218, 6.381, 4.885, 3.365   ), 
    pca3 = c(  6.804, 5.387, 4.078, 6.557, 9.347, 2.788, 9.617, 7.854, 6.552, 4.746, 3.539   ),
    ca1 =  c(  6.512, 5.082, 6.983, 3.188, 9.191, 2.793, 9.261, 7.840, 6.479, 4.650, 3.541   ),
    ca2 =  c(  6.512, 5.082, 6.983, 3.188, 9.191, 2.793, 9.261, 7.840, 6.479, 4.650, 3.541   ),
    ca3 =  c(  6.512, 5.082, 6.983, 3.188, 9.191, 2.793, 9.261, 7.840, 6.479, 4.650, 3.541   )
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

    xydata = speciescomposition_db(p=p0, DS="areal_units_input" )
    xydata = xydata[ which(xydata$yr %in% p$yrs), ]
    sppoly = areal_units( p=p0, xydata=xydata, hull_alpha=20, redo=TRUE )  # to force create
 
    plot(sppoly["npts"])

}
 
# do this once for the default (all years) ... the shorter are subset from it 
# .. if not then this needs to be rerun for the full set if years do not span chosen subset
sppoly = areal_units( p=p0)
M = speciescomposition_db( p=p0, DS="carstm_inputs", sppoly=sppoly , redo=TRUE  )  # will redo if not found .. .
str(M); 
M= NULL; gc()


for ( variabletomodel in c("pca1", "pca2", "pca3")) { #  , "pca3" , "ca1", "ca2",   "ca3"))  {
    
    # variabletomodel = "pca1"
    # variabletomodel = "pca2"
    # variabletomodel = "pca3"
    
    # construct basic parameter list defining the main characteristics of the study
    p0$formula = NULL  # MUST reset to force a new formulae to be created on the fly below 
    p = speciescomposition_parameters( 
      p=p0, 
      project_class="carstm", 
      variabletomodel = variabletomodel, 
      yrs=p0$yrs, 
      runlabel=runlabel,
      mc.cores=2, 
      theta=p0$theta[[variabletomodel]]
    )  
    
    # run model and obtain predictions
    fit = carstm_model( 
      p=p, 
      data="speciescomposition_db( p=p, DS='carstm_inputs' ) ", 
      num.threads="6:2",  # adjust for your machine
      # control.inla = list( strategy='laplace'),
      # control.inla = list( strategy='adaptive', int.strategy='eb' ),  # "eb" required for stabilization
      redo_fit=TRUE, # to start optim from a solution close to the final in 2021 ... 
      # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
      # debug = TRUE,
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



    carstm_plotxy( res, vn=c( "res", "random", "time" ), 
      type="b", ylim=c(-0.1, 0.1), xlab="Year", ylab=variabletomodel  )

    carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), 
      type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.025, 0.025),
      xlab="Season", ylab=variabletomodel, h=0.5  )

    carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 9)" ), 
      type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.015, 0.01) ,
      xlab="Bottom temperature (degrees Celcius)", ylab=variabletomodel   )

    carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 9)" ), 
      type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.02, 0.02) ,
      xlab="Depth (m)", ylab=variabletomodel   )
 

    map_centre = c( (p$lon0+p$lon1)/2  , (p$lat0+p$lat1)/2   )
    map_zoom = 7.4
    background = tmap::tm_basemap(leaflet::providers$CartoDB.Positron, alpha=0.8 )

    if (0) {
      # map all :
      vn=c( "random", "space", "combined" )
      vn=c( "random", "spacetime", "combined" )
      vn="predictions"
      tmatch="1999"

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

    outputdir = file.path(p$data_root, "maps", p$carstm_model_label )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


    vn="predictions"

    qn = quantile(  carstm_results_unpack( res, vn )[,,"mean"], probs=c(0.1, 0.9), na.rm=TRUE  )
    brks = pretty( qn ) 

    # graphics.off()

    # SLOW: FASTER TO RUN year of interest and then take a screenshot

    for (y in res$time ){
      tmatch = as.character(y) 
      fn_root = paste( "speciescomposition", variabletomodel, paste0(tmatch, collapse=" - "), sep="_" )
      outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
    
      tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch , 
        palette="RdYlBu",
        breaks = seq(-0.3, 0.3, by=0.1),
        plot_elements=c( "isobaths", "compass", "scale_bar", "legend" ),
        map_mode="view",
        tmap_zoom= c(map_centre, map_zoom),
        background=background, 
        title=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") ) 
      )
      mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 
      print(outfilename)
    }

  
    # pure spatial effect
    vn=c( "random", "space", "combined" )
   
    fn_root = paste( "speciescomposition", variabletomodel, "spatial_effect", sep="_" )
    outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

    tmout = carstm_map(  res=res, vn=vn, tmatch=tmatch , 
        palette="RdYlBu",
        breaks = seq(-0.3, 0.3, by=0.1),
        plot_elements=c( "isobaths", "compass", "scale_bar", "legend" ),
        map_mode="view",
        tmap_zoom= c(map_centre, map_zoom),
        background=background, 
        title=paste("Species composition: ", variabletomodel, "  ", "spatial_effect" ) 
    )
    print(outfilename)
    mapview::mapshot( tmap_leaflet(tmout), file=outfilename, vwidth = 1600, vheight = 1200 )  # very slow: consider 



  }


# end
