 
# -----------------------------
# ordination of all years 1970 to present
year.assessment = 2023

yrs = 1999:year.assessment
 
runlabel="1999_present"
require(aegis)
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
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  areal_units_type = "tesselation",     
  areal_units_constraint="none",
  #areal_units_resolution_km = 1, # km dim of lattice ~ 1 hr
  # areal_units_overlay = "none",
  # spbuffer=5, lenprob=0.95,   # these are domain boundary options for areal_units
  # n_iter_drop=0, sa_threshold_km2=4, 
  # areal_units_constraint_ntarget=10, areal_units_constraint_nmin=1,  # granularity options for areal_units
  carstm_prediction_surface_parameters = list( 
    bathymetry = aegis.bathymetry::bathymetry_parameters( project_class="stmv" ),
    substrate = aegis.substrate::substrate_parameters(   project_class="stmv" ),
    temperature = aegis.temperature::temperature_parameters( project_class="carstm", spatial_domain="canada.east", yrs=1999:year.assessment, carstm_model_label="1999_present" ) 
  ), 
  theta = list(     
    pca1 = c( 6.596, 9.250, 0.013, 8.114, 2.547, 12.028, 0.039, 13.987, 6.422, 12.998, 6.650, 5.292, 3.444  ),   
    pca2 = c( 6.613, 6.758, 2.117, 7.785, 2.909, 13.631, -1.938, 9.406, 5.407, 10.801, 6.426, 6.341, 3.391 ),
    pca3 = c( 7.001, 8.069, 0.073, 3.958, -6.423, 7.670, 4.010, 11.813, 7.701, 12.266, 7.635, 3.808, 2.474  )#, 
    # ca1 =  c(  6.512, 5.082, 6.983, 3.188, 9.191, 2.793, 9.261, 7.840, 12.532, 6.479, 4.650, 3.541   ),
    # ca2 =  c(  6.512, 5.082, 6.983, 3.188, 9.191, 2.793, 9.261, 7.840, 12.532, 6.479, 4.650, 3.541   ),
    # ca3 =  c(  6.512, 5.082, 6.983, 3.188, 9.191, 2.793, 9.261, 7.840, 12.532, 6.479, 4.650, 3.541   )
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
    sppoly = areal_units( p=p0, xydata=xydata, redo=TRUE, verbose=TRUE )  # to force create
 
    plot(sppoly["npts"])

}
 
# do this once for the default (all years) ... the shorter are subset from it 
# .. if not then this needs to be rerun for the full set if years do not span chosen subset
sppoly = areal_units( p=p0)
M = speciescomposition_db( p=p0, DS="carstm_inputs", sppoly=sppoly , redo=TRUE  )  # will redo if not found .. .
str(M); 
M= NULL; gc()

p0$space_name = sppoly$AUID 
p0$space_id = 1:nrow(sppoly)  # must match M$space

p0$time_name = as.character(p0$yrs)
p0$time_id =  1:p0$ny

p0$cyclic_name = as.character(p0$cyclic_levels)
p0$cyclic_id = 1:p0$nw


for ( variabletomodel in c( "pca1", "pca2", "pca3")) { #  , "pca3" , "ca1", "ca2",   "ca3"))  {
    
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
      # required
      runlabel=runlabel
    )  
    
    # run model and obtain predictions
    res = carstm_model( 
      p=p, 
      data="speciescomposition_db( p=p, DS='carstm_inputs' ) ", 
      nposteriors=5000,
      posterior_simulations_to_retain=c(  "random_spatial", "predictions"), 
      theta=p$theta[[variabletomodel]],
      # redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
      num.threads="6:2",  # adjust for your machine
      # debug = TRUE,
      # control.inla = list( strategy='adaptive', int.strategy='eb' ),  # "eb" required for stabilization
      # control.inla = list( strategy='laplace'),
      # control.inla=list(cmin=0),
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
      
      fit = NULL; gc()
    }

}


# bbox = c(-71.5, 41, -52.5,  50.5 )
additional_features = features_to_add( 
    p=p0, 
    isobaths=c( 100, 200, 300, 400, 500  ), 
    xlim=c(-80,-40), 
    ylim=c(38, 60) , redo=TRUE
)


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
      runlabel=runlabel 
    )  
 
    res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results

    dir.create( file.path(p$data_root, "figures"), recursive=TRUE, showWarnings=FALSE)
    fnp = file.path(p$data_root, "figures", paste(variabletomodel, "_timeseries.png", sep="") )
    png(filename=fnp, width=800,height=600, res=144)
    carstm_plotxy( res, vn=c( "res", "random", "time" ), reverse=TRUE, # reverse is to match maps colors .. they are also reversed
      type="b", xlab="Year", ylab=variabletomodel, xv=p$yrs, ylim=c(-0.05, 0.05)  )  # override xv (x-values)
    dev.off()

 
    fnp = file.path(p$data_root, "figures", paste(variabletomodel, "_cyclic.png", sep="") )
    png(filename=fnp, width=800,height=600, res=144)
    carstm_plotxy( res, vn=c( "res", "random", "cyclic" ), reverse=TRUE,
      type="b", col="slategray", pch=19, lty=1, lwd=2.5, ylim=c(-0.065, 0.065),
      xlab="Season", ylab=variabletomodel, h=0.5  )
    dev.off()

    fnp = file.path(p$data_root, "figures", paste(variabletomodel, "_temperature.png", sep="") )
    png(filename=fnp, width=800,height=600, res=144)
    carstm_plotxy( res, vn=c( "res", "random", "inla.group(t, method = \"quantile\", n = 9)" ), reverse=TRUE,
      type="b", col="slategray", pch=19, lty=1, lwd=2.5  ,
      xlab="Bottom temperature (degrees Celsius)", ylab=variabletomodel   )
    dev.off()


    fnp = file.path(p$data_root, "figures", paste(variabletomodel, "_depth.png", sep="") )
    png(filename=fnp, width=800,height=600, res=144)
    carstm_plotxy( res, vn=c( "res", "random", "inla.group(z, method = \"quantile\", n = 9)" ), reverse=TRUE,
      type="b", col="slategray", pch=19, lty=1, lwd=2.5  ,
      xlab="Depth (m)", ylab=variabletomodel   )
    dev.off()


    map_centre = c( (p$lon0+p$lon1)/2  , (p$lat0+p$lat1)/2   )
    map_zoom = 7.4


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
        colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
        breaks = seq(-0.2, 0.1 , by=0.05),
        additional_features=additional_features,
        legend.position=c( 0.1, 0.9 ),
        annotation=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") )  
      )

    }

    
    outputdir = file.path(p$data_root, "maps", p$carstm_model_label )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    vn="predictions"

    qn = quantile(  carstm_results_unpack( res, vn )[,,"mean"], probs=c(0.1, 0.9), na.rm=TRUE  )
    brks = pretty( qn ) 

    # graphics.off()
 
    for (y in res$time_id ){
      tmatch = as.character(y) 
      fn_root = paste( "speciescomposition", variabletomodel, paste0(tmatch, collapse=" - "), sep="_" )
      outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )
    
      plt = carstm_map(  res=res, vn=vn, tmatch=tmatch , 
        breaks = brks,
        annotation=paste("Species composition: ", variabletomodel, "  ", paste0(tmatch, collapse="-") ), 
        legend.position=c( 0.1, 0.9 ),
        colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),  #RdYlBu BuYlRd
        additional_features=additional_features,
        outfilename=outfilename
      )
    }

  
    # pure spatial effect
    vn=c( "random", "space", "combined" )
   
    fn_root = paste( "speciescomposition", variabletomodel, "spatial_effect", sep="_" )
    outfilename = file.path( outputdir, paste(fn_root, "png", sep=".") )

    toplot = carstm_results_unpack( res, vn )
    brks = pretty(  quantile(toplot[,"mean"], probs=c(0.025, 0.975), na.rm=TRUE )  )

  ##NOTE: color scale is reversed for red=hot
    plt = carstm_map(  res=res, vn=vn, 
        sppoly = sppoly, 
        colors= (RColorBrewer::brewer.pal(5, "RdYlBu")),
        breaks = brks,
        annotation=paste("Species composition: ", variabletomodel, "persistent spatial effect" ), 
        legend.position=c( 0.1, 0.9 ),
        additional_features=additional_features,
        outfilename=outfilename
    )
 

    # posterior predictive check
    
    MM = speciescomposition_db( p=p, DS='carstm_inputs' )
    iobs = which(MM$tag == "observations")
    vn = variabletomodel

    fit = NULL; gc()
    fit = carstm_model( p=p, DS="carstm_modelled_fit") #,  sppoly = sppoly )

    pld = data.table(
      observed = MM[iobs , ..vn] , 
      fitted = fit$summary.fitted.values[["mean"]] [iobs]
    )
    names(pld) = c("observed", "fitted")
    anno1 = paste( "Pearson correlation: ", round( cor( pld$fitted, pld$observed, use="pairwise.complete.obs" ), 3))
    # cor( fitted, observed, use="pairwise.complete.obs", "spearman" )

    out = ggplot(pld, aes(x =  observed, y = fitted )) +
      geom_abline(slope=1, intercept=0, color="darkgray", lwd=1.4 ) +
      geom_point(color="slategray") +
      labs(caption=anno1, color="slateblue") +
      theme( plot.caption = element_text(hjust = 0, size=12 ) )# move caption to the left 
   
    outputdir = file.path( p$modeldir, p$carstm_model_label )
    fn = file.path(outputdir, paste("posterior_predictive_check_", vn, ".png", sep="") )
    ggsave(filename=fn, plot=out, device="png", width=12, height = 8)
    print(out)    

    fit  = MM = NULL; gc()


  }


# end
