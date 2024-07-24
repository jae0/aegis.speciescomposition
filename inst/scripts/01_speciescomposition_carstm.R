 
# -----------------------------
# ordination of all years 1970 to present
year.start = 1999
year.assessment = 2023

yrs = year.start:year.assessment
 
carstm_model_label="default"
require(aegis)
require(aegis.speciescomposition)
require(vegan)

p = speciescomposition_parameters( yrs=yrs, carstm_model_label=carstm_model_label )


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
  variabletomodel = "",  # will be over-ridden .. this brings in all pca's and ca's
  carstm_model_label = carstm_model_label,
  year.assessment = max(yrs),
  yrs = yrs, 
  spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
  theta = list(     
    pca1 = c(3.6110,6.1471,0.0875,5.1708,2.3268,10.0566,-0.6076,11.2549,3.5392,9.9422,3.5000,3.7251,3.6416 ),   
    pca2 = c( 4.3355, 4.7796, 1.2823, 6.0580, 2.5184, 9.9576, -0.5588, 7.5028, 3.5334, 9.8473, 3.9602, 5.6906, 3.5997 ),
    pca3 = c( 3.6170, 6.0820, 0.0889, 5.3414, 2.3956, 10.1444, -0.5991, 11.2819, 3.5499, 10.0855, 3.5549, 3.4580, 3.5575 ), 
    ca1 =  c(  3.6170, 6.0820, 0.0889, 5.3414, 2.3956, 10.1444, -0.5991, 11.2819, 3.5499, 10.0855, 3.5549, 3.4580, 3.5575),
    ca2 =  c( 3.6170, 6.0820, 0.0889, 5.3414, 2.3956, 10.1444, -0.5991, 11.2819, 3.5499, 10.0855, 3.5549, 3.4580, 3.5575 ),
    ca3 =  c( 3.6170, 6.0820, 0.0889, 5.3414, 2.3956, 10.1444, -0.5991, 11.2819, 3.5499, 10.0855, 3.5549, 3.4580, 3.5575 )
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
 
    sppoly = areal_units( p=p0, xydata=xydata, redo=TRUE, verbose=TRUE )  # to force create
    dim(sppoly)
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


for ( variabletomodel in c( "pca1", "pca2" )) { # "pca1", "pca2", "pca3" , "ca1", "ca2",   "ca3"))  {
    
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
      carstm_model_label=carstm_model_label
    )  
    
    # run model and obtain predictions
    res = carstm_model( 
      p=p, 
      data="speciescomposition_db( p=p, DS='carstm_inputs' , sppoly=sppoly ) ", 
      # nposteriors=5000,
      sppoly=sppoly,
      toget = c("summary", "random_spatial", "predictions"),
      # posterior_simulations_to_retain = c("predictions"),  # not used at the moment
      theta=p$theta[[variabletomodel]],
      #redo_fit=FALSE, # to start optim from a solution close to the final in 2021 ... 
      num.threads="4:2",  # adjust for your machine
      # debug = TRUE,
      control.inla = list( strategy='adaptive', int.strategy='eb' ),  # "eb" required for stabilization
      # control.inla = list( strategy='laplace'),
      # control.inla=list(cmin=0),
      verbose=TRUE 
    )
 
}


# bbox = c(-71.5, 41, -52.5,  50.5 )
additional_features = features_to_add( 
    p=p0, 
    isobaths=c( 100, 200, 300, 400, 500  ), 
    xlim=c(-80,-40), 
    ylim=c(38, 60) , redo=TRUE
)


for ( variabletomodel in c("pca1", "pca2" )) { #  , "pca3" , "ca1", "ca2",   "ca3"))  {
    
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
      carstm_model_label=carstm_model_label 
    )  
 
    #  extract results
  
    # very large files .. slow 
    # fit = carstm_model( p=p, DS="modelled_fit" )  # extract currently saved model fit
  
    # plot(fit)
    # plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    
    # posterior predictive check
    M = speciescomposition_db( p=p, DS='carstm_inputs' )
    carstm_posterior_predictive_check(p=p, M=M  )
 
    # EXAMINE POSTERIORS AND PRIORS
    res = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and direct summary
 
    names(res$hypers)
    for (i in 1:length(names(res$hypers)) ){
      o = carstm_prior_posterior_compare( hypers=res$hypers, all.hypers=res$all.hypers, vn=names(res$hypers)[i] )  
      dev.new(); print(o)
    }
 
    # oeffdir = file.path(p$data_root, "figures")  # old ... delete files ..todo
    oeffdir = file.path(p$data_root, p$carstm_model_label, "figures") #new ...
    fn_root_prefix = variabletomodel
    carstm_plot_marginaleffects(  p=p, outputdir=oeffdir, fn_root_prefix=fn_root_prefix ) 

    # maps 
    ## >>>>>>>>>>>> reconsider file save locations / names:  <<<<<<<<<<<<<<
    outputdir = file.path(p$data_root, p$carstm_model_label, "maps" )
  
    carstm_plot_map( p=p, outputdir=outputdir, fn_root_prefix=variabletomodel,
      toplot="random_spatial", probs=c(0.025, 0.975),    
      additional_features=additional_features, 
      colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) ) 

    carstm_plot_map( p=p, outputdir=outputdir, fn_root_prefix=variabletomodel,
      toplot="predictions", 
      additional_features=additional_features, 
      colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")) )


  }


# end
