
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
      inputdata_spatial_discretization_planar_km = 1,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
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

        M = speciescomposition_db( p=p, DS="carstm_inputs", redo=TRUE  )  # will redo if not found
        # to extract fits and predictions
        M= NULL
        gc()
 
    }

# Time used:
#     Pre = 19.1, Running = 6828, Post = 125, Total = 6972 
# Fixed effects:
#               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
# (Intercept) -0.003 0.053     -0.109   -0.003      0.102 -0.003   0

# Random effects:
#   Name	  Model
#     dyri AR1 model
#    year AR1 model
#    auid_main BYM2 model
#    inla.group(t, method = "quantile", n = 11) RW2 model
#    inla.group(z, method = "quantile", n = 11) RW2 model
#    inla.group(substrate.grainsize, method = "quantile", n = 11) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                                                mean       sd 0.025quant 0.5quant
# Precision for the Gaussian observations                                     325.306    4.961    315.492  325.338
# Precision for dyri                                                          426.868   43.957    354.301  420.280
# Rho for dyri                                                                  0.864    0.017      0.830    0.863
# Precision for year                                                          399.466  214.599     77.074  371.890
# Rho for year                                                                  0.911    0.026      0.848    0.916
# Precision for auid_main                                                    2325.091 1197.850    821.015 2056.728
# Phi for auid_main                                                               NaN      NaN      0.000    0.000
# Precision for inla.group(t, method = "quantile", n = 11)                   4726.692 5206.293    841.337 3159.026
# Precision for inla.group(z, method = "quantile", n = 11)                    318.533  202.355     56.093  278.016
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 1721.444 3425.397    132.476  807.505
# Precision for auid                                                          200.506   11.638    180.860  199.291
# Phi for auid                                                                  0.996    0.004      0.987    0.997
# GroupRho for auid                                                             0.958    0.004      0.950    0.958
#                                                                            0.975quant     mode
# Precision for the Gaussian observations                                      3.35e+02  325.508
# Precision for dyri                                                           5.33e+02  400.220
# Rho for dyri                                                                 8.97e-01    0.862
# Precision for year                                                           8.60e+02  240.942
# Rho for year                                                                 9.49e-01    0.924
# Precision for auid_main                                                      5.38e+03 1629.575
# Phi for auid_main                                                                 NaN      NaN
# Precision for inla.group(t, method = "quantile", n = 11)                     1.82e+04 1767.671
# Precision for inla.group(z, method = "quantile", n = 11)                     8.13e+02  164.351
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 11)   9.01e+03  307.918
# Precision for auid                                                           2.26e+02  195.592
# Phi for auid                                                                 1.00e+00    1.000
# GroupRho for auid                                                            9.67e-01    0.958

# Expected number of effective parameters(stdev): 1716.00(43.99)
# Number of equivalent replicates : 6.46 

# Deviance Information Criterion (DIC) ...............: -30935.74
# Deviance Information Criterion (DIC, saturated) ....: -28333.01
# Effective number of parameters .....................: 1723.80

# Watanabe-Akaike information criterion (WAIC) ...: -30973.15
# Effective number of parameters .................: 1473.10

# Marginal log-Likelihood:  24146.77 
# Posterior marginals for the linear predictor and
#  the fitted values are computed

# Time used:
#     Pre = 19.3, Running = 10776, Post = 167, Total = 10962 
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 0.008 0.019     -0.029    0.008      0.045 0.008   0

# Random effects:
#   Name	  Model
#     dyri AR1 model
#    year AR1 model
#    auid_main BYM2 model
#    inla.group(t, method = "quantile", n = 11) RW2 model
#    inla.group(z, method = "quantile", n = 11) RW2 model
#    inla.group(substrate.grainsize, method = "quantile", n = 11) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                                                mean       sd 0.025quant 0.5quant
# Precision for the Gaussian observations                                    3.49e+02 5.24e+00    338.686 3.49e+02
# Precision for dyri                                                         1.64e+03 3.00e+02   1008.588 1.66e+03
# Rho for dyri                                                               6.57e-01 1.08e-01      0.375 6.85e-01
# Precision for year                                                         4.20e+03 1.80e+03   1593.119 3.92e+03
# Rho for year                                                               7.90e-01 1.19e-01      0.476 8.21e-01
# Precision for auid_main                                                    6.40e+02 1.13e+02    478.240 6.19e+02
# Phi for auid_main                                                               NaN      NaN      0.000 0.00e+00
# Precision for inla.group(t, method = "quantile", n = 11)                   9.37e+03 1.17e+04     45.279 4.55e+03
# Precision for inla.group(z, method = "quantile", n = 11)                   1.21e+02 6.79e+01     35.096 1.06e+02
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 9.13e+05 8.78e+07    339.693 1.59e+04
# Precision for auid                                                         3.67e+02 3.38e+01    317.243 3.61e+02
# Phi for auid                                                               9.85e-01 1.10e-02      0.959 9.88e-01
# GroupRho for auid                                                          9.37e-01 8.00e-03      0.920 9.37e-01
#                                                                            0.975quant     mode
# Precision for the Gaussian observations                                      3.59e+02  348.495
# Precision for dyri                                                           2.21e+03 1782.133
# Rho for dyri                                                                 7.91e-01    0.750
# Precision for year                                                           8.53e+03 3333.746
# Rho for year                                                                 9.28e-01    0.870
# Precision for auid_main                                                      9.12e+02  569.026
# Phi for auid_main                                                                 NaN      NaN
# Precision for inla.group(t, method = "quantile", n = 11)                     4.09e+04    7.285
# Precision for inla.group(z, method = "quantile", n = 11)                     2.94e+02   79.667
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 11)   3.68e+06  466.496
# Precision for auid                                                           4.47e+02  342.795
# Phi for auid                                                                 9.99e-01    0.998
# GroupRho for auid                                                            9.53e-01    0.937

# Expected number of effective parameters(stdev): 1582.64(10.55)
# Number of equivalent replicates : 7.06 

# Deviance Information Criterion (DIC) ...............: -32132.28
# Deviance Information Criterion (DIC, saturated) ....: -32318.73
# Effective number of parameters .....................: 1582.58

# Watanabe-Akaike information criterion (WAIC) ...: -32107.53
# Effective number of parameters .................: 1408.04

# Marginal log-Likelihood:  24822.72 
# Posterior marginals for the linear predictor and
#  the fitted values are computed




      
      # run model and obtain predictions
      fit = carstm_model( p=p, M="speciescomposition_db( p=p, DS='carstm_inputs' ) "  )
      
        # extract results
      if (0) {
        # very large files .. slow 
        fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      }


    res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
    res$summary$dic$dic
    res$summary$dic$p.eff
    res$dyear
   
    if (0) {
      # map all :

      # variabletomodel = "pca1"
      # variabletomodel = "pca2"

        time_match = list(year="2019" )
      
        vn = paste(p$variabletomodel, "predicted", sep=".")
     #   vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
     #   vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")

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
