  # -----------------------------
  # ordination
  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date())
    year.assessment=lubridate::year(Sys.Date()) -1
  }

  p = aegis.speciescomposition::speciescomposition_parameters( yrs=1999:year.assessment )


  # vn="pca1"  # ~ 8-10 hrs / variable
  # vn="pca2"
  ncpus = parallel::detectCores()

  if (0) {
    # about 700 MB/process ..
    ram_required_per_process = 1  # about 600 MB per process in 2017 GB
    ncpus = min( parallel::detectCores(), floor( ram_local() / ram_required_per_process ) )
    ncpus.covars = min( parallel::detectCores(), floor( ram_local() / .7 ) )  # 700 GB in 2018 .. for prediction of global covars
  }


  for ( vn in p$varstomodel) {
    print(vn)
    p = aegis.speciescomposition::speciescomposition_parameters(
      project_class="stmv",
      data_root = project.datadirectory( "aegis", "speciescomposition" ),
      DATA = 'aegis_db( p=p, DS="stmv_inputs" )',
      stmv_variables=list(Y=vn),
      yrs = c(1999:year.assessment),  # years for modelling and interpolation
      spatial_domain = "SSE",
      aegis_dimensionality="space-year",
      stmv_global_modelengine = "gam",
      stmv_global_family = gaussian(link="identity"),
      stmv_global_modelformula = formula( paste(
        vn,
        ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( tmin, k=3, bs="ts") + s( tmax, k=3, bs="ts") + s( degreedays, k=3, bs="ts")   ',
#        ' + s( t, tsd, tmin, tmax, degreedays, k=50, bs="ts")  ',
        ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
#        ' + s( log(z), log(dZ), log(ddZ), k=30, bs="ts")  ',
        ' + s( substrate.grainsize, k=3, bs="ts") ' )),  # no space
      stmv_local_modelengine ="twostep",
      stmv_twostep_time = "gam",
      stmv_twostep_space = "fft",  # other possibilities: "fft", "tps"
      stmv_fft_filter="matern_tapered",  #  matern, krige (very slow), lowpass, lowpass_matern
      # stmv_lowpass_nu = 0.1,
      # stmv_lowpass_phi = stmv::matern_distance2phi( distance=0.25, nu=0.1, cor=0.5 ), # default p$res = 0.5;
      stmv_autocorrelation_fft_taper = 0.5,  # benchmark from which to taper
      stmv_autocorrelation_localrange=0.1,
      stmv_autocorrelation_basis_interpolation = c(0.5, 0.1, 0.05, 0.01),
      stmv_variogram_method = "fft",
      stmv_filter_depth_m = 0, # the depth covariate is input as log(depth) so, choose stats locations with elevation > log(1 m) as being on land
      stmv_local_model_distanceweighted = TRUE,
      stmv_rsquared_threshold = 0.2, # lower threshold
      stmv_distance_statsgrid = 4, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_distance_scale = c(20, 30, 40, 50), # km ... approx guess of 95% AC range .. data tends to be sprse realtive to pure space models
      stmv_nmin = 8*(year.assessment-1999),# floor( 7 * p$ny ) # min number of data points req before attempting to model timeseries in a localized space
      stmv_nmax = 8*(year.assessment-1999)*11, # max( floor( 7 * p$ny ) * 11, 8000), # no real upper bound
      stmv_runmode = list(
        globalmodel = TRUE,
        scale = rep("localhost", scale_ncpus),
        interpolate = list(
            cor_0.5 = rep("localhost", interpolate_ncpus),
            cor_0.1 = rep("localhost", interpolate_ncpus),
            cor_0.05 = rep("localhost", max(1, interpolate_ncpus-1)),
            cor_0.01 = rep("localhost", max(1, interpolate_ncpus-2))
          ),  # ncpus for each runmode
        interpolate_predictions = list(
          c1 = rep("localhost", max(1, interpolate_ncpus-1)),  # ncpus for each runmode
          c2 = rep("localhost", max(1, interpolate_ncpus-1)),  # ncpus for each runmode
          c3 = rep("localhost", max(1, interpolate_ncpus-2)),
          c4 = rep("localhost", max(1, interpolate_ncpus-3)),
          c5 = rep("localhost", max(1, interpolate_ncpus-4)),
          c6 = rep("localhost", max(1, interpolate_ncpus-4)),
          c7 = rep("localhost", max(1, interpolate_ncpus-5))
        ),
        save_intermediate_results = FALSE,
        save_completed_data = TRUE # just a dummy variable with the correct name
      )  # ncpus for each runmode
    )

    stmv( p=p )

    if (0) {
    # quick view
      predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
      statistics  = stmv_db( p=p, DS="stmv.stats" )
      locations   = spatial_grid( p )

      # comparisons
      dev.new(); surface( as.image( Z=rowMeans(predictions), x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

      # stats
      (p$statsvars)
      dev.new(); levelplot( predictions[,1] ~ locations[,1] + locations[,2], aspect="iso" )
      dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
      dev.new(); levelplot( statistics[,match("sdTot", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
      dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange
    }

    aegis_db( p=p, DS="predictions.redo" ) # warp predictions to other grids
    aegis_db( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
    aegis_db( p=p, DS="complete.redo" )
    aegis_db( p=p, DS="baseline.redo" )
    aegis_db_map( p=p )

    gc()

    if (0) {
      global_model = stmv_global_model( p=p, DS="global_model")
      summary( global_model )
      plot(global_model)
    }

  }




# model testing
if (0) {

  p$stmv_global_modelformula = formula( paste(
    'ca1',
    ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( tmin, k=3, bs="ts") + s( tmax, k=3, bs="ts") + s( degreedays, k=3, bs="ts")   ',
    ' + s( t, tsd, tmin, tmax, degreedays, k=50, bs="ts")  ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    ' + s( log(z), log(dZ), log(ddZ), k=30, bs="ts")  ',
    ' + s( log(substrate.grainsize), k=3, bs="ts") '
  ))

  p$stmv_global_modelformula = formula( paste(
    'ca2',
    ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( tmin, k=3, bs="ts") + s( tmax, k=3, bs="ts") + s( degreedays, k=3, bs="ts")   ',
    ' + s( t, tsd, tmin, tmax, degreedays, k=50, bs="ts")  ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    ' + s( log(z), log(dZ), log(ddZ), k=30, bs="ts")  ',
    ' + s( log(substrate.grainsize), k=3, bs="ts") '
  ))

  require(mgcv)
  o = aegis.db( p=p, DS="stmv_inputs" )  # create fields for

  global_model = bam(
    formula=p$stmv_global_modelformula,
    family=p$stmv_global_family,
    data = o,
    weights=o$wt,
    method="fREML",
    use.chol=TRUE,
    gc.level=2,
    na.action="na.omit"
  )

  global_model = gam(
    formula=p$stmv_global_modelformula,
    family=p$stmv_global_family,
    data = o,
    weights=o$wt,
    optimizer= p$stmv_gam_optimizer,
    na.action="na.omit"
  )

  summary( global_model )
  plot(global_model, all.terms=TRUE, seWithMean=TRUE, jit=TRUE, rug=TRUE )
  # plot(global_model, all.terms=TRUE, trans=bio.snowcrab::inverse.logit, seWithMean=TRUE, jit=TRUE, rug=TRUE )

}


# ----------------

Family: gaussian
Link function: identity

Formula:
pca1 ~ s(t, k=3, bs="ts") + s(tsd, k=3, bs="ts") + s(tmin,
    k=3, bs="ts") + s(tmax, k=3, bs="ts") + s(degreedays,
    k=3, bs="ts") + s(log(z), k=3, bs="ts") + s(log(dZ),
    k=3, bs="ts") + s(log(ddZ), k=3, bs="ts") + s(log(substrate.grainsize),
    k=3, bs="ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.108107   0.000859     126   <2e-16

Approximate significance of smooth terms:
                             edf Ref.df      F p-value
s(t)                        2.00      2 209.65 < 2e-16
s(tsd)                      1.95      2  29.42 7.7e-14
s(tmin)                     0.70      2   5.94 < 2e-16
s(tmax)                     1.96      2  34.33 < 2e-16
s(degreedays)               1.27      2  24.84 < 2e-16
s(log(z))                   1.63      2 205.26 < 2e-16
s(log(dZ))                  1.44      2  17.65 8.6e-10
s(log(ddZ))                 1.90      2   8.39 0.00013
s(log(substrate.grainsize)) 1.98      2 125.18 < 2e-16

R-sq.(adj) =  0.686   Deviance explained = 68.6%
GCV = 0.011143  Scale est. = 0.011131  n = 15083



# -------------------

Family: gaussian
Link function: identity

Formula:
pca2 ~ s(t, k=3, bs="ts") + s(tsd, k=3, bs="ts") + s(tmin,
    k=3, bs="ts") + s(tmax, k=3, bs="ts") + s(degreedays,
    k=3, bs="ts") + s(log(z), k=3, bs="ts") + s(log(dZ),
    k=3, bs="ts") + s(log(ddZ), k=3, bs="ts") + s(log(substrate.grainsize),
    k=3, bs="ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.043390   0.000749    57.9   <2e-16

Approximate significance of smooth terms:
                                 edf Ref.df       F p-value
s(t)                        1.94e+00      2  214.35 < 2e-16
s(tsd)                      1.76e+00      2    6.29 0.00089
s(tmin)                     1.99e+00      2 1055.06 < 2e-16
s(tmax)                     1.95e+00      2  189.32 < 2e-16
s(degreedays)               7.96e-08      2    0.00 0.93740
s(log(z))                   1.73e+00      2 5336.95 < 2e-16
s(log(dZ))                  1.84e+00      2   54.81 < 2e-16
s(log(ddZ))                 8.83e-08      2    0.00 0.87388
s(log(substrate.grainsize)) 1.98e+00      2  214.75 < 2e-16

R-sq.(adj) =  0.685   Deviance explained = 68.6%
GCV = 0.0084705  Scale est. = 0.0084625  n = 15083

