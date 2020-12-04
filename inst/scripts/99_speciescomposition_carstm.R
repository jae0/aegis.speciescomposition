
# species composition analysis via car

if (!exists("year.assessment")) {
  year.assessment=lubridate::year(Sys.Date())
  year.assessment=lubridate::year(Sys.Date()) -1
}

p = aegis.speciescomposition::speciescomposition_parameters( yrs=1999:year.assessment ) # to get var names

# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
p$varstomodel = c("pca1", "pca2")

for ( variabletomodel in p$varstomodel)  {
    # variabletomodel = "pca1"
    p = aegis.speciescomposition::speciescomposition_parameters(
      project_class="carstm",
      data_root = project.datadirectory( "aegis", "speciescomposition" ),
      variabletomodel = variabletomodel,
      carstm_model_label = "default",
      inputdata_spatial_discretization_planar_km = 1,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = 24/365,  # ie., every 2 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      yrs = 1999:2019,
      aegis_dimensionality="space-year",
      spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_overlay = "none"
    )

    # to recreate the underlying data
    sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
    M = speciescomposition_db( p=p, DS="carstm_inputs", redo=TRUE  )  # will redo if not found
    # to extract fits and predictions

    # run model and obtain predictions
    fit = carstm_model( p=p, M=M )
    fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
    res = carstm_summary( p=p ) # to load currently saved sppoly

    plot(fit)
    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    s = summary(fit)
    s$dic$dic
    s$dic$p.eff

    # maps of some of the results
    carstm_plot( p=p, res=res, vn=paste(p$variabletomodel, "predicted", sep=".") )

    vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
    if (exists(vn, res)) carstm_plot( p=p, res=res, vn=vn, time_match=list(year="1950", dyear="0") )

    vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
    if (exists(vn, res)) {
      res_dim = dim( res[[vn]] )
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.8" )
      carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
    }

    vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
    if (exists(vn, res)) {
      res_dim = dim( res[[vn]] )
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.8" )
      carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
    }

  }

}
# end
