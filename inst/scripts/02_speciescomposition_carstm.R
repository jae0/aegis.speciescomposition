
# species composition analysis via car

if (!exists("year.assessment")) {
  year.assessment=lubridate::year(Sys.Date())
  year.assessment=lubridate::year(Sys.Date()) -1
}

p = aegis.speciescomposition::speciescomposition_parameters( yrs=1999:year.assessment )

# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
for ( variabletomodel in p$varstomodel)  {

    p = aegis.speciescomposition::speciescomposition_parameters(
      project_class = "carstm", # defines which parameter set to load
      data_root = project.datadirectory( "aegis", "speciescomposition" ),
      variabletomodel = variabletomodel,
      inputdata_spatial_discretization_planar_km = 1,  # km controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = 24/365,  # ie., every 2 weeks .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      yrs = 1999:2019,
      aegis_dimensionality="space-year",
      spatial_domain = "SSE",  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
    )

    if (0) {
      # to recreate the underlying data
      sppoly = areal_units( p=p, redo=TRUE )  # this has already been done in aegis.polygons::01 polygons.R .. should nto have to redo
      M = speciescomposition.db( p=p, DS="aggregated_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
      M = speciescomposition_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
      # to extract fits and predictions

      # run model and obtain predictions
      res = speciescomposition_carstm( p=p, DS="carstm_modelled", redo=TRUE )

      res = speciescomposition_carstm( p=p, DS="carstm_modelled" ) # to load currently saved res
      fit =  speciescomposition_carstm( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      s = summary(fit)
      s$dic$dic
      s$dic$p.eff

      # maps of some of the results
      carstm_plot( p=p, res=res, vn="speciescomposition.predicted" )
      carstm_plot( p=p, res=res, vn="speciescomposition.random_strata_nonspatial" )
      carstm_plot( p=p, res=res, vn="speciescomposition.random_strata_spatial" )

    }

}
# end



