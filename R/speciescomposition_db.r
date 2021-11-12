
  speciescomposition_db = function( DS="", p=NULL, redo=FALSE ) {

    ddir = project.datadirectory( "aegis", "speciescomposition" )
    dir.create( ddir, showWarnings=FALSE, recursive=TRUE )

    infix = paste( p$spatial_domain,  p$taxa, sep=".")
    if (p$spatial_domain == "snowcrab" ) {
      infix = paste( "SSE",  p$taxa, sep=".")  # just one domain for now
    }

    if (DS %in% c( "speciescomposition.ordination", "speciescomposition.ordination.redo", "pca", "ca") ) {

      fn.set = file.path( ddir, paste( "speciescomposition.by.set", infix, "rdata", sep=".") )
      fn.pca = file.path( ddir, paste( "pca", infix, "rdata", sep=".") )
      fn.ca  = file.path( ddir, paste( "ca",  infix, "rdata", sep=".") )

      if (DS=="speciescomposition.ordination") {
        set = NULL
        if (file.exists( fn.set) ) load( fn.set)
        return ( set )
      }

      if (DS=="pca") {
        pca.out = NULL
        if (file.exists( fn.pca) ) load( fn.pca)
        return ( pca.out )
      }

      if (DS=="ca") {
        ca.out = NULL
        if (file.exists( fn.ca) ) load( fn.ca)
        return ( ca.out )
      }


 


    # # convert from quantile to z-score
    # set$zm = quantile_to_normal( set$qm )
    # set$zn = quantile_to_normal( set$qn )

      require(data.table)

      # catch info
      sc = survey_db( DS="cat" ,p=p)[ , c("id", "totno", "totwgt", "cf_cat", "spec_bio" )]
      sc = setDT( sc, key="id" ) # species catch
 
      # trip/set loc information
      set = survey_db( DS="set", p=p)

      set = set[ which( set$settype %in% c(1,2,4,5,8) ) , ] # == "good"
        # settype:
        # 1=stratified random,
        # 2=regular survey,
        # 3=unrepresentative(net damage),
        # 4=representative sp recorded(but only part of total catch),
        # 5=comparative fishing experiment,
        # 6=tagging,
        # 7=mesh/gear studies,
        # 8=explorartory fishing,
        # 9=hydrography

    
      # filter area
      igood = which( set$lon >= p$corners$lon[1] & set$lon <= p$corners$lon[2]
              &  set$lat >= p$corners$lat[1] & set$lat <= p$corners$lat[2] )
      set = set[igood, ]

      isc = taxonomy.filter.taxa( sc$spec_bio, method=p$taxa, outtype="internalcodes" )
      set = set[ which( set$id %in% unique( sc$id[isc]) ),]

      set = set[ , c("id", "yr", "dyear", "sa", "lon", "lat", "t", "z", "timestamp", "gear", "vessel", "data.source" )]
      i = which(!is.finite(set$z)) 
      if (length(i) > 0 ) {
        set$z[i] = aegis_lookup(  
            parameters= "bathymetry" , 
            LOCS=set[ i, c("lon", "lat")],  
            project_class="core", 
            output_format="points" , 
            DS="aggregated_data", 
            variable_name="z.mean",
            returntype="vector" ,
            yrs=p$yrs
          ) 
      }

      i = which(!is.finite(set$t)) 
      if (length(i) > 0 ) {
        set$t[i] = aegis_lookup(  
            parameters= "temperature" , 
            LOCS=set[ i, c("lon", "lat", "timestamp")],  
            project_class="carstm", 
            output_format="points" , 
            variable_name=list( "predictions"),
            statvars=c("mean"),
            raster_resolution=min(p$gridparams$res)/2,
            returntype="vector" ,
            yrs=p$yrs
          ) 
      }

      sc = merge(sc, set, by="id", all.x=TRUE, all.y=FALSE) 
      sc = na.omit(sc)

      sc$totno_adjusted = trunc( sc$totno * sc$cf_cat )   # correct to catch subsmapling 
      sc$log_sa = log(sc$sa) 
      
      if (0) {

          # way too slow to use a modelled solution:
          # gc()

          # o = as.data.frame(sc) 
          # o$tag ="obs"
          # o$spec_bio = as.numeric( as.factor( o$spec_bio))

          # u=o
          # u$tag ="preds"
          # u$log_sa = 0
          # u$totno_adjusted = NA
          # u$vessel = "CAR"
          # u$gear = "US 4 seam 3 bridle survey trawl"
          # u$data.source = "groundfish"

          # o =  rbind( o, u[ , names(o)])
          # u = NULL; gc()
          # o$survey = paste( o$data.source, o$vessel, o$gear, sep="__")
          # o$yr = o$yr - min(o$yr) + 1
          # o$yr1 = o$yr2 = o$yr3 = o$yr

          # nmod = inla( 
          #     formula = totno_adjusted ~ 1 + offset(log_sa) 
          #       + f( yr, model="ar1")
          #       + f( id, model="iid", group=yr3 )
          #       + f( spec_bio, model="iid", group=yr1 )
          #       + f( survey, model="iid", group=yr2),  
          #     data = o,
          #     family = "poisson",
          #     control.predictor=list(link=1),
          #     control.inla = list( strategy='adaptive', int.strategy="eb" ),
          #     # inla.mode="experimental",
          #     verbose =TRUE
          # )

          # # .. data loss due to truncation is OK
          # # ... smallest abundance adds little information to ordinations
          

          # sc$density =  nmod$summary.fitted.values$mean[ which(o$tag=="preds" ) ]
          # sc$zdensity = NA 

          # sps = unique( sc$spec_bio )
          # for ( s in sps ) {
          #   ii = which( sc$spec_bio==s & sc$density > 0 )
          #   if (length( ii) > 5 ) {
          #     sc$zdensity[ii] = quantile_to_normal( quantile_estimate( sc$density[ii] ) )
          #   }
          # }

      }
 
      
      # the US 4 seam behaves very differently and has a shorter standard tow .. treat as a separate data source 
      # of course if misses out on the high abundance period from the 1980s but no model-based solution possible at this point 
      # due simply to CPU speed issues

      sc$survey = sc$data.source
      i = which( sc$gear=="US 4 seam 3 bridle survey trawl" )
      sc$survey[ i ]  = paste( sc$survey[i], sc$gear[i], sep="__")
      sc$qn = NA 

      surveys = unique(sc$survey)
      for ( s in surveys ) {
        ii = which( sc$survey==s & sc$totno_adjusted > 0 )
        if (length( ii) > 0 ) sc$qn[ii] = quantile_estimate( sc$totno_adjusted[ii] )  # convert to quantiles, by survey
      }

      oo = which( sc$totno_adjusted == 0 )  # retain as zero values
      if (length(oo)>0 ) sc$qn[oo] = 0


      # convert from quantile to z-score
      sc$zn = quantile_to_normal( sc$qn )


      m = xtabs( zn*1e6 ~ as.factor(id) + as.factor(spec_bio), data=sc )  / 1e6
      
      # remove low counts (absence) in the timeseries  .. species (cols) only
      cthreshold = 0.1   # 
      
      finished = F
      while( !(finished) ) {
        i = unique( which(rowSums(m) == 0 ) )
        j = unique( which(colSums(m) <= cthreshold ) )
        if ( ( length(i) == 0 ) & (length(j) == 0 ) ) finished=T
        if (length(i) > 0 ) m = m[ -i , ]
        if (length(j) > 0 ) m = m[ , -j ]
      }

      # m = log(m)

      # PCA
      # no need to correct for gear types/surveys .. assuming no size-specific bias .. perhaps wrong but simpler
      corel = cor( m, use="pairwise.complete.obs" ) # set up a correlation matrix ignoring NAs
      corel[ is.na(corel) ] = 0
      s = svd(corel)  # eigenanalysis via singular value decomposition

      # scores = matrix.multiply (m, s$v)  # i.e., b %*% s$v  .. force a multiplication ignoring NAs
      m[which(!is.finite(m))] = 0
      scores = m %*% s$v  # i.e., b %*% s$v  .. force a multiplication ignoring NAs

      evec = s$v
      ev = s$d
      x = cbind( scores[,1] / sqrt(ev[1] ), scores[,2] / sqrt( ev[2]), scores[,3] / sqrt(ev[3] ) )
      y = cbind( evec[,1] * sqrt(ev[1] ) , evec[,2] * sqrt( ev[2]),  evec[,3] * sqrt( ev[3]) )
      rownames(y) = colnames(m)

      scores = data.frame( id=rownames(m), pca1=as.numeric(x[,1]), pca2=as.numeric(x[,2]), pca3=as.numeric(x[,3]), stringsAsFactors=FALSE )
      set = merge(set, scores, by="id", all.x=T, all.y=F, sort=FALSE)
      pca.out = list( scores=scores, eignenvectors=evec, eigenvalues=ev, cscores=y )
      save( pca.out, file=fn.pca, compress=TRUE)


      # Correpsondence analysis
      require(vegan)
      n = m * 0
      n[ which(m>0) ] = 1
      ord = cca( n )
      sp.sc = scores(ord, choices=c(1:3))$species
      si.sc = scores(ord, choices=c(1:3))$sites
      scores = data.frame( id=as.character(rownames(si.sc)), ca1=as.numeric(si.sc[,1]), ca2=as.numeric(si.sc[,2]), ca3=as.numeric(si.sc[,3]) )
      variances=  ord$CA$eig[1:10]/sum(ord$CA$eig)*100
      set = merge(set, scores, by="id", all.x=T, all.y=F, sort=F)
      ca.out = list( scores=scores, ca=ord, variances=variances )

      save( ca.out, file=fn.ca, compress=T)
      save( set, file=fn.set, compress=T )

      return (fn.set)
    }



    # -----------------------



    if (DS %in% c( "speciescomposition", "speciescomposition.redo" ) ) {
      # remove dups and in planar coords
      fn = file.path( ddir, paste( "speciescomposition", infix, "rdata", sep=".") )

			if (DS=="speciescomposition") {
        SC = NULL
        if (file.exists( fn) ) load( fn )
        return ( SC )
			}

      SC = speciescomposition_db( DS="speciescomposition.ordination", p=p )
      SC = lonlat2planar( SC, proj.type=p$aegis_proj4string_planar_km )
 
      oo = which(!is.finite( SC$plon+SC$plat ) )
      if (length(oo)>0) SC = SC[ -oo , ]  # a required field for spatial interpolation

      yrs = sort( unique( SC$yr ) )
      # check for duplicates
      for ( y in yrs ) {
        yy = which (SC$yr == y)
        ii = which( duplicated( SC$id[yy] ) )

        if (length(ii) > 0) {
          dd = which( SC$id[yy]  %in% SC$id[yy] [which(duplicated(SC$id[yy] ))])
          print( "The following sets have duplicated positions. The first only will be retained" )
          print( SC[yy,] [dd ] )
          SC = SC[ - ii,]
        }
      }

      imperative = c( "pca1", "pca2",  "pca3", "ca1", "ca2", "ca3" )
      ii = which( is.finite( rowSums( SC[,imperative ] ) ) )
      if (length(ii) == 0) stop( "No data .. something went wrong")
      SC = SC[ii,]

      save( SC, file=fn, compress=T )
			return (fn)
		}


    # -------------

    if ( DS=="areal_units_input" ) {

      fn = file.path( p$datadir,  "areal_units_input.rdata" )
      if ( !file.exists(p$datadir)) dir.create( p$datadir, recursive=TRUE, showWarnings=FALSE )

      xydata = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( xydata )
        }
      }
      xydata = speciescomposition_db( p=p, DS="speciescomposition"  )  #

      names(xydata)[which(names(xydata)=="z.mean" )] = "z"
      xydata = xydata[ geo_subset( spatial_domain=p$spatial_domain, Z=xydata ) , ] # need to be careful with extrapolation ...  filter depths
      xydata = xydata[ , c("lon", "lat", "yr" )]
      save(xydata, file=fn, compress=TRUE )
      return( xydata )
    }


    # -------------


    if ( DS=="carstm_inputs") {

      # prediction surface
      crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
      sppoly = areal_units( p=p )  # will redo if not found
      sppoly = st_transform(sppoly, crs=crs_lonlat )
      areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
      
      # over-ride
      p$variabletomodel = "speciescomposition"  # force to be generic to share across variables 

      label =  "carstm_inputs" 
      if (p$carstm_inputs_prefilter =="rawdata") label = "carstm_inputs_rawdata"
      fn = carstm_filenames( p=p, returntype=label, areal_units_fn=areal_units_fn )

      # inputs are shared across various secneario using the same polys
      #.. store at the modeldir level as default
      outputdir = dirname( fn )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      M = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          M = M[ which( M$yr %in% p$yrs), ]
          return( M )
        }
      }
      # message( "Generating carstm_inputs ... ")


      # do this immediately to reduce storage for sppoly (before adding other variables)

      M = speciescomposition_db( p=p, DS="speciescomposition"  )
      setDT(M)
      
      vars_to_retain = c("pca1", "pca2", "pca3", "ca1", "ca2", "ca3" )

      if (p$carstm_inputs_prefilter != "aggregated") {
        if (exists("quantile_bounds", p)) {
          for (vn in vars_to_retain) {
            TR = quantile(M[[vn]], probs=p$quantile_bounds, na.rm=TRUE )
            oo = which( M[[vn]] < TR[1])
            if (length(oo) > 0) M[[vn]][oo] = TR[1]
            oo = which( M[[vn]] > TR[2])
            if (length(oo) > 0) M[[vn]][oo] = TR[2]
          }
        }
      }

      # INLA does not like duplicates ... causes optimizer to crash frequently
      # oo = which(duplicated( M[, p$variabletomodel] ))
      # if ( length(oo)> 0) {
      #   eps = exp( log( .Machine$double.eps ) / 2)  # ~ 1.5e-8
      #   M[oo, p$variabletomodel]  = M[oo, p$variabletomodel]  + runif( length(oo), -eps, eps )
      # }

      M = planar2lonlat(M, proj.type=p$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
      
      ii = which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] )
      M = M[ ii, ]

      names(M)[which(names(M)=="yr") ] = "year"
      M = M[ which(M$year %in% p$yrs), ]
      M$tiyr = lubridate::decimal_date ( M$timestamp )
      M$dyear = M$tiyr - M$year

 
      M = carstm_prepare_inputdata( p=p, M=M, sppoly=sppoly, 
            lookup_parameters = p$carstm_lookup_parameters,
            vars_to_retain=vars_to_retain,
            vars_to_drop ="speciescomposition" )  # drop dummy variable
 
      jj = which(!is.finite(M$t))
      if (length(jj) > 0 ) {
        M$t[jj] = median( M$t[-jj] )
      }

      jj = which(!is.finite(M$substrate.grainsize))
      if (length(jj) > 0 ) {
        M$substrate.grainsize[jj] = median( M$substrate.grainsize[-jj] )
      }

      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

      save( M, file=fn, compress=TRUE )

      M = M[ which( M$yr %in% p$yrs), ]

      return( M )
    }



  } # end function



