
  speciescomposition_db = function( DS="", p=NULL, redo=FALSE, species=NULL, cm=NULL, rotate="none", nfactors=2, cthreshold = 0.005, ev_template=NULL, sppoly=NULL ) {

    ddir = project.datadirectory( "aegis", "speciescomposition" )
    dir.create( ddir, showWarnings=FALSE, recursive=TRUE )

    infix = paste( p$spatial_domain,  p$taxa, p$carstm_model_label, sep=".")
    if (p$spatial_domain == "snowcrab" ) {
      infix = paste( "SSE",  p$taxa, sep=".")  # just one domain for now
    }


    if (DS %in% c( "speciescomposition.ordination", "speciescomposition.ordination.redo", "pca", "ca") ) {

      fn.set = file.path( ddir, paste( "speciescomposition.by.set", infix, "rdz", sep=".") )
      fn.pca = file.path( ddir, paste( "pca", infix, "rdz", sep=".") )
      fn.ca  = file.path( ddir, paste( "ca",  infix, "rdz", sep=".") )

      if (DS=="speciescomposition.ordination") {
        set = NULL
        if (file.exists( fn.set) ) set = read_write_fast( fn.set)
        return ( set )
      }

      if (DS=="pca") {
        pca.out = NULL
        if (file.exists( fn.pca) ) pca.out = read_write_fast( fn.pca)
        return ( pca.out )
      }

      if (DS=="ca") {
        ca.out = NULL
        if (file.exists( fn.ca) ) ca.out = read_write_fast( fn.ca)
        return ( ca.out )
      }

 
      require(data.table)
 
      res = survey_data_prepare(p=p, cthreshold = cthreshold)  # drop species that are rare (frequency of occurance in surveys < 0.005)
      
      # unpack
      m = res$m 
      set = res$set
      res = NULL

      # PCA
      # no need to correct for gear types/surveys .. assuming no size-specific bias .. perhaps wrong but simpler
      message("TODO: currently using simple Pearson, using a model-based AC-adjusted cor is better")
      cm = cor( ifelse(m > cthreshold, 1, NA) * m , use="pairwise.complete.obs" ) # set up a correlation matrix ignoring NAs (and low incidence)
 
      cm[ is.na(cm) ] = 0  # reset to 0
      m2 = (m - mean(m, na.rm=TRUE))   
      m2 = m2 / sd(m2, na.rm=TRUE)
      pca.out = pca_basic( cm=cm, indat=m2, nfactors=3 )

      scores = data.frame( id=rownames(m), pca1=pca.out$scores[, "PC1"], pca2=pca.out$scores[, "PC2"], pca3=pca.out$scores[, "PC3"], stringsAsFactors=FALSE )
      set = merge(set, scores, by="id", all.x=T, all.y=F, sort=FALSE)

      read_write_fast( pca.out, file=fn.pca)

      if (0) {

        toplot = as.data.frame( pca.out$loadings )
        toplot$vern = taxonomy.recode( from="spec", to="taxa", tolookup=rownames( toplot ) )$vern

        plot( PC2 ~ PC1, toplot, type="n")
        text( PC2 ~ PC1, labels=vern, data=toplot )

      }

      # Correpsondence analysis (on presence-absence)
      require(vegan)
      n = m[] * 0
      n[ which( m > cthreshold ) ] = 1
      ord = cca( n )
      sp.sc = scores(ord, choices=c(1:3))$species
      si.sc = scores(ord, choices=c(1:3))$sites
      scores = data.frame( id=as.character(rownames(si.sc)), ca1=as.numeric(si.sc[,1]), ca2=as.numeric(si.sc[,2]), ca3=as.numeric(si.sc[,3]) )
      variances=  ord$CA$eig[1:10]/sum(ord$CA$eig)*100
      set = merge(set, scores, by="id", all.x=T, all.y=F, sort=F)
      ca.out = list( scores=scores, ca=ord, variances=variances )

      read_write_fast( ca.out, file=fn.ca)
      read_write_fast( set, file=fn.set )

      return (fn.set)
    }



    # -----------------------



    if (DS %in% c( "speciescomposition", "speciescomposition.redo" ) ) {
      # remove dups and in planar coords
      fn = file.path( ddir, paste( "speciescomposition", infix, "rdz", sep=".") )

			if (DS=="speciescomposition") {
        SC = NULL
        if (file.exists( fn) ) loSC = read_write_fast( fn )
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

      ii = which( is.finite( rowSums( SC[,c( "pca1", "pca2",  "pca3", "ca1", "ca2", "ca3" ) ] ) ) )
      if (length(ii) == 0) stop( "No data .. something went wrong")
      SC = SC[ii,]

      read_write_fast( SC, file=fn )
			return (fn)
		}


    # -------------

    if ( DS=="areal_units_input" ) {

      
      outdir = file.path( p$data_root, "modelled", p$carstm_model_label ) 
      fn = file.path( outdir, "areal_units_input.rdz"  )
      if ( !file.exists(outdir)) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )

      xydata = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          xydata = read_write_fast( fn)
          return( xydata )
        }
      }
      xydata = speciescomposition_db( p=p, DS="speciescomposition"  )  #

      names(xydata)[which(names(xydata)=="z.mean" )] = "z"
      xydata = xydata[ geo_subset( spatial_domain=p$spatial_domain, Z=xydata ) , ] # need to be careful with extrapolation ...  filter depths
      xydata = xydata[ , c("lon", "lat", "yr" )]
      read_write_fast(xydata, file=fn )
      return( xydata )
    }


    # -------------


    if ( DS=="carstm_inputs") {

      # prediction surface
      crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
      if (is.null(sppoly))  sppoly = areal_units( p=p )  # will redo if not found
      sppoly = st_transform(sppoly, crs=crs_lonlat )
      areal_units_fn = attributes(sppoly)[["areal_units_fn"]]
      
      # over-ride
      p$variabletomodel = "speciescomposition"  # force to be generic to share across variables 

      fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs", areal_units_fn, sep="_") )
      if (p$carstm_inputs_prefilter =="rawdata") {
        fn = file.path( p$modeldir, p$carstm_model_label, paste("carstm_inputs_rawdata", areal_units_fn, sep="_") )
      }

      # inputs are shared across various secneario using the same polys
      #.. store at the modeldir level as default
      outputdir = dirname( fn )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      M = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          M = read_write_fast( fn)
          M = M[ which( M$yr %in% p$yrs), ]
          return( M )
        }
      }
      # message( "Generating carstm_inputs ... ")


      # do this immediately to reduce storage for sppoly (before adding other variables)

      M = speciescomposition_db( p=p, DS="speciescomposition"  )
      setDT(M)
      
      vars_to_retain = c("pca1", "pca2", "pca3", "ca1", "ca2", "ca3", "gear", "data.source", "vessel" , 
        "t", "z", "substrate.grainsize", "id") # , "sal", "oxyml" )


      # # INLA does not like duplicates ... causes optimizer to crash frequently
      oo = which(duplicated( M[, p$variabletomodel] ))
      if ( length(oo)> 0) {
        eps = exp( log( .Machine$double.eps ) / 2)  # ~ 1.5e-8
        M[oo, p$variabletomodel]  = M[oo, p$variabletomodel]  + runif( length(oo), -eps, eps )
      }

      M = planar2lonlat(M, proj.type=p$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
      
      ii = which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] )
      M = M[ ii, ]

      names(M)[which(names(M)=="yr") ] = "year"
      M = M[ which(M$year %in% p$yrs), ]
      M$tiyr = lubridate::decimal_date ( M$timestamp )
      M$dyear = M$tiyr - M$year
      
      M = carstm_prepare_inputdata( p=p, M=M, sppoly=sppoly, 
        NA_remove=FALSE,
        vars_to_retain=vars_to_retain,
        vars_to_drop ="speciescomposition" )  # drop dummy variable

      jj = which(is.na(M$id))
      M$id[jj] = "dummyvalue"

      jj = which(!is.finite(M$t))
      if (length(jj) > 0 ) {
        M$t[jj] = median( M$t[-jj] )
      }

      jj = which(!is.finite(M$substrate.grainsize))
      if (length(jj) > 0 ) {
        M$substrate.grainsize[jj] = median( M$substrate.grainsize[-jj] )
      }
      M$log.substrate.grainsize = log( M$substrate.grainsize )

      M$space = match( M$AUID, sppoly$AUID) # for bym/car .. must be numeric index matching neighbourhood graphs
      M$space_time = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 
      M$space_cyclic = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 

      M$time = match( M$year, p$yrs ) # copy for space_time component .. for groups, must be numeric index
      M$time_space = M$time    

      M$cyclic = match( M$dyri, discretize_data( span=c( 0, 1, p$nw) )  ) # match midpoints
      M$cyclic_space = M$cyclic # copy cyclic for space - cyclic component .. for groups, must be numeric index
    
      attr( M, "proj4string_planar" ) =  p$aegis_proj4string_planar_km
      attr( M, "proj4string_lonlat" ) =  projection_proj4string("lonlat_wgs84")

      read_write_fast( M, file=fn )

      M = M[ which( M$yr %in% p$yrs), ]

      return( M )
    }



  } # end function



