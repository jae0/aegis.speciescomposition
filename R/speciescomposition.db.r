
  speciescomposition.db = function( DS="", p=NULL ) {

    ddir = project.datadirectory( "aegis", "speciescomposition" )
    dir.create( ddir, showWarnings=FALSE, recursive=TRUE )

    infix = paste( p$spatial_domain,  p$taxa, sep=".")

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

      sc = survey.db( DS="cat" ,p=p)  # species catch

      sc = sc[ which(is.finite( sc$zn ) ), ]
      sc = sc[ , c("id", "spec_bio", "zn" ) ]  # zscore-transformed into 0,1
      sc = sc[ , c("id", "zn","spec_bio" ) ]  # zscore-transformed into 0,1

      set = survey.db( DS="set", p=p) # trip/set loc information

      set = set[ ,  c("id", "yr", "dyear", "sa", "lon", "lat", "t", "z", "timestamp" ) ]
      set = na.omit( set ) # all are required fields

      # filter area
      igood = which( set$lon >= p$corners$lon[1] & set$lon <= p$corners$lon[2]
              &  set$lat >= p$corners$lat[1] & set$lat <= p$corners$lat[2] )
      set = set[igood, ]

      # filter species
      # sc$spec = taxonomy.parsimonious( spec=sc$spec )

     # browser()

      isc = taxonomy.filter.taxa( sc$spec_bio, method=p$taxa, outtype="internalcodes" )
      set = set[ which( set$id %in% unique( sc$id[isc]) ),]

      # .. data loss due to truncation is OK
      # ... smallest abundance adds little information to ordinations
      k = 1e3         # a large constant number to make xtabs work  but not too large as truncation is desired
      sc$zn = as.integer( sc$zn*k )
      m = xtabs( zn ~ as.factor(id) + as.factor(spec_bio), data=sc ) /k

      # remove low counts (absence) in the timeseries  .. species (cols) only
      cthreshold = 0.05 * k  # quantiles to be removed

      finished = F
      while( !(finished) ) {
        i = unique( which(rowSums(m) == 0 ) )
        j = unique( which(colSums(m) <= cthreshold ) )
        if ( ( length(i) == 0 ) & (length(j) == 0 ) ) finished=T
        if (length(i) > 0 ) m = m[ -i , ]
        if (length(j) > 0 ) m = m[ , -j ]
      }

      # PCA
      # no need to correct for gear types/surveys .. assuming no size-specific bias .. perhaps wrong but simpler
      corel = cor( m, use="pairwise.complete.obs" ) # set up a correlation matrix ignoring NAs
      corel[ is.na(corel) ] = 0
      s = svd(corel)  # eigenanalysis via singular value decomposition

      # matrix.multiply = function (x, y, nfac=2){
      #   ndat = dim(x)[1]
      #   z = matrix(0, nrow=ndat, ncol = nfac)
      #   for (j in 1:nfac) {
      #     for (i in 1:ndat) {
      #       z[i,j] = sum ( x[i,] * t(y[,j]), na.rm=T )
      #     }
      #   }
      #   return (z)
      # }

      # scores = matrix.multiply (m, s$v)  # i.e., b %*% s$v  .. force a multiplication ignoring NAs
      m[which(!is.finite(m))] = 0
      scores = m %*% s$v  # i.e., b %*% s$v  .. force a multiplication ignoring NAs

      evec = s$v
      ev = s$d
      x = cbind( scores[,1] / sqrt(ev[1] ), scores[,2] / sqrt( ev[2]) )
      y = cbind( evec[,1] * sqrt(ev[1] ) , evec[,2] * sqrt( ev[2]) )
      rownames(y) = colnames(m)

      scores = data.frame( id=rownames(m), pca1=as.numeric(x[,1]), pca2=as.numeric(x[,2]), stringsAsFactors=F )
      set = merge(set, scores, by="id", all.x=T, all.y=F, sort=F)
      pca.out = list( scores=scores, eignenvectors=evec, eigenvalues=ev, cscores=y )
      save( pca.out, file=fn.pca, compress=T)


      # Correpsondence analysis
      require(vegan)
        n = m * 0
        n[ which(m>0) ] = 1
        ord = cca( n )
        sp.sc = scores(ord)$species
        si.sc = scores(ord)$sites
        scores = data.frame( id=as.character(rownames(si.sc)), ca1=as.numeric(si.sc[,1]), ca2=as.numeric(si.sc[,2]) )
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

      SC = speciescomposition.db( DS="speciescomposition.ordination", p=p )
      SC = lonlat2planar( SC, proj.type=p$aegis_proj4string_planar_km )
      SC$lon = SC$lat = NULL

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

      save( SC, file=fn, compress=T )
			return (fn)
		}



    # -----------------------




    if ( DS=="carstm_inputs") {

      fn = file.path( p$modeldir, paste( "speciescomposition", "carstm_inputs", p$auid,
        p$inputdata_spatial_discretization_planar_km,
        round(p$inputdata_temporal_discretization_yr, 6),
        "rdata", sep=".") )

      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( M )
        }
      }
      message( "Generating carstm_inputs ... ")

      # prediction surface
      sppoly = areal_units( p=p )  # will redo if not found

      crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))

      # do this immediately to reduce storage for sppoly (before adding other variables)
      M = speciescomposition.db( p=p, DS="speciescomposition" )

      M = M[ which(M$yr %in% p$yrs), ]
      M$tiyr = lubridate::decimal_date ( M$date )

      # globally remove all unrealistic data
  #    keep = which( M[,p$variabletomodel] >= -3 & M[,p$variabletomodel] <= 25 ) # hard limits
  #    if (length(keep) > 0 ) M = M[ keep, ]

      # p$quantile_bounds_data = c(0.0005, 0.9995)
      if (exists("quantile_bounds_data", p)) {
        TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
        keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
      }

      M$dyear = M$tiyr - M$yr
      M$dyear = discretize_data( M$dyear, seq(0, 1, by=p$inputdata_temporal_discretization_yr), digits=6 )

      # reduce size
      M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
      # levelplot(z.mean~plon+plat, data=M, aspect="iso")

      M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      M$lon = NULL
      M$lat = NULL
      M = M[ which(is.finite(M$StrataID)),]
      M$StrataID = as.character( M$StrataID )  # match each datum to an area

      names(M)[which(names(M)==paste(p$variabletomodel, "mean", sep=".") )] = p$variabletomodel

      M$tag = "observations"

      APS = as.data.frame(sppoly)
      APS$StrataID = as.character( APS$StrataID )
      APS$tag ="predictions"
      APS[,p$variabletomodel] = NA


      pB = aegis.bathymetry::bathymetry_parameters( p=p, project_class="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
      pS = aegis.substrate::substrate_parameters( p=p, project_class =="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
      pT = aegis.temperature::temperature_parameters( p=p, project_class =="carstm_auid" ) # transcribes relevant parts of p to load       TI = carstm_model ( p=pT, DS="carstm_modelled" )  # unmodeled!
     
    
      M$plon = round(M$plon / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km
      M$plat = round(M$plat / p$inputdata_spatial_discretization_planar_km + 1 ) * p$inputdata_spatial_discretization_planar_km

      BI = carstm_model ( p=pB, DS="aggregated_data" )  # unmodeled!
      jj = match( as.character( M$StrataID), as.character( BI$StrataID) )
      M[, pB$variabletomodel] = BI[jj, paste(pB$variabletomodel,"predicted",sep="." )]
      jj =NULL
      BI = NULL

      SI = carstm_model ( p=pS, DS="aggregated_data" )  # unmodeled!
      jj = match( as.character( M$StrataID), as.character( SI$StrataID) )
      M[, pS$variabletomodel] = SI[jj, paste(pS$variabletomodel,"predicted",sep="." )]
      jj =NULL
      SI = NULL

      TI = carstm_model ( p=pT, DS="aggregated_data" )  # unmodeled!
      jj = match( as.character( M$StrataID), as.character( TI$StrataID) )
      M[, pT$variabletomodel] = TI[jj, paste(pT$variabletomodel,"predicted",sep="." )]
      jj =NULL
      TI = NULL



      BI = carstm_model ( p=pB, DS="carstm_modelled" )  # unmodeled!
      jj = match( as.character( APS$StrataID), as.character( BI$StrataID) )
      APS[, pB$variabletomodel] = BI[jj, paste(pB$variabletomodel,"predicted",sep="." )]
      jj =NULL
      BI = NULL

      SI = carstm_model ( p=pS, DS="carstm_modelled" )  # unmodeled!
      jj = match( as.character( APS$StrataID), as.character( SI$StrataID) )
      APS[, pS$variabletomodel] = SI[jj, paste(pS$variabletomodel,"predicted",sep="." )]
      jj =NULL
      SI = NULL

      TI = carstm_model ( p=pT, DS="carstm_modelled" )  # unmodeled!
      jj = match( as.character( APS$StrataID), as.character( TI$StrataID) )
      APS[, pT$variabletomodel] = TI[jj, paste(pT$variabletomodel,"predicted",sep="." )]
      jj =NULL
      TI = NULL

      vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel,  pT$variabletomodel, "tag", "StrataID" )
      APS = APS[, vn]

      # expand APS to all time slices
      n_aps = nrow(APS)
      APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
      names(APS) = c(vn, "tiyr")

      M$tiyr = M$yr + M$dyear
      M = rbind( M[, names(APS)], APS )
      APS = NULL

      M$strata  = as.numeric( M$StrataID)
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
      M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
      M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

      M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
      M$tiyr2 = M$tiyr  # a copy

      M$year = floor(M$tiyr)
      M$dyear  =  factor( as.character( trunc(  (M$tiyr - M$year )/ p$tres )*p$tres), levels=p$dyears)

      save( M, file=fn, compress=TRUE )
      return( M )
    }


  } # end function



