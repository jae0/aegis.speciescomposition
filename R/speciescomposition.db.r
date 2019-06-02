
  speciescomposition.db = function( DS="", p=NULL ) {

    ddir = project.datadirectory( "aegis", "speciescomposition" )
    dir.create( ddir, showWarnings=FALSE, recursive=TRUE )

    infix = paste( p$spatial.domain,  p$taxa, sep=".")

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
      SC = lonlat2planar( SC, proj.type=p$internal.crs )
      SC$lon = SC$lat = NULL

      oo = which(!is.finite( SC$plon+SC$plat ) )
      if (length(oo)>0) SC = SC[ -oo , ]  # a required field for spatial interpolation

      yrs = sort( unique( SC$yr ) )
      # check for duplicates
      for ( y in yrs ) {
        yy = which (SC$yr == y)
        ii = which( duplicated( SC$id[yy] ) )

        if (length(ii) > 0) {
          print( "The following sets have duplicated positions. The first only will be retained" )
          print( SC[yy,] [ duplicates.toremove( SC$id[yy] ) ] )
          SC = SC[ - ii,]
        }
      }
  
      save( SC, file=fn, compress=T )
			return (fn)
		}

  } # end function


