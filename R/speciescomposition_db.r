
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

      sc = survey_db( DS="cat" ,p=p)  # species catch

      sc = sc[ which(is.finite( sc$zn ) ), ]
      sc = sc[ , c("id", "spec_bio", "zn" ) ]  # zscore-transformed into 0,1
      sc = sc[ , c("id", "zn","spec_bio" ) ]  # zscore-transformed into 0,1

      set = survey_db( DS="set", p=p) # trip/set loc information

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

      SC = speciescomposition_db( DS="speciescomposition.ordination", p=p )
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
      xydata = planar2lonlat(xydata, p$areal_units_proj4string_planar_km)  # should not be required but to make sure
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

      fn = carstm_filenames( p=p, returntype="carstm_inputs", areal_units_fn=areal_units_fn )
      if (!p$carstm_inputs_aggregated) {
        fn = carstm_filenames( p=p, returntype="carstm_inputs_rawdata", areal_units_fn=areal_units_fn )
      }

      # inputs are shared across various secneario using the same polys
      #.. store at the modeldir level as default
      outputdir = dirname( fn )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      M = NULL
      if (!redo)  {
        if (file.exists(fn)) {
          load( fn)
          return( M )
        }
      }
      # message( "Generating carstm_inputs ... ")


      # do this immediately to reduce storage for sppoly (before adding other variables)

      M = speciescomposition_db( p=p, DS="speciescomposition"  )

      # globally remove all unrealistic data
        # p$quantile_bounds_data = c(0.0005, 0.9995)
      if (exists("quantile_bounds_data", p)) {
        TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds_data, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
        keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
        if (length(keep) > 0 ) M = M[ keep, ]
          # this was -1.7, 21.8 in 2015
      }

      M = planar2lonlat(M, proj.type=p$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
      M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]

      names(M)[which(names(M)=="yr") ] = "year"
      M = M[ which(M$year %in% p$yrs), ]
      M$tiyr = lubridate::decimal_date ( M$timestamp )
      M$dyear = M$tiyr - M$year

      M$AUID = st_points_in_polygons(
        pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname = "AUID"
      )
      M = M[ which(!is.na(M$AUID)),]
      M$AUID = as.character( M$AUID )  # match each datum to an area


      # --------------------------
      # bathymetry observations  lookup
      pB = bathymetry_parameters( project_class="core"  )
      vnB = pB$variabletomodel
      if ( !(exists(vnB, M ))) {
        vnB2 = paste(vnB, "mean", sep=".")
        if ((exists(vnB2, M ))) {
          names(M)[which(names(M) == vnB2 )] = vnB
        } else {
          M[,vnB] = NA
        }
      }
      iM = which(!is.finite( M[, vnB] ))
      if (length(iM > 0)) {
        M[iM, vnB] = bathymetry_lookup( LOCS=M[iM, c("lon", "lat")], spatial_domain=p$spatial_domain, lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
      }

      M = M[ is.finite(M[ , vnB]  ) , ]


      # must go after depths have been finalized
      if (p$carstm_inputs_aggregated) {
        if ( exists("spatial_domain", p)) {
          M = M[ geo_subset( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...  filter depths
        }
      }

 
      # --------------------------
      # substrate observations  lookup
      pS = substrate_parameters( project_class="core"  )
      if (!(exists(pS$variabletomodel, M ))) M[,pS$variabletomodel] = NA
      iM = which(!is.finite( M[, pS$variabletomodel] ))
      if (length(iM > 0)) {
        M[iM, pS$variabletomodel] = substrate_lookup( LOCS=M[iM, c("lon", "lat")], spatial_domain=p$spatial_domain, lookup_from="core", lookup_to="points" , lookup_from_class="aggregated_data" ) # core=="rawdata"
      }

      M = M[ is.finite(M[ , pS$variabletomodel]  ) , ]


      # --------------------------
      # temperature observations lookup
      pT = temperature_parameters( project_class="core"  )
      if (!(exists(pT$variabletomodel, M ))) M[,pT$variabletomodel] = NA
      iM = which(!is.finite( M[, pT$variabletomodel] ))
      if (length(iM > 0)) {
        M[iM, pT$variabletomodel] = temperature_lookup_rawdata( spatial_domain=p$spatial_domain, M=M[iM, c("lon", "lat", "timestamp")], sppoly=sppoly )
      }


      M$lon = M$lat = M$plon = M$plat = NULL

      M = M[ which(is.finite(M[, pB$variabletomodel] )),]
      M = M[ which(is.finite(M[, pS$variabletomodel] )),]
      M = M[ which(is.finite(M[, pT$variabletomodel] )),]

      M$tag = "observations"


      # end observations
      # ----------

      # predicted locations (APS)


      region.id = slot( slot(sppoly, "nb"), "region.id" )
      APS = st_drop_geometry(sppoly)

      APS$AUID = as.character( APS$AUID )
      APS$tag ="predictions"
      APS[,p$variabletomodel] = NA


      vnmod = pB$variabletomodel
      vnp = paste(vnmod, "predicted", sep=".")
      # vnps = paste(vnmod, "predicted_se", sep=".")


      

      if (p$carstm_inputdata_model_source$bathymetry=="carstm") {
        LU = carstm_model( p=pB, DS="carstm_modelled_summary" ) # to load exact sppoly, if present
        LU_sppoly = areal_units( p=pB )  # default poly

        if (is.null(LU)) {
          message("Exactly modelled surface not found, estimating from default run...")
          pBD = bathymetry_parameters( project_class="carstm" ) # choose "default" full bathy carstm run and re-estimate:
          LU = carstm_model( p=pBD, DS="carstm_modelled_summary" )
          LU_sppoly = areal_units( p=pBD )  # default poly
        }

        bm = match( LU_sppoly$AUID, LU$AUID )

        LU_sppoly[[ vnp ]] = LU[[ vnp ]] [ bm ]
        # LU_sppoly[[ vnps ]] = LU[[ vnps ]] [ bm ]
        bm = NULL
        LU = NULL

        # now rasterize and estimate
        raster_template = raster( sppoly, res=p$areal_units_resolution_km, crs=st_crs( sppoly ) ) # +1 to increase the area
        # transfer the coordinate system to the raster
        LU_sppoly = sf::st_transform( as( LU_sppoly, "sf" ), crs=st_crs(LU_sppoly) )  # B is a carstm LU
        LU_sppoly = fasterize::fasterize( LU_sppoly, raster_template, field=vnp )
        sppoly[[ vnp ]] = sp::over( sppoly, LU_sppoly[, vnp ], fn=mean, na.rm=TRUE )
        # sppoly[[ vnp s]] = sp::over( sppoly, LU_sppoly[, vnp ], fn=sd, na.rm=TRUE )
        LU_sppoly = NULL
        raster_template = NULL
      }


      if (p$carstm_inputdata_model_source$bathymetry %in% c("stmv", "hybrid")) {
        pBD = bathymetry_parameters( project_class=p$carstm_inputdata_model_source$bathymetry )  # full default
        vnmod = pBD$variabletomodel
        vnp = paste(vnmod, "predicted", sep=".")
        # vnps = paste(vnmod, "predicted_se", sep=".")

        LU = bathymetry_db( p=pBD, DS="baseline", varnames=c(vnmod, "plon", "plat") )
        LU = planar2lonlat(LU, pBD$aegis_proj4string_planar_km)
        LU = sf::st_as_sf( LU[, c( vnmod, "lon", "lat")], coords=c("lon", "lat") )
        st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
        LU = sf::st_transform( LU, crs=st_crs(sppoly) )
        sppoly[[ vnp ]] = aggregate( LU[, vnmod], sppoly, mean, na.rm=TRUE )[[vnmod]]
        # sppoly[[ vnps ]] = aggregate( LU[, vnmod], sppoly, sd, na.rm=TRUE )[[vnmod]]
      }

      iAS = match( as.character( APS$AUID), as.character( sppoly$AUID ) )
      APS[, pB$variabletomodel] = sppoly[[ paste(pB$variabletomodel, "predicted", sep=".") ]] [iAS]
      iAS =NULL
      #sppoly = NULL
      gc()


      if (p$carstm_inputdata_model_source$substrate=="carstm") {
        LU = carstm_model( p=pS, DS="carstm_modelled_summary" ) # to load exact sppoly, if present
        # sppoly = areal_units( p=pS )  # default poly no need to reload

        if (is.null(LU)) {
          message("Exactly modelled surface not found, estimating from default run... not fully tested")
          pSD = substrate_parameters( project_class="carstm" ) # choose "default" full bathy carstm run and re-estimate:
          vnmod = pSD$variabletomodel
          vnp = paste(vnmod, "predicted", sep=".")
          # vnps = paste(vnmod, "predicted_se", sep=".")

          LU = carstm_model( p=pSD, DS="carstm_modelled_summary" )
          LU_sppoly = areal_units( p=pSD )  # default poly
          iLU = match( LU_sppoly$AUID, LU$AUID )
          LU_sppoly[[vnp]] = LU[[vnp]][ iLU ]

          # LU_sppoly[[vnps]] = LU[[vnps]][ iLU ]
          # now rasterize and restimate to desired sppoly
          raster_template = raster( sppoly, res=p$areal_units_resolution_km, crs=st_crs( sppoly ) ) # +1 to increase the area
          # transfer the coordinate system to the raster
          LU = sf::st_transform( as(LU, "sf"), crs=st_crs(sppoly) )  # B is a carstm sppoly
          LU = fasterize::fasterize( LU, raster_template, field=vnp )
          sppoly[[vnp]] = sp::over( sppoly, LU[, vnp ], fn=mean, na.rm=TRUE )
          # sppoly[[vnps]] = sp::over( sppoly, LU[, vnp ], fn=sd, na.rm=TRUE )
        }
      }

      if (p$carstm_inputdata_model_source$substrate %in% c("stmv", "hybrid")) {
        pSD = substrate_parameters( project_class=p$carstm_inputdata_model_source$substrate )  # full default
        vnmod = pSD$variabletomodel
        vnp = paste(vnmod, "predicted", sep=".")
        # vnps = paste(vnmod, "predicted_se", sep=".")

        LU = substrate_db( p=pSD, DS="complete"  )
        LU_coords = bathymetry_db( p=bathymetry_parameters( spatial_domain=p$spatial_domain, project_class=p$carstm_inputdata_model_source$substrate  ), DS="baseline"  )
        LU = cbind(LU, LU_coords)
        LU = planar2lonlat(LU, pSD$aegis_proj4string_planar_km)
        LU = sf::st_as_sf(  LU[,  c(vnmod, "lon", "lat") ], coords=c("lon", "lat") )
        st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
        LU = sf::st_transform( LU, crs=st_crs(sppoly) )
        sppoly[[ vnp ]] = aggregate( LU[, vnmod], sppoly, mean, na.rm=TRUE )[[vnmod]]
        # sppoly$[[ vnps ]] = aggregate( LU[, vnmod], sppoly, sd, na.rm=TRUE )[[vnmod]]
      }

      jj = match( as.character( APS$AUID), as.character( LU$AUID) )
      APS[, pS$variabletomodel] = LU[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
      jj =NULL
      LU = NULL

      # to this point APS is static, now add time dynamics (teperature)
      # ---------------------

      vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "AUID" )
      APS = APS[, vn]

      # expand APS to all time slices
      n_aps = nrow(APS)
      APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
      names(APS) = c(vn, "tiyr")
      APS$year = aegis_floor( APS$tiyr)
      APS$dyear = APS$tiyr - APS$year

      if (p$carstm_inputdata_model_source$temperature=="carstm") {

        LU = carstm_model( p=pT, DS="carstm_modelled_summary" ) # to load exact sppoly, if present
        LU_sppoly = areal_units( p=pT )  # default poly

        if (is.null(LU)) {
          message("Exactly modelled surface not found, estimating from default run...")
          pTD = temperature_parameters( project_class="carstm" ) # choose "default" full bathy carstm run and re-estimate:
          LU = carstm_model( p=pTD, DS="carstm_modelled_summary" )
          LU_sppoly = areal_units( p=pTD )  # default poly
        }

        bm = match( LU_sppoly$AUID, LU$AUID )
        LU_sppoly$t.predicted = LU$t.predicted[ bm ]
        # LU_sppoly$t.predicted_se = LU$t.predicted_se[ bm ]
        bm = NULL
        LU = NULL

          if (! "POSIXct" %in% class(timestamp)  ) timestamp = as.POSIXct( timestamp, tz=tz, origin=lubridate::origin  )
          tstamp = data.frame( yr = lubridate::year(timestamp) )
          tstamp$dyear = lubridate::decimal_date( timestamp ) - tstamp$yr
          timestamp_map = array_map( "ts->2", tstamp[, c("yr", "dyear")], dims=c(p$ny, p$nw), res=c( 1, 1/p$nw ), origin=c( min(p$yrs), 0) )
    # TODO
    stop("not finished ... must addd time lookup")  # TODO
    # TODO

          dindex = cbind(locs_index, timestamp_map ) # check this

          # convert to raster then match
          require(raster)
          raster_template = raster(extent(locs)) # +1 to increase the area
          res(raster_template) = p$areal_units_resolution_km  # crs usually in meters, but aegis's crs is in km
          crs(raster_template) = projection(locs) # transfer the coordinate system to the raster
          Bsf = sf::st_transform( as(LU_sppoly, "sf"), crs=CRS(proj4string(locs)) )  # LU_sppoly is a carstm sppoly
          Boo = as(LU_sppoly, "SpatialPolygonsDataFrame")
          for (vn in Bnames) {
            # Bf = fasterize::fasterize( Bsf, raster_template, field=vn )
  #          vn2 = paste(vn, "sd", sep="." )
            locs[,vn] = st_polygons_in_polygons( locs, Bf[,vn], fn=mean, na.rm=TRUE )
  #          locs[,vn2] = st_polygons_in_polygons( locs, Bf[,vn], fn=sd, na.rm=TRUE )
          }
          vnames = intersect( names(LU_sppoly), vnames )
          if ( length(vnames) ==0 ) vnames=names(LU_sppoly) # no match returns all
          return(locs[,vnames])

      }

      if (p$carstm_inputdata_model_source$temperature %in% c("stmv", "hybrid")) {
        pTD = temperature_parameters( project_class=p$carstm_inputdata_model_source$temperature )  # full default
        vnmod = pTD$variabletomodel
        vnp = paste(vnmod, "predicted", sep=".")
        # vnps = paste(vnmod, "predicted_se", sep=".")

for (ti in p$nt){

}

        LU = temperature_db( p=pTD, DS="complete"  )
        LU_coords = bathymetry_db( p=bathymetry_parameters( spatial_domain=p$spatial_domain, project_class=p$carstm_inputdata_model_source$temperature  ), DS="baseline"  )

        LU = cbind(LU, LU_coords)
        LU = planar2lonlat(LU, pTD$aegis_proj4string_planar_km)
        LU = sf::st_as_sf(  LU[,  c(vnmod, "lon", "lat") ], coords=c("lon", "lat") )
        st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
        LU = sf::st_transform( LU, crs=st_crs(sppoly) )
        sppoly[[ vnp ]] = aggregate( LU[, vnmod], sppoly, mean, na.rm=TRUE )[[vnmod]]
        # sppoly$[[ vnps ]] = aggregate( LU[, vnmod], sppoly, sd, na.rm=TRUE )[[vnmod]]
      }


      LU = LU[[ paste(pT$variabletomodel,"predicted",sep="." )]]

      auid_map = match( APS$AUID, dimnames(LU)$AUID )
      year_map = match( as.character(APS$year), dimnames(LU)$year )

      dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
      dyear_map = as.numeric( cut( APS$dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )

      dindex = cbind(auid_map, year_map, dyear_map )

      APS[, pT$variabletomodel] = LU[ dindex]
      jj =NULL
      LU = NULL

      M = rbind( M[, names(APS)], APS )
      APS = NULL

      M$auid = match( M$AUID, region.id )

      M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

      M$year = aegis_floor( M$tiyr)
      M$year_factor = as.numeric( factor( M$year, levels=p$yrs))
      M = M[ is.finite(M$year_factor), ]
      M$dyear =  M$tiyr - M$year   # revert dyear to non-discretized form

      if (0) {
        M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
        M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
        M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )
        M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )
        M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)
      }

      save( M, file=fn, compress=TRUE )
      return( M )
    }



  } # end function



