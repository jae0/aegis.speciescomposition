speciescomposition_lookup = function( M, spatial_domain=NULL, sppoly=NULL,   tz="America/Halifax", lookup_mode="stmv", vnames=NULL ) {
  # lookup from rawdata

  if (is.null(spatial_domain))  {
    pPC = speciescomposition_parameters(  project_class="core"  )
  } else {
    pPC = speciescomposition_parameters( spatial_domain=spatial_domain, project_class="core"  )
  }

  crs_lonlat =  st_crs(projection_proj4string("lonlat_wgs84"))

  M = as.data.frame( M )
  names(M) = c("lon", "lat", "timestamp")
  M = lonlat2planar(M, proj.type=pPC$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
  if (! "POSIXct" %in% class(M$timestamp)  ) M$timestamp = lubridate::date_decimal( M$timestamp, tz=tz  )
  M$yr = lubridate::year(M$timestamp)
  M$dyear = lubridate::decimal_date( M$timestamp ) - M$yr

  LU = speciescomposition_db ( p=pPC, year.assessment=max(pPC$yrs), DS="aggregated_data" )  # raw data
  names(LU)[ which(names(LU) =="speciescomposition.mean") ] = vnames
  LU = LU[ which( LU$lon > pPC$corners$lon[1] & LU$lon < pPC$corners$lon[2]  & LU$lat > pPC$corners$lat[1] & LU$lat < pPC$corners$lat[2] ), ]
  LU = lonlat2planar(LU, proj.type=pPC$aegis_proj4string_planar_km)

  LUT_map = array_map( "ts->1", LU[,c("yr", "dyear")], dims=c(pPC$ny, pPC$nw), res=c( 1, 1/pPC$nw ), origin=c( min(pPC$yrs), 0) )
  LUS_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=pPC$gridparams )

  T_map = array_map( "ts->1", M[, c("yr", "dyear")], dims=c(pPC$ny, pPC$nw), res=c( 1, 1/pPC$nw ), origin=c( min(pPC$yrs), 0) )
  M_map = array_map( "xy->1", M[, c("plon","plat")], gridparams=pPC$gridparams )

  iLM = match( paste(M_map, T_map, sep="_"), paste(LUS_map, LUT_map, sep="_") )
  M[ , vnames ] = LU[ iLM, paste(vnames, "mean", sep="." ) ]

  gc()

  if (!is.null(sppoly)) {
        # if any still missing then use a mean by AUID
    ii = NULL
    ii =  which( !is.finite(M[ , vnames ]))
    if (length(ii) > 0) {
      if (!exists("AUID", M)) {
        M_AUID = st_points_in_polygons(
          pts = st_as_sf( M[ii,], coords=c("lon","lat"), crs=crs_lonlat ),
          polys = sppoly[, "AUID"],
          varname = "AUID"
        )
        M_AUID = as.character( M_AUID )  # match each datum to an area
      }
      T_map = array_map( "ts->1", M[ii, c("yr", "dyear")], dims=c(pPC$ny, pPC$nw), res=c( 1, 1/pPC$nw ), origin=c( min(pPC$yrs), 0) )
      M_uid =  paste(M_AUID,  T_map, sep=".")

      LU$AUID = st_points_in_polygons(
        pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname="AUID"
      )
      
      LUT_map = array_map( "ts->1", LU[,c("yr", "dyear")], dims=c(pPC$ny, pPC$nw), res=c( 1, 1/pPC$nw ), origin=c( min(pPC$yrs), 0) )
      LU$uid = paste(LU$AUID, LUT_map, sep=".")
      LU = tapply( LU[, paste(vnames, "mean", sep="." )], LU$uid, FUN=median, na.rm=TRUE )
      jj = match( as.character( M_uid), as.character( names(LU )) )
      
      M[ ii, vnames ] = LU[jj]
    }
  }
 

  return( M[ , vnames ] )

}
