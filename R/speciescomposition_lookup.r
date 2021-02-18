speciescomposition_lookup = function( LOCS=NULL, AU_target=NULL, AU=NULL, 
  lookup_from="core", lookup_to="points", 
  FUNC=mean,  vnames="pca1", vnames_from=paste(vnames, "mean", sep="."), 
  lookup_from_class="aggregated_data", tz="America/Halifax" , year.assessment=NULL, pO=NULL ) {
  # lookup from rawdata
 
  message("need to check::  [match( APS$AUID, as.character( sppoly$AUID ) )] ")

  if (is.null(year.assessment)) year.assessment = max( lubridate::year(LOCS$timestamp) )

  if (is.null(pO)) pO = speciescomposition_parameters(  project_class=lookup_from, variabletomodel=vnames, year.assessment=year.assessment  )

  crs_lonlat =  st_crs(projection_proj4string("lonlat_wgs84"))

  if ( lookup_from %in% c("core") & lookup_to == "points" )  {
    # matching to point (LU ) to point (LOCS)
    # if any still missing then use stmv depths
    
    LOCS = lonlat2planar(LOCS, proj.type=pO$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    
    if (! "POSIXct" %in% class(LOCS$timestamp)  ) LOCS$timestamp =  lubridate::date_decimal( LOCS$timestamp, tz=tz )
    LOCS$yr = lubridate::year( LOCS$timestamp ) 
    LOCS$dyear = lubridate::decimal_date( LOCS$timestamp ) - LOCS$yr

    LU = speciescomposition_db ( p=pO, DS="speciescomposition" )  # raw data
    # LU = planar2lonlat(LU, proj.type=pO$aegis_proj4string_planar_km)
    # LU = lonlat2planar(LU, proj.type=pO$aegis_proj4string_planar_km)
    names(LU)[ which(names(LU) == vnames_from ) ] =  vnames

    LU_map = paste( 
      array_map( "xy->1", LU[,c("plon","plat")], gridparams=pO$gridparams ), 
      array_map( "ts->1", LU[,c("yr", "dyear")], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) ), 
      sep="_" 
    )

    LOCS_map = paste(
      array_map( "xy->1", LOCS[, c("plon","plat")], gridparams=pO$gridparams ), 
      array_map( "ts->1", LOCS[, c("yr", "dyear")], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) ), 
      sep="_"
    )

    LOCS[ , vnames ] = LU[ match( LOCS_map, LU_map ), vnames ]

    return( LOCS[ , vnames ] )
  }


  if ( lookup_from %in% c("core") & lookup_to == "areal_units" )  {
    # point (LU) -> areal unit (LOCS: AU/timestamp)
    #   $ pcaX.predicted             : num [1:190, 1:20, 1:10] 2.37 2.01 5.27 1.98 1.56 ...
    # ..- attr(*, "dimnames")=List of 3
    # .. ..$ AUID : chr [1:190] "106" "107" "108" "121" ...
    # .. ..$ year : chr [1:20] "1999" "2000" "2001" "2002" ...
    # .. ..$ dyear: chr [1:10] "0.05" "0.15" "0.25" "0.35" ...

    if (!exists("AU")) AU = areal_units( p=pO )
    if (!exists("AUID", AU)) AU$AUID = as.character(1:nrow(AU))
    if (!exists("AUID", LOCS) | !exists("timestamp", LOCS) ) stop("require AUID and timestamp in LOCS") 

    if (! "POSIXct" %in% class(LOCS$timestamp)  ) LOCS$timestamp =  lubridate::date_decimal( LOCS$timestamp, tz=tz )
    LOCS$yr = lubridate::year(LOCS$timestamp) 
    LOCS$dyear = lubridate::decimal_date( LOCS$timestamp ) - LOCS$yr


    LU = speciescomposition_db ( p=pO, DS="speciescomposition" )  # raw data
    # LU = lonlat2planar(LU, proj.type=pO$aegis_proj4string_planar_km)
    names(LU)[ which(names(LU) ==vnames_from) ] =  vnames
    LU = sf::st_as_sf( LU, coords=c("lon", "lat") )
    st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
    LU = sf::st_transform( LU, crs=st_crs(AU) )
    LU_map = paste( 
      st_points_in_polygons( pts=LU, polys=AU[, "AUID"], varname= "AUID" ), 
      array_map( "ts->1", st_drop_geometry( LU) [,c("yr", "dyear")], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) ), 
      sep="_"
    )
    
    LOCS = st_as_sf( LOCS, coords=c("lon","lat") )
    st_crs(LOCS) = st_crs( projection_proj4string("lonlat_wgs84") )
    LOCS = sf::st_transform( LOCS, crs=st_crs(AU) )
    LOCS_map =  paste( 
      st_points_in_polygons( pts=LOCS, polys = AU[, "AUID"], varname= "AUID" ),  
      array_map( "ts->1", st_drop_geometry(LOCS)[ , c("yr", "dyear") ], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) ), 
      sep="_"
    )

    LOCS_regridded = tapply( st_drop_geometry(LU)[, vnames], LU_map, FUN=FUNC, na.rm=TRUE )

    LOCS[ , vnames ] = LOCS_regridded[ match( LOCS_map, as.character( names( LOCS_regridded )) ) ]

    return( LOCS[, vnames] )
  }


  if ( lookup_from %in% c("stmv", "hybrid") & lookup_to == "points" )  {
    

    # matching to point (LU) to points (LOCS)
    LU = aegis_db( p=pO, DS="complete" )  # raw data
    LU = planar2lonlat(LU, proj.type=pO$aegis_proj4string_planar_km)
    LU_map = array_map( "xy->1", LU[, c("plon","plat")], gridparams=pO$gridparams )
    
    LOCS = lonlat2planar(LOCS, proj.type=pO$aegis_proj4string_planar_km) # get planar projections of lon/lat in km
    LOCS_map = array_map( "xy->1", LOCS[, c("plon","plat")], gridparams=pO$gridparams )
    LOCS_index = match( LOCS_map, LU_map )

    if (! "POSIXct" %in% class(LOCS$timestamp)  ) LOCS$timestamp =  lubridate::date_decimal( LOCS$timestamp, tz=tz )
    LOCS$yr = lubridate::year(LOCS$timestamp) 
    LOCS$dyear = lubridate::decimal_date( LOCS$timestamp ) - LOCS$yr

    # TIMESTAMP_index = array_map( "ts->2", LOCS [, c("yr", "dyear")], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) )
    TIMESTAMP_index = array_map( "ts->year_index", LOCS [, c("yr" )], dims=c(pO$ny ), res=c( 1  ), origin=c( min(pO$yrs) ) )


    return( LOCS[ cbind( LOCS_index, TIMESTAMP_index ) ] )
  }


   if ( lookup_from %in% c("stmv", "hybrid") & lookup_to == "areal_units" )  {
     # points (LU) -> areal units (LOCS)
     LU = aegis_db( p=pO, DS="complete" )  # raw data

     LU = planar2lonlat(LU, pO$aegis_proj4string_planar_km)
    
     LU = sf::st_as_sf( LU, coords=c("lon", "lat") )
     st_crs(LU) = st_crs( projection_proj4string("lonlat_wgs84") )
     LU = sf::st_transform( LU, crs=st_crs(LOCS) )

     if (!exists("AUID", AU_target)) AU_target$AUID = as.character(1:nrow(AU_target))

    if (! "POSIXct" %in% class(LOCS$timestamp)  ) LOCS$timestamp =  lubridate::date_decimal( LOCS$timestamp, tz=tz )
    LOCS$yr = lubridate::year(LOCS$timestamp) 
    LOCS$dyear = lubridate::decimal_date( LOCS$timestamp ) - LOCS$yr
     # TIMESTAMP_index = array_map( "ts->2", LOCS [, c("yr", "dyear")], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) )
    TIMESTAMP_index = array_map( "ts->year_index", LOCS [, c("yr" )], dims=c(pO$ny ), res=c( 1  ), origin=c( min(pO$yrs) ) )

 
 
    # now rasterize and re-estimate
    LOCS$AU_index = match( LOCS$AUID, LU$AUID  )    # assuming AUID's are consistent
    
    # AU_target .... must be sent ... <---------
    AU_target = sf::st_transform( AU_target, crs=st_crs(pO$aegis_proj4string_planar_km) )


  #   LOCS_AUID = st_points_in_polygons( pts=st_as_sf( LOCS, coords=c("lon","lat"), crs=crs_lonlat ), polys = LU[, "AUID"], varname= "AUID" )
    
  #   LOCS_map =  paste( 
  #     as.character( LOCS_AUID ),  
  #     array_map( "ts->1", LOCS[ , c("yr", "dyear") ], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) ), 
  #     sep="_"
  #   )

     LU_map = paste( 
       st_points_in_polygons( pts=LU, polys=AU_target[, "AUID"], varname= "AUID" ),
       array_map( "ts->1", LU[,c("yr", "dyear")], dims=c(pO$ny, pO$nw), res=c( 1, 1/pO$nw ), origin=c( min(pO$yrs), 0) ), 
       sep="_"
     )
    

  #   for (vn in vnames) {
  #     LOCS_regridded = tapply( LU[, vn], LU_uid, FUN=FUNC, na.rm=TRUE )
  #     LOCS[ , vn ] = LOCS_regridded[ match( LOCS_map, as.character( names( LOCS_regridded )) ) ]
  #   }
  #   return( st_drop_geometry(LOCS)[,vnames] )
   }



  if ( lookup_from %in% c("carstm" ) & lookup_to == "points" )  {
     # areal unit (LU) to points (LOCS)
    LU = carstm_model( p=pO, DS="carstm_modelled_summary" ) 
    if (is.null(LU)) stop("Carstm predicted fields not found")
    if (!exists(vnames_from, LU)) {
      message( "vnames_from: ", vnames_from, "not found. You probably want: pcaX.predicted" )
      stop()
    }
    names(LU)[ which(names(LU) == vnames_from ) ] =  vnames

#    nw = length(LU$dyear)
    ny = length(LU$year)
    yr0 = min(as.numeric(LU$year))


    AU = sf::st_transform( LU$sppoly, crs=st_crs(pO$aegis_proj4string_planar_km) )
    AU$au_index = 1:nrow(AU)
    AU = st_cast(AU, "POLYGON")
           
    LOCS = sf::st_as_sf( LOCS, coords=c("lon", "lat") )
    st_crs(LOCS) = st_crs( projection_proj4string("lonlat_wgs84") )
    LOCS = sf::st_transform( LOCS, crs=st_crs(pO$aegis_proj4string_planar_km) )
    LOCS$AUID = st_points_in_polygons( pts=LOCS, polys = AU[, "AUID"], varname= "AUID" )   
    LOCS$AU_index = match( LOCS$AUID, LU$AUID  )    
    
    if (! "POSIXct" %in% class(LOCS$timestamp)  ) LOCS$timestamp =  lubridate::date_decimal( LOCS$timestamp, tz=tz )
    LOCS$yr = lubridate::year(LOCS$timestamp) 
    LOCS$dyear = lubridate::decimal_date( LOCS$timestamp ) - LOCS$yr

    # TIMESTAMP_index = array_map( "ts->2", LOCS [, c("yr", "dyear")], dims=c(ny, nw), res=c( 1, 1/nw ), origin=c( yr0, 0) )
    TIMESTAMP_index = array_map( "ts->year_index", LOCS [, c("yr" )], dims=c(ny ), res=c( 1  ), origin=c( yr0 ) )

  
    LOCS[,vnames] = LU[[vnames_from]][ cbind( LOCS$AU_index, TIMESTAMP_index )]

    return( LOCS[,vnames] )  
  
  }


  if ( lookup_from %in% c("carstm") & lookup_to == "areal_units" )  {
    # areal unit (LU) to areal units (AU/LOCS) 

    # from source data: LU = modelled predictions; AU are associated areal units linked by "AUID" 
    LU = carstm_model( p=pO, DS="carstm_modelled_summary" ) 
    if (is.null(LU)) stop("Carstm predicted fields not found")
    vnames_from = paste( vnames, "predicted", sep=".")

    if (!exists(vnames_from, LU)) {
      message( "vnames_from: ", vnames_from, " not found. You probably want: pcaX.predicted" )
      stop()
    }
    names(LU)[ which(names(LU) == vnames_from ) ] =  vnames

#    nw = length(LU$dyear)
    ny = length(LU$year)
    yr0 = min(as.numeric(LU$year))

    AU = sf::st_transform( LU$sppoly, crs=st_crs(pO$aegis_proj4string_planar_km) )
    AU = st_cast(AU, "POLYGON")
    AU$au_uid = 1:nrow(AU)
    
    # au_uid is internal index of AU /LU
    AU_raster = fasterize::fasterize( AU, raster::raster( AU, res=min(pO$gridparams$res)/2, crs=st_crs( AU ) ), field="au_uid" )  
    AU_pts = sf::st_as_sf( as.data.frame( raster::rasterToPoints(AU_raster)), coords=c("x", "y") )
    st_crs(AU_pts) = st_crs( AU ) 

    pts_AU = match(AU_pts$layer, AU$au_uid[match(LU$AUID, AU$AUID)] ) ## (layer==au_uid) -- to -- LU

    # to target locations: AU_target 
    AU_target = sf::st_transform( AU_target, crs=st_crs(AU) )  # AU_target .... must be sent ... <---------
    AU_target = st_cast( AU_target, "POLYGON")

    if (! "POSIXct" %in% class(LOCS$timestamp)  ) LOCS$timestamp =  lubridate::date_decimal( LOCS$timestamp, tz=tz )
    LOCS$yr = lubridate::year(LOCS$timestamp) 
    # LOCS$dyear = lubridate::decimal_date( LOCS$timestamp ) - LOCS$yr
    
    # TIMESTAMP_index = array_map( "ts->2", LOCS [, c("yr", "dyear")], dims=c(ny, nw), res=c( 1, 1/nw ), origin=c( yr0, 0) )
    TIMESTAMP_index = array_map( "ts->year_index", LOCS [, c("yr" )], dims=c(ny ), res=c( 1  ), origin=c( yr0 ) )

    # id membership in AU_target
    pts_AUID = st_points_in_polygons( pts=AU_pts, polys=AU_target[,"AUID"], varname="AUID" ) 


    for (nn in 1:length(vnames)) {
      vn = vnames[nn]
      LOCS_regridded = apply( 
        LU[[vn]], 
#        MARGIN=c(2,3),
        MARGIN=c(2), 
        FUN=function( LUV ) {
          tapply(X=LUV[ pts_AU ], INDEX=pts_AUID, FUN=FUNC, na.rm=TRUE) 
        }
      )
      space_index = match( LOCS$AUID, as.numeric(as.character(dimnames(LOCS_regridded )[[1]] ))  )   # AUID of AU_target (from pts_AUID)
      LOCS[ , vn ] = LOCS_regridded[ cbind( space_index, TIMESTAMP_index ) ]
    }

    return( LOCS[,vnames] )
 
  } 


}
