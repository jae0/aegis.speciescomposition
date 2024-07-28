
survey_data_prepare = function(p,  cthreshold = 0.005, ci=0.999){
    require(data.table)

    # catch info by species and set ; NOTE: id = paste( x$trip, x$set, sep="." )
    sc = survey_db( p=p, DS="cat.filter"  ) 
    sc = sc[ , c("id", "totno_adjusted", "totwgt_adjusted",  "spec_bio" )]
    sc_species = unique(  sc$id[ taxonomy.filter.taxa( sc$spec_bio, method=p$taxa, outtype="internalcodes" ) ] )

    # trip/set loc information by set
    set = survey_db( p=p, DS="filter"  ) 
    set = set[ which( set$id %in% sc_species), ]
    set = set[ , c("id", "yr", "dyear", "sa", "sa_towdistance", "lon", "lat", "z", "t",
      "timestamp", "gear", "vessel", "data.source" )]

    sc = merge(sc, set, by="id", all.x=TRUE, all.y=FALSE) 
      
    # NOTE:: the US 4 seam behaves very differently and has a shorter standard tow .. treat as a separate data source 
    # of course it misses out on the high abundance period from the 1980s but no model-based solution possible at this point 
    # due simply to CPU speed issues

    sc$qscore = NA 
    sc$density = sc$totwgt_adjusted / sc$sa_towdistance  ### <<< --- using towed distance due to inconsistant sa for groundfish from RV surveys


    if (0) {
      qi = quantile(sc$density, ci, na.rm=TRUE)
      ii = which( sc$density > qi )
      hist( log(sc$density[ii]) )
    }

    # NOTE:: non-zero catches are not recorded in cat, also 
    # distribution of density is very long tailed .. 
    # use quantile -> normal transformatiom to fit 
    # the 99.9% CI (ci=0.999 probability) to be between (0, 1)
    gears = unique(sc$gear)
    taxa = unique(sc$spec_bio)
    for ( g in gears ) {
      for ( tx in taxa ) {
        ii = which( sc$gear==g & sc$spec_bio== tx & sc$density > 0 )
        if (length(ii) > 30 ) {
          qnts = quantile_estimate( sc$density[ii] )
          sc$qscore[ii] = quantile_to_normal( qnts, mean=0.5, ci=ci ) # convert to quantiles, by survey to have prob=ci to fit between 0,1
        }
      }
    }
 
    m = data.table::dcast( setDT(sc), 
      formula =  id ~ spec_bio, value.var="qscore", 
      fun.aggregate=mean, fill=NA, drop=FALSE, na.rm=TRUE
    )  # mean is just to keep dcast happy

    id = m$id 
    m$id = NULL
    sps = colnames(m)

    m = as.matrix(m[]) 
    dimnames(m)[[1]] = id
    # remove low counts (absence) in the timeseries  .. species (cols) only
    # cthreshold = 0.005  # 0.01%  -> 97;  0.05 => 44 species; 0.001 => 228 species; 0.005 -> 146 
    m [ m < cthreshold ] = 0
    m [ !is.finite(m) ] = 0

    # reduce:: remove taxa/locations with very low counts (~cthreshold)
    finished = FALSE
    while( !(finished) ) {
      oo = colSums( ifelse( is.finite(m), 1, 0 ), na.rm=TRUE)
      i = unique( which(rowSums(m, na.rm=TRUE) == 0 ) )
      uu = colSums(m, na.rm=TRUE)/oo
      j = unique( which(!is.finite( uu ) | uu < cthreshold ) )
      if ( ( length(i) == 0 ) & (length(j) == 0 ) ) finished=TRUE
      if (length(i) > 0 ) m = m[ -i , ]
      if (length(j) > 0 ) m = m[ , -j ]
    }

    return( list(set=set, sc=sc, m=m ) )
}
