

    # -----------------------------
    # ordination
    year.assessment = 2021

    require(aegis.speciescomposition)
   
    p = speciescomposition_parameters( yrs=1970:year.assessment )


    spatial_domain = "SSE"
    yrs = 1970:2021


  speciescomposition_db( DS="speciescomposition.ordination.redo", p=p )  # analsysis

  if (0) { 
    # extract summaries and plot
      pca = speciescomposition_db( DS="pca", p=p )  # analsysis
      ca  = speciescomposition_db( DS="ca", p=p )  # analsysis

      toplot = as.data.frame( pca$cscores )
      toplot$vern = taxonomy.recode( from="spec", to="taxa", tolookup=rownames( toplot ) )$vern

      plot(V2 ~ V1, toplot, type="n")
      text( V2 ~ V1, labels=vern, data=toplot )


    plot(V3 ~ V1, toplot, type="n")
    text( V3 ~ V1, labels=vern, data=toplot )

  }

  speciescomposition_db( DS="speciescomposition.redo", p=p ) # compute planar coords and remove dups

   
 