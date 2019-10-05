

  # -----------------------------
  # ordination
  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date())
    year.assessment=lubridate::year(Sys.Date()) -1
  }

  p = aegis.speciescomposition::speciescomposition_parameters( yrs=1999:year.assessment )

  speciescomposition.db( DS="speciescomposition.ordination.redo", p=p )  # analsysis
  pca = speciescomposition.db( DS="pca", p=p )  # analsysis
  ca  = speciescomposition.db( DS="ca", p=p )  # analsysis

  speciescomposition.db( DS="speciescomposition.redo", p=p ) # compute planar coords and remove dups

