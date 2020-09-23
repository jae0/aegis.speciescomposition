

  # -----------------------------
  # ordination
  year.assessment = 2020

  p = aegis.speciescomposition::speciescomposition_parameters( yrs=1999:year.assessment )

  speciescomposition_db( DS="speciescomposition.ordination.redo", p=p )  # analsysis
  pca = speciescomposition_db( DS="pca", p=p )  # analsysis
  ca  = speciescomposition_db( DS="ca", p=p )  # analsysis

  speciescomposition_db( DS="speciescomposition.redo", p=p ) # compute planar coords and remove dups

