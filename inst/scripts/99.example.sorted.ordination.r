
project.library( "aegis.base", "aegis" )


# load data of interest (annual time-series)
# .. an example dataset from the Eastern Scotian Shelf
data( "ess.data" )

# choose subset of data of interest and scale/center
rownames( ess.data ) = ess.data$yr  # the rownames must be set
ess.data$yr = NULL                  # remove year as we do not want to use it as a variable
ess.data.scaled = scale( ess.data, center=T, scale=T )

# final call for analysis
ess.pca.stats = pca.analyse.data( ess.data.scaled, plotdata=FALSE, scaledata=FALSE ) # look in work directory for figures and data outputs
str( ess.pca.stats ) # look at contents of the variable
