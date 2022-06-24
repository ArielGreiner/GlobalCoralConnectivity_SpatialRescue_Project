# Calculate the spatial distance (Euclidean) between connected reefs
# Eventually, we could replace it with the marine distance
# But that calculation is long and computationally intensive: it's a least-cost algorithm;
# in addition, little bias is expected if we screen at a distance as short as 40 km

# We do not calculate the distance among all reefs because with >12000 reefs the distance matrix would be huge

library(Matrix)
library(sp)

setwd("~/Github/GlobalCoralConnectivity_SpatialRescue_Project/")

#load connectivity matrix
load('OriginalDataConversion/ConvertingWoodetal2014Data/connmat_reduced.RData')
a <- connmat_reduced
rm(connmat_reduced)
dim(a)
# 12292 reefs, i.e. the ones that have a reefsurface > 0

nonzero <- summary(a) # Gives a dataframe with the indices (row number and column number) of the nonzero entries and the value of the entry
head(nonzero)
# These indices counts row starting from 1, so it is convenient for us
# Instead, the internal representation of the sparse matrix objects count row starting from 0. 
# This will be important later, when converting values to sparse matrix


nonzero <- nonzero[,c(1:2)] # We do not need the values so we get rid of them
# A bunch of statistics:
    num.nonzero <- dim(nonzero)[1] # There are 1,523,284 nonzero connections
    # The theoretical number of possible connection is 12922 * 12922 = 151,093,264, so the eonnectance is:
    num.nonzero / (12292^2)
    # Still, these connections involve all reefs:
    length(unique(nonzero[,1])) # They involve all the 12292 reefs as starting point
    length(unique(nonzero[,2])) # They involve all the 12292 reefs as destination point

# now we need to calculate the spatial Euclidean distance corresponding to these 1,523,284 nonzero connections
load('GeneratingFigures_Code/AdditionalFiles/centroids_reduced.RData') #if you are interested in generating this file, contact me.
euclid.dist <- rep(NA,num.nonzero)
# The following for loop takes about 5 minutes on a decent computer
for (k in 1 : num.nonzero ) {
    if (k %% 10000 == 0) {
        cat(k,"of",num.nonzero,"\n")
        flush.console()
    }
    origin_reef_coord <- as.matrix(centroids_reduced[nonzero$i[k],])    # the first argument of spDistsN1 must be a matric
    destination_reef_coord <- as.numeric(centroids_reduced[nonzero$j[k],]) # the second argument of spDistsN1 must be numeric
    euclid.dist[k] <- spDistsN1(origin_reef_coord, destination_reef_coord, longlat=T) # Euclidean distance in km
}

euclid.dist.orig <- euclid.dist
euclid.dist <- round(euclid.dist)

# Now we create a new sparse matrix using the indices of the nonzero entries and the calculated euclidean distances
euclid.dist.mat.T <- new("dgTMatrix",
                        i = as.integer(nonzero$i - 1),          # We subtract 1 because the internal representation of this class counts rows and columns starting from 0
                        j = as.integer(nonzero$j - 1),
                        x = euclid.dist,
                        Dim = dim(a) )              # Same dimensions as the connectivity matrix a

# Convert it into the dgCMatrix class (more efficient)
euclid.dist.mat <- as(euclid.dist.mat.T, "dgCMatrix")

# Save it
save(euclid.dist.mat,file='GeneratingFigures_Code/AdditionalFiles/Euclid.dist.mat.12292.RData')

# Some statistics
load('GeneratingFigures_Code/AdditionalFiles/Euclid.dist.mat.12292.RData') 
hist(euclid.dist.mat@x,breaks=seq(0,8000,40))
length(which(euclid.dist.mat@x>40)) / length(euclid.dist.mat@x)  # Fraction of connections that are > 40 km
