library(rgdal)
library(Matrix)
library(igraph)
library(fields)

#load in the shape file with all of the data
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/cells12292.RData') #need to load rgdal library first (maybe)
#this ^ file has 'b' (shape file) and 'id.noReef' which are the rows that were removed for having no reef SA in them
head(b@data) #removed the rows with no reef
reefdata <- b@data

#load in Sally's connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Conn_Mat_Sum_sparse.RData')
head(Conn_Mat_Sum_sparse.RData)
connmat <- Conn_Mat_Sum_sparse

#remove the rows and columns that have no reef SA
range(reefdata$PolyNo) #1 to 12397
#recall: id.noReef <- which(is.na(a@data$ReefArea))
connmat_reduced <- connmat[-id.noReef,-id.noReef]
dim(connmat_reduced) #12292 12292
save(connmat_reduced, file = "connmat_reduced.RData")