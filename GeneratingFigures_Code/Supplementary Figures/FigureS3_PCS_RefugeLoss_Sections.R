library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(colorRamps)
library(maps)
setwd("~/Github/GlobalCoralConnectivity_SpatialRescue_Project/")

#function for plotting Pacific centred maps written by Cole Brookson, used with permission.
plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#load in scores and locations of reef cells
reefdata_full <- read.csv(file = "OriginalDataConversion/GeneratingScoresfrom50ReefsScores/scaled_scores_latlong_reproj_SJ.csv")
#head(reefdata_full)
#dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)] #removing rows that aren't needed + re-organizing row orders
#head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")


#load connectivity matrix
load('OriginalDataConversion/ConvertingWoodetal2014Data/connmat_reduced.RData')


connmat_reduced_list <- list()

#thresholding by 19.5% (average dispersal distance = 10km)
connmat_reduced_list[[1]] <- connmat_reduced
connmat_reduced_list[[1]]@x[connmat_reduced_list[[1]]@x <= 0.195] <- 0
#check
#which(connmat_reduced_list[[2]][which(connmat_reduced_list[[2]] <= 0.195)] > 0) #yes good :)

#thresholding by 9% (average dispersal distance = 20km)
connmat_reduced_list[[2]] <- connmat_reduced
connmat_reduced_list[[2]]@x[connmat_reduced_list[[2]]@x <= 0.09] <- 0

#thresholding by 5.3% (average dispersal distance = 30km)
connmat_reduced_list[[3]] <- connmat_reduced
connmat_reduced_list[[3]]@x[connmat_reduced_list[[3]]@x <= 0.053] <- 0

#thresholding by 3.5% (average dispersal distance = 40km)
connmat_reduced_list[[4]] <- connmat_reduced
connmat_reduced_list[[4]]@x[connmat_reduced_list[[4]]@x <= 0.035] <- 0

#thresholding by 2.6% (average dispersal distance = 50km)
connmat_reduced_list[[5]] <- connmat_reduced
connmat_reduced_list[[5]]@x[connmat_reduced_list[[5]]@x <= 0.026] <- 0

#thresholding by 2% (average dispersal distance = 60km)
connmat_reduced_list[[6]] <- connmat_reduced
connmat_reduced_list[[6]]@x[connmat_reduced_list[[6]]@x <= 0.02] <- 0

#thresholding by 1.6% (average dispersal distance = 70km)
connmat_reduced_list[[7]] <- connmat_reduced
connmat_reduced_list[[7]]@x[connmat_reduced_list[[7]]@x <= 0.016] <- 0

#thresholding by 1.3% (average dispersal distance = 80km)
connmat_reduced_list[[8]] <- connmat_reduced
connmat_reduced_list[[8]]@x[connmat_reduced_list[[8]]@x <= 0.013] <- 0

#thresholding by 1.07% (average dispersal distance = 90km)
connmat_reduced_list[[9]] <- connmat_reduced
connmat_reduced_list[[9]]@x[connmat_reduced_list[[9]]@x <= 0.0107] <- 0

#thresholding by 0.89% (average dispersal distance = 100km)
connmat_reduced_list[[10]] <- connmat_reduced
connmat_reduced_list[[10]]@x[connmat_reduced_list[[10]]@x <= 0.0089] <- 0
#check
#which(connmat_reduced_list[[10]][which(connmat_reduced_list[[10]] <= 0.0089)] > 0) #yes good :)

#thresholding by 0.7526% (average dispersal distance = 110km)
connmat_reduced_list[[11]] <- connmat_reduced
connmat_reduced_list[[11]]@x[connmat_reduced_list[[11]]@x <= 0.007526] <- 0

#thresholding by 0.6451% (average dispersal distance = 120km)
connmat_reduced_list[[12]] <- connmat_reduced
connmat_reduced_list[[12]]@x[connmat_reduced_list[[12]]@x <= 0.006451] <- 0

#thresholding by 0.55% (average dispersal distance = 130km)
connmat_reduced_list[[13]] <- connmat_reduced
connmat_reduced_list[[13]]@x[connmat_reduced_list[[13]]@x <= 0.0055] <- 0

#thresholding by 0.4732% (average dispersal distance = 140km)
connmat_reduced_list[[14]] <- connmat_reduced
connmat_reduced_list[[14]]@x[connmat_reduced_list[[14]]@x <= 0.004732] <- 0

#thresholding by 0.4% (average dispersal distance = 150km)
connmat_reduced_list[[15]] <- connmat_reduced
connmat_reduced_list[[15]]@x[connmat_reduced_list[[15]]@x <= 0.004] <- 0

#thresholding by 0.35% (average dispersal distance = 160km)
connmat_reduced_list[[16]] <- connmat_reduced
connmat_reduced_list[[16]]@x[connmat_reduced_list[[16]]@x <= 0.0035] <- 0

#thresholding by 0.3011% (average dispersal distance = 170km)
connmat_reduced_list[[17]] <- connmat_reduced
connmat_reduced_list[[17]]@x[connmat_reduced_list[[17]]@x <= 0.003011] <- 0

#thresholding by 0.26% (average dispersal distance = 180km)
connmat_reduced_list[[18]] <- connmat_reduced
connmat_reduced_list[[18]]@x[connmat_reduced_list[[18]]@x <= 0.0026] <- 0

#thresholding by 0.2259% (average dispersal distance = 190km)
connmat_reduced_list[[19]] <- connmat_reduced
connmat_reduced_list[[19]]@x[connmat_reduced_list[[19]]@x <= 0.002259] <- 0

#thresholding by 0.2% (average dispersal distance = 200km)
connmat_reduced_list[[20]] <- connmat_reduced
connmat_reduced_list[[20]]@x[connmat_reduced_list[[20]]@x <= 0.002] <- 0

#thresholding by 0.1721% (average dispersal distance = 210km)
connmat_reduced_list[[21]] <- connmat_reduced
connmat_reduced_list[[21]]@x[connmat_reduced_list[[21]]@x <= 0.001721] <- 0

#thresholding by 0.1506% (average dispersal distance = 220km)
connmat_reduced_list[[22]] <- connmat_reduced
connmat_reduced_list[[22]]@x[connmat_reduced_list[[22]]@x <= 0.001506] <- 0

#thresholding by 0.1299% (average dispersal distance = 230km)
connmat_reduced_list[[23]] <- connmat_reduced
connmat_reduced_list[[23]]@x[connmat_reduced_list[[23]]@x <= 0.001299] <- 0

#thresholding by 0.1182% (average dispersal distance = 240km)
connmat_reduced_list[[24]] <- connmat_reduced
connmat_reduced_list[[24]]@x[connmat_reduced_list[[24]]@x <= 0.001182] <- 0

#thresholding by 0.1% (average dispersal distance = 250km)
connmat_reduced_list[[25]] <- connmat_reduced
connmat_reduced_list[[25]]@x[connmat_reduced_list[[25]]@x <= 0.001] <- 0

#thresholding by 0.09% (average dispersal distance = 260km)
connmat_reduced_list[[26]] <- connmat_reduced
connmat_reduced_list[[26]]@x[connmat_reduced_list[[26]]@x <= 0.0009] <- 0

#thresholding by 3.3% (average dispersal distance = 42km)
connmat_reduced_list[[27]] <- connmat_reduced
connmat_reduced_list[[27]]@x[connmat_reduced_list[[27]]@x <= 0.033] <- 0
#what is the distribution of the scores?
#hist(reefdata$score) #centred on 0, -1 to 0.2

#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre

world <- readOGR(dsn = "GeneratingFigures_Code/AdditionalFiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)



plateaupoint_pr_list <- list()
numreseeded_pr_list <- list()
plateaupoint_cr_list <- list()
numreseeded_cr_list <- list()
numclust_pr <- rep(NA,27)
gmmean_pr <- rep(NA,27)
largest_clust_pr <- rep(NA,27)
avgss_pr <- rep(NA,27)
avgnd_pr <- rep(NA,27)
numclust_cr <- rep(NA,27)
gmmean_cr <- rep(NA,27)
largest_clust_cr <- rep(NA,27)
avgss_cr <- rep(NA,27)
avgnd_cr <- rep(NA,27)

for(j in 1:27){

print(paste("j =", j))
    
connmat_reduced <- connmat_reduced_list[[j]]
avgdist <- j*10
if(j==27){
  avgdist <- 42
}

print("starting PCS scenario")
#RESEEDING FROM PCS SCENARIO
ord <- order(reefdata$score) #lowest to highest
numsteps <- 50
num_remove <- 434
num_iter <- 20
final_reef <- (num_iter*num_remove)
reefdata$num <- seq(1,12292,by=1)
starterreefs <- reefdata$num[-ord[1:final_reef]]
reefdata$num <- NULL
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = (numsteps+1)) #50 time steps
abundance_mat[,1][starterreefs] <- 100 #time step 0: assign all the starterreefs to an abundance of 100

#calculate abundance after 1 time step
abundance_mat[,2] <- as.vector(t(abundance_mat[,1]) %*% connmat_reduced) #need to do this because it thinks of the right side as a matrix
#re-set abundance of all non-vacant reefs to 100
dummy <- abundance_mat[,2]
dummy[dummy > 0] <- 100 
#calculate abundance for the next time steps - note: the total abundance/sum of the subsequent columns should decrease every time step
for(i in 3:(numsteps+1)){
  abundance_mat[,i] <- as.vector(t(dummy) %*% connmat_reduced)
  dummy <- abundance_mat[,i]
  dummy[dummy > 0] <- 100
}


#Plot the number that were re-seeded over time
stop <- "no"
numreseeded_pr <- rep(NA,(numsteps+1))
numreseeded_pr[1] <- 0
for(i in 2:(numsteps+1)){
  numreseeded_pr[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[starterreefs,i][abundance_mat[starterreefs,i]>0])	
  if((numreseeded_pr[i] <= numreseeded_pr[i-1]) & stop == "no"){ #added 'or' condition to be consistent with the other scenarios
    plateaupoint <- i-1
    stop <- "yes"
  }
}
plateaupoint_pr_list[[j]] <- plateaupoint
numreseeded_pr_list[[j]] <- numreseeded_pr
print("finished calculating abundance_mat")



id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)


print("final connectivity metric calculating")
#final connectivity metrics
connmat_PCS_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_PCS <- graph.adjacency(as.matrix(connmat_PCS_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]

#number of networks
cl_PCS <- clusters(g_PCS,mode="weak") #number of networks, size of each network, membership of each network
numclust_pr[j] <- cl_PCS$no #99
#geom_mean of network size
gmmean_pr[j] <- gm_mean(cl_PCS$csize) #8.917658
#largest network
largest_clust_pr[j] <- max(cl_PCS$csize) #1852
#average source strength
so_st_PCS <- ego_size(graph = g_PCS, order = 40, nodes = V(g_PCS), mode = "out")
avgss_pr[j] <- mean(so_st_PCS) #114.0211
#average node degree
dg_PCS <- degree(g_PCS) #getting node degree
avgnd_pr[j] <- mean(dg_PCS) #10.22306

abundance_mat <- NA
plateaupoint <- 500

print("starting refuge-loss scenario")
#RESEEDING FROM REFUGE-LOSS SCENARIO
ord <- order(reefdata$score, decreasing = TRUE) #highest to lowest
numsteps <- 50
num_remove <- 434
num_iter <- 20
final_reef <- (num_iter*num_remove)
reefdata$num <- seq(1,12292,by=1)
starterreefs <- reefdata$num[-ord[1:final_reef]]
reefdata$num <- NULL
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = (numsteps+1)) #10 time steps
abundance_mat[,1][starterreefs] <- 100 #time step 0: assign all the starterreefs to an abundance of 100

#calculate abundance after 1 time step
abundance_mat[,2] <- as.vector(t(abundance_mat[,1]) %*% connmat_reduced) #need to do this because it thinks of the right side as a matrix
#re-set abundance of all non-vacant reefs to 100
dummy <- abundance_mat[,2]
dummy[dummy > 0] <- 100 
#calculate abundance for the next 9 time steps - note: the total abundance/sum of the subsequent columns should decrease every time step
for(i in 3:(numsteps+1)){
  abundance_mat[,i] <- as.vector(t(dummy) %*% connmat_reduced)
  dummy <- abundance_mat[,i]
  dummy[dummy > 0] <- 100
}


stop <- "no"
numreseeded_cr <- rep(NA,(numsteps+1))
numreseeded_cr[1] <- 0
for(i in 2:(numsteps+1)){
  numreseeded_cr[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[starterreefs,i][abundance_mat[starterreefs,i]>0])	
  if((numreseeded_cr[i] <= numreseeded_cr[i-1]) & stop == "no"){ #added in the 'or previous = bigger' condition bc of '0' issues
    plateaupoint <- i-1
    stop <- "yes"
  }
}
plateaupoint_cr_list[[j]] <- plateaupoint
numreseeded_cr_list[[j]] <- numreseeded_cr
print("finished calculating abundance matrix")


id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)



print("starting final connectivity metrics")
connmat_RefugeLoss_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_RefugeLoss <- graph.adjacency(as.matrix(connmat_RefugeLoss_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]


#number of networks
cl_RefugeLoss <- clusters(g_RefugeLoss,mode="weak") #number of networks, size of each network, membership of each network
numclust_cr[j] <- cl_RefugeLoss$no #175
#geom_mean of network size
gmmean_cr[j] <- gm_mean(cl_RefugeLoss$csize) #5.803611
#largest cluster
largest_clust_cr[j] <- max(cl_RefugeLoss$csize) #1038
#average source strength
so_st_RefugeLoss <- ego_size(graph = g_RefugeLoss, order = 40, nodes = V(g_RefugeLoss), mode = "out")
avgss_cr[j] <- mean(so_st_RefugeLoss) #115.4481
#average node degree
dg_RefugeLoss <- degree(g_RefugeLoss) #getting node degree
avgnd_cr[j] <- mean(dg_RefugeLoss) #9.710095

abundance_mat <- NA

initialreefs = 12292 - final_reef


}

save(plateaupoint_pr_list, numreseeded_pr_list, plateaupoint_cr_list, numreseeded_cr_list, numclust_pr, gmmean_pr, largest_clust_pr, avgss_pr, avgnd_pr, numclust_cr, gmmean_cr, largest_clust_cr, avgss_cr, avgnd_cr, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/AllReseedingresults_norandom.RData")


