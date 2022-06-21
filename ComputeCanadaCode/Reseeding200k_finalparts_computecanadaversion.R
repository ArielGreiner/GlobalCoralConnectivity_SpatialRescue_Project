#library(rgdal) #hopefully don't need
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(ggplot2)
#library(colorRamps) #don't think i need
#library(maps) #hopefully don't need
print("loaded libraries")

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

setwd("/home/agreiner/scratch/worldconn_randomrep")
#setwd("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/R_amalgamations/")

####3.3% Removal Scenario####
#load connectivity matrix
load('/home/agreiner/scratch/worldconn_randomrep/connmat_reduced.RData')
#load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA
print(paste0("connmat_reduced first row, first element = ", connmat_reduced[1,1]))

#thresholding by 3.3% probability (see 'Comparing_Options' file, note: connmat_reduced@x turns the matrix into a vector)
connmat_reduced@x[connmat_reduced@x <= 0.033] <- 0
#check - okay it's actually fine
print("reduced conn mat")
#which(connmat_reduced[which(connmat_reduced <= 0.033)] > 0)

#load in the centroids, NOTE: NEED TO FIX THE LONGITUDE IN THESE ONES (>180 ones need to have -360 subtracted)
load('/home/agreiner/scratch/worldconn_randomrep/centroids.RData')
#load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/centroids.RData')
#head(centroids) 

load('/home/agreiner/scratch/worldconn_randomrep/cells12292.RData')
#load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/cells12292.RData') #this includes 'id.noReef' and 'b' (shapefile), GOING TO IGNORE B

reefdata_full <- read.csv(file = "/home/agreiner/scratch/worldconn_randomrep/scaled_scores_latlong_reproj_SJ.csv")
#reefdata_full <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/RecalculatingResultswithNewScores_2.2020/scaled_scores_latlong_reproj_SJ.csv")
head(reefdata_full)
dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)]
head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")
reefdata2 <- reefdata[,c(2,3,4)]
print("loaded and edited reefdata and reefdata2")

print(paste0("centroids first element = ", centroids[1,1], "reefdata first element = ", reefdata[1,1]))

print("remove the grid cells with no reef SA")
#remove the rows corresponding to the grid cells that have no reef SA
centroids_reduced <- centroids[-id.noReef,]

#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360
centroids_reduced$Longitude[centroids_reduced$Longitude > 180] <- centroids_reduced$Longitude[centroids_reduced$Longitude > 180] - 360

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre
#reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] <- reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] - 360
print("edited longitude, fixed centroids")
#probs don't need the centroids here but ah well

#need to load in id.totalnotseeded
load("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/reseedingoutputset1.RData")

random_iter <- 1000
numsteps <- 50
#NEW 7.26.2020
#final connectivity metrics
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#connmat_Random_final <- g_Random <- cl_Random <- so_st_Random <- dg_Random <- list() #don't have enough memory to save these
meanclust <- maxclust <- avgsost <- avgdeg <- numclust <- rep(NA,(random_iter*200))
print("entering final connectivity metrics for-loop")
j = 11
for(i in 50001:(random_iter*200)){
  print(paste("i = ",i))
  #connmat_Random_final[[i]] <- connmat_reduced[-id.totalnotseeded[[i]], -id.totalnotseeded[[i]]]
  connmat_Random_final <- connmat_reduced[-id.totalnotseeded[[i]], -id.totalnotseeded[[i]]]
  #g_Random[[i]] <- graph.adjacency(as.matrix(connmat_Random_final[[i]]), weighted = TRUE)
  g_Random <- graph.adjacency(as.matrix(connmat_Random_final), weighted = TRUE) 
  #connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
  #df[[i]] <- reefdata1[-ord[1:final_reef],]
  
  #number of clusters
  #cl_Random[[i]] <- clusters(g_Random[[i]],mode="weak") #number of clusters, size of each cluster, membership of each cluster
  cl_Random <- clusters(g_Random,mode="weak") #number of clusters, size of each cluster, membership of each cluster
  #numclust[i] <- cl_Random[[i]]$no ...all of the instances below with cl_Random were once cl_Random[[i]]
  numclust[i] <- cl_Random$no
  #geom_mean of cluster size
  meanclust[i] <- gm_mean(cl_Random$csize) 
  #largest cluster
  maxclust[i] <- max(cl_Random$csize) 
  #average source strength
  #so_st_Random[[i]] <- ego_size(graph = g_Random[[i]], order = 40, nodes = V(g_Random[[i]]), mode = "out")
  so_st_Random <- ego_size(graph = g_Random, order = 40, nodes = V(g_Random), mode = "out")
  #avgsost[i] <- mean(so_st_Random[[i]])
  avgsost[i] <- mean(so_st_Random)
  #average node degree
  #dg_Random[[i]] <- degree(g_Random[[i]]) #getting node degree
  dg_Random <- degree(g_Random) #getting node degree
  avgdeg[i] <- mean(dg_Random, na.rm = T) #WHY AM I GETTING NAs? #12.16.2020...it was dg_Random[i] not dg_Random[[i]] and maybe that was the problem?
  if(i == 5000*j){
    save(numclust,meanclust,maxclust,avgsost,avgdeg, file = paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/meanfinalconnectivitymetrics_all_",5000*j,".RData"))
    j = j + 1
  }
}
mean_numclust <- mean(numclust) 
mean_gmeanclust <- mean(meanclust) 
mean_maxclust <- mean(maxclust)
mean_avgsost <- mean(avgsost) 
mean_avgdeg <- mean(avgdeg, na.rm = T) #this one was causing NaN issues before
save(mean_numclust,mean_gmeanclust,mean_maxclust,mean_avgsost, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/meanfinalconnectivitymetrics_noavgdeg_done.RData")
save(mean_numclust,mean_gmeanclust,mean_maxclust,mean_avgsost,mean_avgdeg, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/meanfinalconnectivitymetrics_all_done.RData")
