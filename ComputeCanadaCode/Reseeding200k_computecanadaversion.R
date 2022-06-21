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

#RESEEDING FROM RANDOM SCENARIO
#set.seed(2)
abundmatlist <- starterreefslist <- numreseededrrlist <- id.totalnotseeded <- perreseededrr <- list()
numsteps <- 50
num_remove <- 434
num_iter <- 20
random_iter <- 1000
final_reef <- (num_iter*num_remove)
plateaupoint <- 500
plateaupoint_rr <- rep(NA,(random_iter*200))
reefdata$num <- seq(1,12292,by=1)
print("entering re-seeding for loop")
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/ordr/RandomScenario_3.3_ordr_1kreps_",k,".RData"))
  for(j in 1:random_iter){
    print(paste("k = ",k, "j = ",j))
    #ord <- sample(seq(1,length(reefdata$score),by=1))
    starterreefs <- reefdata$num[-ordr[[j]][1:final_reef]]
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
    #abundmatlist[[j]] <- abundance_mat
    #starterreefslist[[j]] <- starterreefs
    
    #REMOVING THIS BECAUSE NOT ENOUGH MEMORY TO HOLD THEM ON COMPUTECANADA
    #abundmatlist[[((k-1)*1000 + j)]] <- abundance_mat
    #REMOVING THIS BECAUSE DONT NEED TO SAVE IT
    #starterreefslist[[((k-1)*1000 + j)]] <- starterreefs
    
    
    
    #Plot the number that were re-seeded over time
    stop <- "no"
    numreseeded_rr <- rep(NA,(numsteps+1))
    numreseeded_rr[1] <- 0
    for(i in 2:(numsteps+1)){
      numreseeded_rr[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[starterreefs,i][abundance_mat[starterreefs,i]>0])	
      if((numreseeded_rr[i] <= numreseeded_rr[i-1]) & stop == "no"){ #OR if more reseeded in the previous step (to deal with the flip-flopping)
        plateaupoint <- i-1
        stop <- "yes"
      }
    }
    plateaupoint_rr[((k-1)*1000 + j)] <- plateaupoint
    numreseededrrlist[[((k-1)*1000 + j)]] <- numreseeded_rr
    
    id.totalnotseeded[[((k-1)*1000 + j)]] <- which(abundance_mat[,(numsteps+1)] == 0)
  }
}
meanplateauptrr <- mean(plateaupoint_rr)
save(plateaupoint_rr, meanplateauptrr, numreseededrrlist,id.totalnotseeded, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/reseedingoutputset1.RData")
save(plateaupoint_rr, meanplateauptrr, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/plateaupoint_all.RData")
save(numreseededrrlist,id.totalnotseeded, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/numreseeded.RData")

#save(abundmatlist, starterreefslist, plateaupoint_rr, meanplateauptrr, numreseededrrlist,id.totalnotseeded, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/reseedingoutputset1.RData")
#save(plateaupoint_rr, meanplateauptrr, numreseededrrlist,id.totalnotseeded, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/reseedingoutputset1_abr.RData")
#save(abundmatlist, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/abundmatlist.RData")
#save(starterreefslist, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/starterreefslist.RData")
#save(plateaupoint_rr, meanplateauptrr, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/plateaupoint_all.RData")
#save(numreseededrrlist,id.totalnotseeded, file  = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/numreseeded.RData")

print("finished re-seeding for loop")
print(paste("some sample outputs, meanplateauptrr = ", meanplateauptrr, "id.totalnotseeded[[1]] =",id.totalnotseeded[[1]]))

print("entering numreseeded for-loop")
numreseededrr_mean <- rep(NA,(numsteps+1))
numreseededrr_sd <- rep(NA,(numsteps+1))
set <- rep(NA,(random_iter*200))
for(i in 1:(numsteps+1)){ #can probably change this to be a loop over j in 1:(random_iter*200)
  for(k in 1:200){
    print(paste("i = ",i,"k = ",k))
    for(j in 1:random_iter){
      set[((k-1)*1000 + j)] <- numreseededrrlist[[((k-1)*1000 + j)]][i]
    }
  }
  numreseededrr_mean[i] <- mean(set)
  numreseededrr_sd[i] <- sd(set)
  set <- rep(NA,(random_iter*200))	
}
save(numreseededrr_mean,numreseededrr_sd, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/numreseededrr_all.RData")
print("finished numreseeded for-loop")
print(paste("sample output: numreseededrr_mean[1]=", numreseededrr_mean, "numreseededrr_sd[1] =", numreseededrr_sd[1]))

png(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/Realized_CU_Random_200k_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseededrr_mean,main="Realized Random Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs", type = "n")
for(i in 1:(random_iter*200)){
  lines(x = seq(1,(numsteps+1),by=1),y = numreseededrrlist[[i]],col = alpha("blue", 0.2))
}
points(x=seq(1,(numsteps+1),by=1),y=numreseededrr_mean, pch=20)
dev.off()

#NEW 7.26.2020
#final connectivity metrics
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
connmat_Random_final <- g_Random <- cl_Random <- so_st_Random <- dg_Random <- list()
meanclust <- maxclust <- avgsost <- avgdeg <- numclust <- rep(NA,(random_iter*200))
print("entering final connectivity metrics for-loop")
for(i in 1:(random_iter*200)){
  print(paste("i = ",i))
  connmat_Random_final[[i]] <- connmat_reduced[-id.totalnotseeded[[i]], -id.totalnotseeded[[i]]]
  g_Random[[i]] <- graph.adjacency(as.matrix(connmat_Random_final[[i]]), weighted = TRUE) 
  #connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
  #df[[i]] <- reefdata1[-ord[1:final_reef],]
  
  #number of clusters
  cl_Random[[i]] <- clusters(g_Random[[i]],mode="weak") #number of clusters, size of each cluster, membership of each cluster
  numclust[i] <- cl_Random[[i]]$no 
  #geom_mean of cluster size
  meanclust[i] <- gm_mean(cl_Random[[i]]$csize) 
  #largest cluster
  maxclust[i] <- max(cl_Random[[i]]$csize) 
  #average source strength
  so_st_Random[[i]] <- ego_size(graph = g_Random[[i]], order = 40, nodes = V(g_Random[[i]]), mode = "out")
  avgsost[i] <- mean(so_st_Random[[i]]) 
  #average node degree
  dg_Random[[i]] <- degree(g_Random[[i]]) #getting node degree
  avgdeg[i] <- mean(dg_Random[i], na.rm = T) #WHY AM I GETTING NAs?
}
mean_numclust <- mean(numclust) 
mean_gmeanclust <- mean(meanclust) 
mean_maxclust <- mean(maxclust)
mean_avgsost <- mean(avgsost) 
mean_avgdeg <- mean(avgdeg, na.rm = T) #this one was causing NaN issues before
save(mean_numclust,mean_gmeanclust,mean_maxclust,mean_avgsost, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/meanfinalconnectivitymetrics_noavgdeg.RData")
save(mean_numclust,mean_gmeanclust,mean_maxclust,mean_avgsost,mean_avgdeg, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/meanfinalconnectivitymetrics_all.RData")

#NOT SURE IF THIS MAKES SENSE TO PLOT...NEED TO THINK ABOUT
#MOVED BELOW THE OTHER STUFF BECAUSE I CARE LESS ABOUT THIS RUNNING
#Histogram of percentage of each of the initial clusters that have been rebuilt
#load in original cluster designations
#reefdata2 <- reefdata[,c(2,3,4)]
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
#range(clusters(g_orig, mode = "weak")$membership)[2] #604 clusters

#reefdata is in the same order as abundance_mat (see: use + construction of id.reefs)
#so can add in another column to reefdata2 with abundance y/n re: re-seeding
#reefdata2 <- reefdata[,c(2,3,4)]
print("entering the optional section, per-cluster dealios")
print("entering 1st per-reseeded for-loop")
for(j in 1:(random_iter*200)){
  print(paste("j = ",j))
  reefdata2$final_abund_rr <- 1
  reefdata2$final_abund_rr[id.totalnotseeded[[j]]] <- 0
  head(reefdata2$final_abund_rr)
  
  #reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)
  for(i in 1:(clusters(g_orig, mode = "weak")$no)){
    reseededclusters$num_reseeded_rr[i] <- sum(reefdata2$final_abund_rr[reefdata2$cl.m == i])
    reseededclusters$per_reseeded_rr[i] <- (reseededclusters$num_reseeded_rr[i]/reseededclusters$clustersize[i])*100
  }
  perreseededrr[[j]] <- reseededclusters$per_reseeded_rr
}
print("ending of first per-reseeded for loop")

print("entering 2nd per-reseeded for-loop")
perreseededrr_mean <- rep(NA,(numsteps+1))
perreseededrr_sd <- rep(NA,(numsteps+1))
set <- rep(NA,(random_iter*200))
for(i in 1:(numsteps+1)){
  print(paste("i =", i))
  for(j in 1:(random_iter*200)){
    set[j] <- perreseededrr[[j]][i]
  }
  perreseededrr_mean[i] <- mean(set)
  perreseededrr_sd[i] <- sd(set)
  set <- rep(NA,(random_iter*200))	
}
save(perreseededrr, perreseededrr_mean, perreseededrr_sd, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/per-reseededmetrics_all.RData")
print("ended both per-reseeded for-loops")
print(paste("some sample outputs: perreseededrr_mean[1] = ", perreseededrr_mean[1], "perreseededrr_sd[1] =", perreseededrr_sd[1]))

png(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/reseeding/Realized_CU_MeanRandom_200k_HistPerReseeded_",numsteps,"steps.png"))
hist(perreseededrr_mean, main = "Mean Percentage of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

#abundance_mat <- NA
#plateaupoint <- 500

