library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(colorRamps)
library(maps)

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

###REALIZED CONNECTIVITY###
#load in the shape file with all of the data
reefdata_full <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/RecalculatingResultswithNewScores_2.2020/scaled_scores_latlong_reproj_SJ.csv")
head(reefdata_full)
dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)]
head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")

load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/cells12292.RData') #this includes 'id.noReef' and 'b' (shapefile)
reefdata_old <- b@data
head(reefdata_old)

#load connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA

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

#for plotting pacific centred maps
reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
#range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre
#reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] <- reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] - 360


world <- readOGR(dsn = "~/Dropbox/University of Toronto/Research Related/Code from Marco 6.2018 mediterranean larval connectivity/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

#note: plateaupoint_cr = 26, plateaupoint_50r = 23, plateaupoint_pr = 23 so final plateaupoint = 23
ultimateplateau <- 26
reefdata2 <- reefdata[,c(2,3,4)]
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)
#^ not sure what i need of this section above

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

for(j in 15:27){

print(paste("j =", j))
    
connmat_reduced <- connmat_reduced_list[[j]]
avgdist <- j*10
if(j==27){
  avgdist <- 42
}

print("starting predicted scenario")
#RESEEDING FROM PREDICTED SCENARIO
setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/PCS/Reseeding"))
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

png(paste0("Avg",avgdist,"km_CU_Predicted_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseeded_pr,main="Predicted Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs")
dev.off()



id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)
#heatcols <- diverge_hcl(10)
#heatcols <- heat_hcl(12,h=c(0,-100), l = c(75,40), c = c(40,80), power = 1)
#heatcols <- heat_hcl(12,h=c(0,-100),l=c(75,40),c=c(40,80), power=1)
heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue


#pacific centred versions
png(paste0("PCentre_Avg",avgdist,"km_CU_Predicted_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Predicted Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
  points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=20,cex=0.2)
  if(i <= plateaupoint){
    points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=20,cex=0.2)	
  }
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

png(paste0("PCentre_Avg",avgdist,"km_CU_Predicted_Reseed_BRB_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Predicted Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="light blue",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
  points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)
  if(i <= plateaupoint){
    points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)	
  }
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

png(paste0("PCentre_Avg",avgdist,"km_Predicted_allblack.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Predicted Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()



png(paste0("PCentre_Avg",avgdist,"km_CU_Predicted_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Predicted Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=15,cex=0.03) #cex=0.05
for(i in 2:(numsteps+1)){
  points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=15,cex=0.03)
  if(i <= plateaupoint){
    points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=15,cex=0.03)	
  }
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=15,cex=0.03) #not great but maybe better?
dev.off()

#Histogram of percentage of each of the initial clusters that have been rebuilt
#load in original cluster designations
#reefdata2 <- reefdata[,c(2,3,4)]
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
#range(clusters(g_orig, mode = "weak")$membership)[2] #604 clusters

#reefdata is in the same order as abundance_mat (see: use + construction of id.reefs)
#so can add in another column to reefdata2 with abundance y/n re: re-seeding
#reefdata2 <- reefdata[,c(2,3,4)]
reefdata2$final_abund_pr <- 1
reefdata2$final_abund_pr[id.totalnotseeded] <- 0
head(reefdata2$final_abund_pr)

#reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)
for(i in 1:(clusters(g_orig, mode = "weak")$no)){
  reseededclusters$num_reseeded_pr[i] <- sum(reefdata2$final_abund_pr[reefdata2$cl.m == i])
  reseededclusters$per_reseeded_pr[i] <- (reseededclusters$num_reseeded_pr[i]/reseededclusters$clustersize[i])*100
}
png(paste0("Avg",avgdist,"km_CU_Predicted_HistPerReseeded_",numsteps,"steps.png"))
hist(reseededclusters$per_reseeded_pr, main = "Percentages of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

print("final connectivity metric calculating")
#NEW 7.26.2020
#final connectivity metrics
connmat_PCS_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_PCS <- graph.adjacency(as.matrix(connmat_PCS_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]

#number of clusters
cl_PCS <- clusters(g_PCS,mode="weak") #number of clusters, size of each cluster, membership of each cluster
numclust_pr[j] <- cl_PCS$no #99
#geom_mean of cluster size
gmmean_pr[j] <- gm_mean(cl_PCS$csize) #8.917658
#largest cluster
largest_clust_pr[j] <- max(cl_PCS$csize) #1852
#average source strength
so_st_PCS <- ego_size(graph = g_PCS, order = 40, nodes = V(g_PCS), mode = "out")
avgss_pr[j] <- mean(so_st_PCS) #114.0211
#average node degree
dg_PCS <- degree(g_PCS) #getting node degree
avgnd_pr[j] <- mean(dg_PCS) #10.22306

abundance_mat <- NA
plateaupoint <- 500

print("starting catastrophe scenario")
#RESEEDING FROM CATASTROPHE SCENARIO
setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Refuge_Loss/Reseeding"))
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


#Plot the number that were re-seeded over time
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

png(paste0("Avg",avgdist,"km_CU_Catastrophe_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseeded_cr,main="Catastrophe Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs")
dev.off()

#starterreefs + first set seeded - colour the starterreefs green, the rest colour red
#for(i in 2:(numsteps+1)){
#id.abund <- which(abundance_mat[,i] > 0)

#png(paste0("Realized_CU_Catastrophe_Reseed_TS",i-1,"_GR.png"),width=10,height=10,units="cm",res=1000,pointsize=4)
#plot(world, main = paste("Catastrophe Scenario Final Reefs +",i-1,"Time Step"))
#points(reefdata$Longitude,reefdata$Latitude,col=alpha("grey",0.2),pch=20,cex=0.2)
#points(reefdata$Longitude[id.abund],reefdata$Latitude[id.abund],col="red",pch=20,cex=0.2) 
#points(reefdata$Longitude[starterreefs],reefdata$Latitude[starterreefs],col="green",pch=20,cex=0.2)
#dev.off()
#}
#dev.off()

id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)

heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue


#pacific centred versions 
png(paste0("PCentre_Avg",avgdist,"km_CU_Catastrophe_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Catastrophe Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
  points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=20,cex=0.2)
  if(i <= plateaupoint){
    points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=20,cex=0.2)	
  }
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

png(paste0("PCentre_Avg",avgdist,"km_CU_Catastrophe_Reseed_BRB_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Catastrophe Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="light blue",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
  points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)
  if(i <= plateaupoint){
    points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)	
  }
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

png(paste0("PCentre_Avg",avgdist,"km_Catastrophe_allblack.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Catastrophe Scenario"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()


png(paste0("PCentre_Avg",avgdist,"km_CU_Catastrophe_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Catastrophe Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=15,cex=0.03) #cex=0.05
for(i in 2:(numsteps+1)){
  points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=15,cex=0.03)
  if(i <= plateaupoint){
    points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=15,cex=0.03)	
  }
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=15,cex=0.03) #not great but maybe better?
dev.off()


#Histogram of percentage of each of the initial clusters that have been rebuilt
#load in original cluster designations
#reefdata2 <- reefdata[,c(2,3,4)]
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
#range(clusters(g_orig, mode = "weak")$membership)[2] #604 clusters

#reefdata is in the same order as abundance_mat (see: use + construction of id.reefs)
#so can add in another column to reefdata2 with abundance y/n re: re-seeding
#reefdata2 <- reefdata[,c(2,3,4)]
reefdata2$final_abund_cr <- 1
reefdata2$final_abund_cr[id.totalnotseeded] <- 0
head(reefdata2$final_abund_cr)

#reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)
for(i in 1:(clusters(g_orig, mode = "weak")$no)){
  reseededclusters$num_reseeded_cr[i] <- sum(reefdata2$final_abund_cr[reefdata2$cl.m == i])
  reseededclusters$per_reseeded_cr[i] <- (reseededclusters$num_reseeded_cr[i]/reseededclusters$clustersize[i])*100
}
png(paste0("Avg",avgdist,"km_CU_Catastrophe_HistPerReseeded_",numsteps,"steps.png"))
hist(reseededclusters$per_reseeded_cr, main = "Percentages of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

print("starting final connectivity metrics")
#NEW 7.26.2020
#final connectivity metrics
connmat_Catastrophe_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_Catastrophe <- graph.adjacency(as.matrix(connmat_Catastrophe_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]


#number of clusters
cl_Catastrophe <- clusters(g_Catastrophe,mode="weak") #number of clusters, size of each cluster, membership of each cluster
numclust_cr[j] <- cl_Catastrophe$no #175
#geom_mean of cluster size
gmmean_cr[j] <- gm_mean(cl_Catastrophe$csize) #5.803611
#largest cluster
largest_clust_cr[j] <- max(cl_Catastrophe$csize) #1038
#average source strength
so_st_Catastrophe <- ego_size(graph = g_Catastrophe, order = 40, nodes = V(g_Catastrophe), mode = "out")
avgss_cr[j] <- mean(so_st_Catastrophe) #115.4481
#average node degree
dg_Catastrophe <- degree(g_Catastrophe) #getting node degree
avgnd_cr[j] <- mean(dg_Catastrophe) #9.710095

abundance_mat <- NA

initialreefs = 12292 - final_reef
png(paste0("Avg",avgdist,"km_CU_ReseedFull_",numsteps,"steps_allscenarios_totalpercentages_norandom_no50r.png"))
plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100),main="Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "#1d457f", type = 'l', ylim = c(20,100), lwd = 3) #blue
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_cr + initialreefs)/12292)*100), col  = "#f9ad2a", lwd = 3) #darkgoldenrod1
legend("bottomright", c("PCS Scenario", "Refuge-Loss Scenario"), col = c("#1d457f","#f9ad2a"), lty = c(1,1), lwd = c(3,3))
dev.off()

}

save(plateaupoint_pr_list, numreseeded_pr_list, plateaupoint_cr_list, numreseeded_cr_list, numclust_pr, gmmean_pr, largest_clust_pr, avgss_pr, avgnd_pr, numclust_cr, gmmean_cr, largest_clust_cr, avgss_cr, avgnd_cr, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/AllReseedingresults_norandom.RData")


