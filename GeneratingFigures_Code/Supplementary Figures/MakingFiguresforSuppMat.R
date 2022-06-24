library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(ggplot2)
library(maps)

setwd("~/Github/GlobalCoralConnectivity_SpatialRescue_Project/")

#load connectivity matrix
load('OriginalDataConversion/ConvertingWoodetal2014Data/connmat_reduced.RData')


#what do i need to threshold it by for the average dispersal distance to be 10km? 100km?

#matrix of euclidean distances as calculated from the original Sally Wood connectivity matrix
load('GeneratingFigures_Code/AdditionalFiles/Euclid.dist.mat.12292.RData')
#histogram of euclidean distances 
hist(euclid.dist.mat@x,breaks=seq(0,8000,40))

# Check if the two matrices are in the same order (they should because one has been built from the other!)
all(connmat_reduced@i == euclid.dist.mat@i)
all(connmat_reduced@p == euclid.dist.mat@p)
#yes :)

#Method 1: Choose a probability threshold for the connectivity matrix
euclid_dist_threshold <- euclid.dist.mat@x
#probabilities negatively correlated with distances (lower probabilities associated with larger distances)
quantile(connmat_reduced@x, c(0.05,0.25,0.5,0.75,0.9))
#          5%          25%          50%          75%          90% 
#0.0001075269 0.0002150538 0.0006451613 0.0032258065 0.0121505376 
connmat_reduced_threshold <- connmat_reduced@x

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

num_clust_presentday <- rep(NA,27)
size_clust_presentday <- list()


#how many present-day networks

#250km present-day connectivity map
avgdist <- 250
probthresh <- 0.1
connmat_reduced <- connmat_reduced_list[[25]]

#load in scores and locations of reef cells
reefdata_full <- read.csv(file = "OriginalDataConversion/GeneratingScoresfrom50ReefsScores/scaled_scores_latlong_reproj_SJ.csv")
#head(reefdata_full)
#dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)] #removing rows that aren't needed + re-organizing row orders
#head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")


#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre


###Make Figure S1
library(rgdal)
library(scales)
library(maps)

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

g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
set.seed(3) 
cols <- sample(rainbow(clusters(g_orig, mode = "weak")$no)) 
#colouring some of them specific colours 
cols[4] <- "brown"
cols[5] <- "cyan"
cols[6] <- "darkgreen"
cols[8] <- "thistle1"
cols[9] <- "turquoise1"
cols[10] <- "yellow1"
cols[11] <- "sienna1"


#map the networks when no reefs removed
reefdata$cl.m <- clusters(g_orig, mode = "weak")$membership
reefdata$cl.c <- cols[reefdata$cl.m]

png("GeneratingFigures_Code/Supplementary Figures/FigureS1a.png",width=20,height=20,units="cm",res=1500, pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=20,cex=0.2) #cex=0.05
dev.off()

g_250 <- graph.adjacency(as.matrix(connmat_reduced_list[[25]]), weighted = TRUE)
size_clust_presentday_250km <- clusters(g_250, mode = "weak")$csize 

#250km
#num_clust_presentday[25] #15
png("GeneratingFigures_Code/Supplementary Figures/FigureS1b.png")
hist(size_clust_presentday_250km, xlab = "Size of Network", main = " ", breaks = 20, cex.lab = 1.5)
dev.off()


###Generating Figure S2
oddonly <- seq(3,25,by=2)
load("GeneratingFigures_Code/Supplementary Figures/ExtraFiles/finalremovallists_norandom.RData")
load("GeneratingFigures_Code/Supplementary Figures/ExtraFiles/FullRandom_multiprob_removal.RData")
load("GeneratingFigures_Code/Supplementary Figures/ExtraFiles/FullRandom_sd_multiprob_removal.RData") #final_sd_sost_random_list, final_sd_nclo_random_list, final_sd_avgsizecluster_random_list, final_sd_avgnd_random_list

#number of networks 
final_nclo_pcs <- rep(NA,26)
final_nclo_rl <- rep(NA,26)
final_nclo_random <- rep(NA,26)
final_nclo_sd_random <- rep(NA,26)
for(k in 1:26){ 
  final_nclo_pcs[k] <- final_nclo_pcs_list[[k]]
  final_nclo_rl[k] <- final_nclo_rl_list[[k]]
  final_nclo_random[k] <- final_nclo_random_list[[k]]
  final_nclo_sd_random[k] <- final_sd_nclo_random_list[[k]]
}

#odd only with 42
final_nclo_pcs_realized <- final_nclo_pcs_list[[27]]
final_nclo_rl_realized <- final_nclo_rl_list[[27]]
final_nclo_random_realized <- final_nclo_random_list[[27]]
final_nclo_sd_random_realized <- final_sd_nclo_random_list[[27]]
png("GeneratingFigures_Code/Supplementary Figures/FigureS2a.png")
plot(x=seq(30,250, by = 20), y = final_nclo_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Number of networks", col = "#1d457f",pch=20, main = "Number of Networks in Relict Reefs", xaxt='n', ylim = c(0,2500), cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_nclo_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20), y = final_nclo_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_nclo_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20), y = final_nclo_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_nclo_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = final_nclo_random[oddonly]-final_nclo_sd_random[oddonly], y1 = final_nclo_random[oddonly]+final_nclo_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_nclo_random_realized-final_nclo_sd_random_realized, y1 = final_nclo_random_realized+final_nclo_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
legend("topright", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)"), col = c("#1d457f","#f9ad2a","#cc5c76"), lty = c(1,1,1), lwd = c(3,3,3))
dev.off()


#average network size 
final_avgsizecluster_pcs <- rep(NA,26)
final_avgsizecluster_rl <- rep(NA,26)
final_avgsizecluster_random <- rep(NA,26)
final_avgsizecluster_sd_random <- rep(NA,26)
for(k in 1:26){ 
  final_avgsizecluster_pcs[k] <- final_avgsizecluster_pcs_list[[k]]
  final_avgsizecluster_rl[k] <- final_avgsizecluster_rl_list[[k]]
  final_avgsizecluster_random[k] <- final_avgsizecluster_random_list[[k]]
  final_avgsizecluster_sd_random[k] <- final_sd_avgsizecluster_random_list[[k]]
}

#odd with 42km
final_avgsizecluster_pcs_realized <- final_avgsizecluster_pcs_list[[27]]
final_avgsizecluster_rl_realized <- final_avgsizecluster_rl_list[[27]]
final_avgsizecluster_random_realized <- final_avgsizecluster_random_list[[27]]
final_avgsizecluster_sd_random_realized <- final_sd_avgsizecluster_random_list[[27]]
png("GeneratingFigures_Code/Supplementary Figures/FigureS2b.png")
plot(x=seq(30,250, by = 20), y = final_avgsizecluster_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Geometric Mean Network Size", col = "#1d457f",pch=20, main = "Geometric Mean Network Size of Relict Reefs", xaxt = 'n', ylim = c(1,10), cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_avgsizecluster_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgsizecluster_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_avgsizecluster_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgsizecluster_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_avgsizecluster_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = final_avgsizecluster_random[oddonly]-final_avgsizecluster_sd_random[oddonly], y1 = final_avgsizecluster_random[oddonly]+final_avgsizecluster_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_avgsizecluster_random_realized-final_avgsizecluster_sd_random_realized, y1 = final_avgsizecluster_random_realized+final_avgsizecluster_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
legend("topleft", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)"), col = c("#1d457f","#f9ad2a","#cc5c76"), lty = c(1,1,1), lwd = c(3,3,3))
dev.off()

#average source strength  
final_avg_sost_pcs <- rep(NA,26)
final_avg_sost_rl <- rep(NA,26)
final_avg_sost_random <- rep(NA,26)
final_avg_sost_sd_random <- rep(NA,26)
for(k in 1:26){ 
  final_avg_sost_pcs[k] <- final_avg_sost_pcs_list[[k]]
  final_avg_sost_rl[k] <- final_avg_sost_rl_list[[k]]
  final_avg_sost_random[k] <- final_avg_sost_random_list[[k]]
  final_avg_sost_sd_random[k] <- final_sd_sost_random_list[[k]]
}

#odd with 42km
final_avg_sost_pcs_realized <- final_avg_sost_pcs_list[[27]]
final_avg_sost_rl_realized <- final_avg_sost_rl_list[[27]]
final_avg_sost_random_realized <- final_avg_sost_random_list[[27]]
final_avg_sost_sd_random_realized <- final_sd_sost_random_list[[27]]
png("GeneratingFigures_Code/Supplementary Figures/FigureS2c.png")
plot(x=seq(30,250,20), y = final_avg_sost_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Average Source Strength", col = "#1d457f",pch=20, main = "Average Source Strength in Remaining Reefs", ylim = c(0,2300), xaxt = 'n', cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_avg_sost_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250,20), y = final_avg_sost_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_avg_sost_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250,20), y = final_avg_sost_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_avg_sost_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250,20), y0 = final_avg_sost_random[oddonly]-final_avg_sost_sd_random[oddonly], y1 = final_avg_sost_random[oddonly]+final_avg_sost_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_avg_sost_random_realized-final_avg_sost_sd_random_realized, y1 = final_avg_sost_random_realized+final_avg_sost_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250,20),labels = T)
legend("topleft", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)"), col = c("#1d457f","#f9ad2a","#cc5c76"), lty = c(1,1,1), lwd = c(3,3,3))
dev.off()

#average node degree 
final_avgnd_pcs <- rep(NA,26)
final_avgnd_rl <- rep(NA,26)
final_avgnd_random <- rep(NA,26)
final_avgnd_sd_random <- rep(NA,26)
for(k in 1:26){ 
  final_avgnd_pcs[k] <- final_avgnd_pcs_list[[k]]
  final_avgnd_rl[k] <- final_avgnd_rl_list[[k]]
  final_avgnd_random[k] <- final_avgnd_random_list[[k]]
  final_avgnd_sd_random[k] <- final_sd_avgnd_random_list[[k]]
}

#odd with 42km
final_avgnd_pcs_realized <- final_avgnd_pcs_list[[27]]
final_avgnd_rl_realized <- final_avgnd_rl_list[[27]]
final_avgnd_random_realized <- final_avgnd_random_list[[27]]
final_avgnd_sd_random_realized <- final_sd_avgnd_random_list[[27]]
png("GeneratingFigures_Code/Supplementary Figures/FigureS2d.png")
plot(x=seq(30,250, by = 20), y = final_avgnd_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Average Node Degree", col = "#1d457f",pch=20, main = "Average Node Degree in Relict Reefs", xaxt = 'n', ylim = c(0,70), cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_avgnd_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgnd_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_avgnd_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgnd_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_avgnd_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = final_avgnd_random[oddonly]-final_avgnd_sd_random[oddonly], y1 = final_avgnd_random[oddonly]+final_avgnd_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_avgnd_random_realized-final_avgnd_sd_random_realized, y1 = final_avgnd_random_realized+final_avgnd_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
legend("topleft", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)"), col = c("#1d457f","#f9ad2a","#cc5c76"), lty = c(1,1,1), lwd = c(3,3,3))
dev.off()

###final plot - putting the four together to make Figure S2
png("GeneratingFigures_Code/Supplementary Figures/FigureS2.png", width = 1000, height = 1000)
par(mfrow = c(2,2))

plot(x=seq(30,250, by = 20), y = final_nclo_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Number of networks", col = "#1d457f",pch=20, main = "(a) Number of Networks in Relict Reefs", xaxt='n', ylim = c(0,2500), cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_nclo_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20), y = final_nclo_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_nclo_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20), y = final_nclo_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_nclo_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = final_nclo_random[oddonly]-final_nclo_sd_random[oddonly], y1 = final_nclo_random[oddonly]+final_nclo_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_nclo_random_realized-final_nclo_sd_random_realized, y1 = final_nclo_random_realized+final_nclo_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)

plot(x=seq(30,250, by = 20), y = final_avgsizecluster_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Geometric Mean Network Size", col = "#1d457f",pch=20, main = "(b) Geometric Mean Network Size of Relict Reefs", xaxt = 'n', ylim = c(1,10), cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_avgsizecluster_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgsizecluster_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_avgsizecluster_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgsizecluster_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_avgsizecluster_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = final_avgsizecluster_random[oddonly]-final_avgsizecluster_sd_random[oddonly], y1 = final_avgsizecluster_random[oddonly]+final_avgsizecluster_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_avgsizecluster_random_realized-final_avgsizecluster_sd_random_realized, y1 = final_avgsizecluster_random_realized+final_avgsizecluster_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)

plot(x=seq(30,250,20), y = final_avg_sost_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Average Source Strength", col = "#1d457f",pch=20, main = "(c) Average Source Strength in Remaining Reefs", ylim = c(0,2300), xaxt = 'n', cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_avg_sost_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250,20), y = final_avg_sost_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_avg_sost_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250,20), y = final_avg_sost_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_avg_sost_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250,20), y0 = final_avg_sost_random[oddonly]-final_avg_sost_sd_random[oddonly], y1 = final_avg_sost_random[oddonly]+final_avg_sost_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_avg_sost_random_realized-final_avg_sost_sd_random_realized, y1 = final_avg_sost_random_realized+final_avg_sost_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")

axis(side = 1, at = seq(30,250,20),labels = T)

plot(x=seq(30,250, by = 20), y = final_avgnd_pcs[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Average Node Degree", col = "#1d457f",pch=20, main = "(d) Average Node Degree in Relict Reefs", xaxt = 'n', ylim = c(0,70), cex.lab = 1.5, cex.main = 1.5)
points(x = 42, y = final_avgnd_pcs_realized, col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgnd_rl[oddonly], col = "#f9ad2a", pch = 20)
points(x = 42, y = final_avgnd_rl_realized, col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20), y = final_avgnd_random[oddonly], col = "#cc5c76", pch = 20)
points(x = 42, y = final_avgnd_random_realized, col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = final_avgnd_random[oddonly]-final_avgnd_sd_random[oddonly], y1 = final_avgnd_random[oddonly]+final_avgnd_sd_random[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = final_avgnd_random_realized-final_avgnd_sd_random_realized, y1 = final_avgnd_random_realized+final_avgnd_sd_random_realized, lwd = 2, lty = 1, col = "#cc5c76")

axis(side = 1, at = seq(30,250, by = 20),labels = T)

dev.off()

###Making Figure S3
#reseeding plots
evenonly <- seq(2,24,by=2)
numsteps <- 50
random_iter <- 100
num_remove <- 434
num_iter <- 20
final_reef <- (num_iter*num_remove)
#plateau point: number re-seeded, generation step
load("GeneratingFigures_Code/Supplementary Figures/ExtraFiles/AllReseedingresults_norandom.RData") #plateaupoint_pr_list etc
numreseeded_pr_final <- rep(NA,26)
numreseeded_cr_final <- rep(NA,26)


#plateaupoint_rr_list <- list()
#numreseeded_rr_list <- list()
load("GeneratingFigures_Code/Supplementary Figures/ExtraFiles/FullRandom_multiprob_reseeding.RData")

#numreseeded
numreseeded_rr_final <- rep(NA,26)
numreseeded_rr_sdfinal <- rep(NA,26)
numreseeded_rr_lowerfinal <- rep(NA,26)
numreseeded_rr_upperfinal <- rep(NA,26)

numreseededrr_mean <- rep(NA,(numsteps+1))
numreseededrr_sd <- rep(NA,(numsteps+1))
numreseededrr_lower <- rep(NA,(numsteps+1))
numreseededrr_upper <- rep(NA,(numsteps+1))

#numreseeded_rr_list
set <- rep(NA,random_iter)

for(k in 1:26){
  numreseededrrlist <- numreseeded_rr_list[[k]]
  for(i in 1:(numsteps+1)){
    for(j in 1:random_iter){
      set[j] <- numreseededrrlist[[j]][i]
    }
    numreseededrr_mean[i] <- mean(set)
    numreseededrr_sd[i] <- sd(set)
    numreseededrr_lower[i] <- range(set)[1]
    numreseededrr_upper[i] <- range(set)[2]
    set <- rep(NA,random_iter)	
  }
  numreseeded_pr_final[k] <- numreseeded_pr_list[[k]][numsteps+1]
  numreseeded_cr_final[k] <- numreseeded_cr_list[[k]][numsteps+1]
  numreseeded_rr_final[k] <- numreseededrr_mean[numsteps+1]
  numreseeded_rr_sdfinal[k] <- numreseededrr_sd[numsteps+1]
  numreseeded_rr_lowerfinal[k] <- numreseededrr_lower[numsteps+1]
  numreseeded_rr_upperfinal[k] <- numreseededrr_upper[numsteps+1]
}

initialreefs = 12292 - final_reef

k = 27
numreseededrrlist <- numreseeded_rr_list[[k]]
for(i in 1:(numsteps+1)){
  for(j in 1:random_iter){
    set[j] <- numreseededrrlist[[j]][i]
  }
  numreseededrr_mean[i] <- mean(set)
  numreseededrr_sd[i] <- sd(set)
  numreseededrr_lower[i] <- range(set)[1]
  numreseededrr_upper[i] <- range(set)[2]
  set <- rep(NA,random_iter)	
}
numreseeded_pr_final_realized <- numreseeded_pr_list[[k]][numsteps+1]
numreseeded_cr_final_realized <- numreseeded_cr_list[[k]][numsteps+1]
numreseeded_rr_final_realized <- numreseededrr_mean[numsteps+1]
numreseeded_rr_sdfinal_realized <- numreseededrr_sd[numsteps+1]
numreseeded_rr_lowerfinal_realized <- numreseededrr_lower[numsteps+1]
numreseeded_rr_upperfinal_realized <- numreseededrr_upper[numsteps+1]

#odd with 42km
png("GeneratingFigures_Code/Supplementary Figures/FigureS3a.png")
plot(x=seq(30,250, by = 20),y=(((numreseeded_pr_final[oddonly] + initialreefs)/12292)*100), xlab = "Mean Dispersal Distance Threshold", ylab = "%  of Reefs Reseeded", main = "Percent of Reefs Re-seeded from Relict Reefs", col = "#1d457f", pch = 20, ylim = c(20,100), xaxt = 'n', cex.lab = 1.5, cex.main = 1.5)
points(x=42,y=(((numreseeded_pr_final_realized + initialreefs)/12292)*100),col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20),y=(((numreseeded_cr_final[oddonly] + initialreefs)/12292)*100),col = "#f9ad2a", pch = 20)
points(x=42,y=(((numreseeded_cr_final_realized + initialreefs)/12292)*100),col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20),y=(((numreseeded_rr_final[oddonly] + initialreefs)/12292)*100),col = "#cc5c76", pch = 20)
points(x=42,y=(((numreseeded_rr_final_realized + initialreefs)/12292)*100),col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = (((numreseeded_rr_lowerfinal[oddonly] + initialreefs)/12292)*100), y1 = (((numreseeded_rr_upperfinal[oddonly] + initialreefs)/12292)*100), lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = (((numreseeded_rr_lowerfinal_realized + initialreefs)/12292)*100), y1 = (((numreseeded_rr_upperfinal_realized + initialreefs)/12292)*100), lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
legend("bottomright", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)"), col = c("#1d457f","#f9ad2a","#cc5c76"), lty = c(1,1,1), lwd = c(3,3,3))
dev.off()



#plateau point
plateaupoint_pr <- rep(NA,26)
plateaupoint_cr <- rep(NA,26)
plateaupoint_rr_mean <- rep(NA,26)
plateaupoint_rr_sd <- rep(NA,26)
plateaupoint_rr_lower <- rep(NA,26)
plateaupoint_rr_upper <- rep(NA,26)

numsteps <- 50
random_iter <- 100

for(k in 1:26){
  plateaupointrrlist <- plateaupoint_rr_list[[k]]
  plateaupoint_pr[k] <- plateaupoint_pr_list[[k]]
  plateaupoint_cr[k] <- plateaupoint_cr_list[[k]]
  plateaupoint_rr_mean[k] <- mean(plateaupointrrlist)
  plateaupoint_rr_sd[k] <- sd(plateaupointrrlist)
  plateaupoint_rr_lower[k] <- range(plateaupointrrlist)[1]
  plateaupoint_rr_upper[k] <- range(plateaupointrrlist)[2]
  set <- rep(NA,random_iter)
}


k = 27
plateaupointrrlist <- plateaupoint_rr_list[[k]]
plateaupoint_pr_realized <- plateaupoint_pr_list[[k]]
plateaupoint_cr_realized <- plateaupoint_cr_list[[k]]
plateaupoint_rr_mean_realized <- mean(plateaupointrrlist)
plateaupoint_rr_sd_realized <- sd(plateaupointrrlist)
plateaupoint_rr_lower_realized <- range(plateaupointrrlist)[1]
plateaupoint_rr_upper_realized <- range(plateaupointrrlist)[2]


#odd with 42km
png("GeneratingFigures_Code/Supplementary Figures/FigureS3b.png")
plot(x=seq(30,250, by = 20),y=plateaupoint_pr[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Plateau Point", main = "# of Generations to Re-seed Reachable Reefs from Relict Reefs", col = "#1d457f", pch = 20, ylim = c(0,45), xaxt = 'n', cex.lab = 1.5, cex.main = 1.5)
points(x=42,y=plateaupoint_pr_realized,col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20),y=plateaupoint_cr[oddonly],col = "#f9ad2a", pch = 20)
points(x=42,y=plateaupoint_cr_realized,col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20),y=plateaupoint_rr_mean[oddonly],col = "#cc5c76", pch = 20)
points(x=42,y=plateaupoint_rr_mean_realized,col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = plateaupoint_rr_lower[oddonly], y1 = plateaupoint_rr_upper[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
segments(x0 = 42, y0 = plateaupoint_rr_lower_realized, y1 = plateaupoint_rr_upper_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
legend("topright", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)"), col = c("#1d457f","#f9ad2a","#cc5c76"), lty = c(1,1,1), lwd = c(3,3,3))
dev.off()

##Full Figure S3
png("GeneratingFigures_Code/Supplementary Figures/FigureS3.png", width = 600, height = 1000)
par(mfrow = c(2,1))

plot(x=seq(30,250, by = 20),y=(((numreseeded_pr_final[oddonly] + initialreefs)/12292)*100), xlab = "Mean Dispersal Distance Threshold", ylab = "%  of Reefs Reseeded", main = "(a) Percent of Reefs Re-seeded from Relict Reefs", col = "#1d457f", pch = 20, ylim = c(20,100), xaxt = 'n', cex.lab = 1.5, cex.main = 1.5)
points(x=42,y=(((numreseeded_pr_final_realized + initialreefs)/12292)*100),col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20),y=(((numreseeded_cr_final[oddonly] + initialreefs)/12292)*100),col = "#f9ad2a", pch = 20)
points(x=42,y=(((numreseeded_cr_final_realized + initialreefs)/12292)*100),col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20),y=(((numreseeded_rr_final[oddonly] + initialreefs)/12292)*100),col = "#cc5c76", pch = 20)
points(x=42,y=(((numreseeded_rr_final_realized + initialreefs)/12292)*100),col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = (((numreseeded_rr_lowerfinal[oddonly] + initialreefs)/12292)*100), y1 = (((numreseeded_rr_upperfinal[oddonly] + initialreefs)/12292)*100), lwd = 2, lty = 1, col = "#cc5c76")
segments(x0 = 42, y0 = (((numreseeded_rr_lowerfinal_realized + initialreefs)/12292)*100), y1 = (((numreseeded_rr_upperfinal_realized + initialreefs)/12292)*100), lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)

plot(x=seq(30,250, by = 20),y=plateaupoint_pr[oddonly], xlab = "Mean Dispersal Distance Threshold", ylab = "Plateau Point", main = "(b) Number of Generations to Re-seed Reachable Reefs from Relict Reefs", col = "#1d457f", pch = 20, ylim = c(0,45), xaxt = 'n', cex.lab = 1.5, cex.main = 1.5)
points(x=42,y=plateaupoint_pr_realized,col = "#1d457f", pch = 10)
points(x=seq(30,250, by = 20),y=plateaupoint_cr[oddonly],col = "#f9ad2a", pch = 20)
points(x=42,y=plateaupoint_cr_realized,col = "#f9ad2a", pch = 10)
points(x=seq(30,250, by = 20),y=plateaupoint_rr_mean[oddonly],col = "#cc5c76", pch = 20)
points(x=42,y=plateaupoint_rr_mean_realized,col = "#cc5c76", pch = 10)
segments(x0 = seq(30,250, by = 20), y0 = plateaupoint_rr_lower[oddonly], y1 = plateaupoint_rr_upper[oddonly], lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
segments(x0 = 42, y0 = plateaupoint_rr_lower_realized, y1 = plateaupoint_rr_upper_realized, lwd = 2, lty = 1, col = "#cc5c76")
axis(side = 1, at = seq(30,250, by = 20),labels = T)
dev.off()
