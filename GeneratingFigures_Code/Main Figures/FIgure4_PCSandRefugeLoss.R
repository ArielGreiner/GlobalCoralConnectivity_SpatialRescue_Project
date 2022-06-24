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

#thresholding by 3.3% probability, note: connmat_reduced@x turns the matrix into a vector
connmat_reduced@x[connmat_reduced@x <= 0.033] <- 0

#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre

world <- readOGR(dsn = "GeneratingFigures_Code/AdditionalFiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

reefdata2 <- reefdata[,c(2,3,4)]
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)


#PCS SCENARIO RESEEDING 
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
plateaupoint_pr <- plateaupoint

id.abund <- list()
for(i in 2:(numsteps+1)){
id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)

png("GeneratingFigures_Code/Main Figures/Figure4a.png",width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="light blue",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)
if(i <= plateaupoint){
	points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)	
}
}
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) 
dev.off()


##final connectivity metrics for Table S2
connmat_PCS_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_PCS <- graph.adjacency(as.matrix(connmat_PCS_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]

#number of networks
cl_PCS <- clusters(g_PCS,mode="weak") #number of networks, size of each network, membership of each network
cl_PCS$no #99
#geom_mean of network size
gm_mean(cl_PCS$csize) #8.917658
#largest network
max(cl_PCS$csize) #1852
#average source strength
so_st_PCS <- ego_size(graph = g_PCS, order = 40, nodes = V(g_PCS), mode = "out")
mean(so_st_PCS) #114.0211
#average node degree
dg_PCS <- degree(g_PCS) #getting node degree
mean(dg_PCS) #10.22306

#just to reset things
abundance_mat <- NA
plateaupoint <- 500


#RESEEDING FROM REFUGE LOSS SCENARIO
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
plateaupoint_cr <- plateaupoint

id.abund <- list()
for(i in 2:(numsteps+1)){
id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)


png("GeneratingFigures_Code/Main Figures/Figure4c.png",width=20,height=20,units="cm",res=1500) #res=1000
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

##final connectivity metrics for Table S2
connmat_RefugeLoss_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_RefugeLoss <- graph.adjacency(as.matrix(connmat_RefugeLoss_final), weighted = TRUE) 

#number of networks
cl_RefugeLoss <- clusters(g_RefugeLoss,mode="weak") #number of networks, size of each network, membership of each network
cl_RefugeLoss$no #175
#geom_mean of network size
gm_mean(cl_RefugeLoss$csize) #5.803611
#largest network
max(cl_RefugeLoss$csize) #1038
#average source strength
so_st_RefugeLoss <- ego_size(graph = g_RefugeLoss, order = 40, nodes = V(g_RefugeLoss), mode = "out")
mean(so_st_RefugeLoss) #115.4481
#average node degree
dg_RefugeLoss <- degree(g_RefugeLoss) #getting node degree
mean(dg_RefugeLoss) #9.710095

abundance_mat <- NA



