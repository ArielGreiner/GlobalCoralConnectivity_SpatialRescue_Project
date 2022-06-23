library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(colorRamps)
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

connmat_thresholder <- connmat_reduced
#test <- matrix(c(1,3,2,4), nrow = 2, ncol = 2)
#test[test <= 2] <- 0
#     [,1] [,2]
#[1,]    0    0
#[2,]    3    4
connmat_thresholder[connmat_thresholder <= 0.033] <- 0
connmat_reduced <- connmat_thresholder

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

#POSSIBLE THING TO DO: Cluster size of rebuilt clusters? like, size of the clusters as they exist at the end...not the amount of each of the former clusters that has been rebuilt but the size of the clusters as they are (i think the percentage rebuilt might be more informative? idk)

#CHANGE TO TIM.COLOURS?

#RESEEDING FROM 50 REEFS
#Start with just the 50 reefs
numsteps <-50
id.50reefs <- which(reefdata$BCU_ID > 0) #this IDs the rows with 50 reefs
length(id.50reefs) #3334
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = numsteps+1) #10 time steps
abundance_mat[,1][id.50reefs] <- 100 #time step 0: assign all the 50 reefs to an abundance of 100

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
#Plot the number that were re-seeded over time
numreseeded_50r <- rep(NA,(numsteps+1))
numreseeded_50r[1] <- 0
for(i in 2:(numsteps + 1)){
#line below could be simplified to just be length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[id.50reefs,i]) but whatever, same for all of the other scenarios
numreseeded_50r[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[id.50reefs,i][abundance_mat[id.50reefs,i]>0])	
if((numreseeded_50r[i] <= numreseeded_50r[i-1]) & stop == "no"){ #added the 'or' condition bc of some of the other cases needing it
	plateaupoint <- i-1
	stop <- "yes"
	}
}
plateaupoint_50r <- plateaupoint



png(paste0("Realized_CU_50R_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseeded_50r,main="Realized 50 Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs")
dev.off()

#starterreefs + first set seeded - colour the starterreefs green, the rest colour red

#LOTS of plots for 50 steps...
#for(i in 2:101){
#id.abund <- which(abundance_mat[,i] > 0)

#png(paste0("Realized_CU_50R_Reseed_TS",i-1,"_GR.png"),width=10,height=10,units="cm",res=1000,pointsize=4)
#plot(world, main = paste("50 Reefs +",i-1,"Time Step"))
#points(reefdata$Longitude,reefdata$Latitude,col=alpha("grey",0.2),pch=20,cex=0.2)
#points(reefdata$Longitude[id.abund],reefdata$Latitude[id.abund],col="red",pch=20,cex=0.2) 
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2)
#dev.off()
#}
#dev.off()

id.abund <- list()
for(i in 2:(numsteps+1)){
id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
id.abund[[i]] <- id.abund[[i]][-id.50reefs]
}
#heatcols <- diverge_hcl(10)
#heatcols <- heat_hcl(12,h=c(0,-100), l = c(75,40), c = c(40,80), power = 1)
#heatcols <- heat_hcl(12,h=c(0,-100),l=c(75,40),c=c(40,80), power=1)
heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)


png(paste0("Realized_CU_50R_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot(world, main = paste("50 Reefs +", numsteps, " Time Steps"), lwd = 0.000002, col = "gainsboro") #lwd = 0.2
points(reefdata$Longitude[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=15,cex=0.03) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=15,cex=0.03)
if(i <= plateaupoint){
	points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=15,cex=0.03)	
}
}
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2) #misleading
points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="black",pch=15,cex=0.03) #not great but maybe better?
dev.off()

png(paste0("Realized_CU_50R_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot(world, main = paste("50 Reefs +", numsteps, " Time Steps"), lwd = 0.000002, col = "gainsboro") #lwd = 0.2
points(reefdata$Longitude[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=20,cex=0.2)
if(i <= plateaupoint){
	points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=20,cex=0.2)	
}
}
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2) #misleading
points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

#pacific centred versions of the above
png(paste0("PCentre_Realized_CU_50R_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("50 Reefs +", numsteps, " Time Steps"),bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=20,cex=0.2)
if(i <= plateaupoint){
	points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=20,cex=0.2)	
}
}
points(reefdata$Longitude_corrected[id.50reefs],reefdata$Latitude[id.50reefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

png(paste0("PCentre_Realized_CU_50R_Reseed_BRB_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("50 Reefs +", numsteps, " Time Steps"),bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="light blue",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)
if(i <= plateaupoint){
	points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col="red",pch=20,cex=0.2)	
}
}
points(reefdata$Longitude_corrected[id.50reefs],reefdata$Latitude[id.50reefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

png(paste0("PCentre_Realized_CU_50R_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("50 Reefs +", numsteps, " Time Steps"),bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=15,cex=0.03) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=15,cex=0.03)
if(i <= plateaupoint){
	points(reefdata$Longitude_corrected[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=15,cex=0.03)	
}
}
points(reefdata$Longitude_corrected[id.50reefs],reefdata$Latitude[id.50reefs],col="black",pch=15,cex=0.03) #not great but maybe better?
dev.off()


##Proportion of the initial clusters that have been rebuilt?

#Histogram of percentage of each of the initial clusters that have been rebuilt
#load in original cluster designations
#reefdata2 <- reefdata[,c(2,3,4)]
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
#range(clusters(g_orig, mode = "weak")$membership)[2] #604 clusters

#reefdata is in the same order as abundance_mat (see: use + construction of id.reefs)
#so can add in another column to reefdata2 with abundance y/n re: re-seeding
#reefdata2 <- reefdata[,c(2,3,4)]
reefdata2$final_abund_50r <- 1
reefdata2$final_abund_50r[id.totalnotseeded] <- 0
head(reefdata2$final_abund_50r)

#reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)
for(i in 1:(clusters(g_orig, mode = "weak")$no)){
	reseededclusters$num_reseeded_50r[i] <- sum(reefdata2$final_abund_50r[reefdata2$cl.m == i])
	reseededclusters$per_reseeded_50r[i] <- (reseededclusters$num_reseeded_50r[i]/reseededclusters$clustersize[i])*100
}
png(paste0("Realized_CU_50R_HistPerReseeded_",numsteps,"steps.png"))
hist(reseededclusters$per_reseeded_50r, main = "Percentages of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

#NEW 7.26.2020
#final connectivity metrics
connmat_50r_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_50r <- graph.adjacency(as.matrix(connmat_50r_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]

#number of clusters
cl_50r <- clusters(g_50r,mode="weak") #number of clusters, size of each cluster, membership of each cluster
cl_50r$no #41
#geom_mean of cluster size
gm_mean(cl_50r$csize) #21.04574
#largest cluster
max(cl_50r$csize) #1247
#average source strength
so_st_50r <- ego_size(graph = g_50r, order = 40, nodes = V(g_50r), mode = "out")
mean(so_st_50r) #137.6775
#average node degree
dg_50r <- degree(g_50r) #getting node degree
mean(dg_50r) #10.66146

abundance_mat <- NA
plateaupoint <- 500


#RESEEDING FROM PREDICTED SCENARIO
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

png(paste0("Realized_CU_Predicted_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseeded_pr,main="Realized Predicted Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs")
dev.off()


#starterreefs + first set seeded - colour the starterreefs green, the rest colour red
#for(i in 2:(numsteps+1)){
#id.abund <- which(abundance_mat[,i] > 0)

#png(paste0("Realized_CU_Predicted_Reseed_TS",i-1,"_GR.png"),width=10,height=10,units="cm",res=1000,pointsize=4)
#plot(world, main = paste("Predicted Scenario Final Reefs +",i-1,"Time Step"))
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
#heatcols <- diverge_hcl(10)
#heatcols <- heat_hcl(12,h=c(0,-100), l = c(75,40), c = c(40,80), power = 1)
#heatcols <- heat_hcl(12,h=c(0,-100),l=c(75,40),c=c(40,80), power=1)
heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue


png(paste0("Realized_CU_Predicted_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot(world, main = paste("Predicted Scenario +", numsteps, "Time Steps"), lwd = 0.000002, col = "gainsboro") #lwd = 0.2
points(reefdata$Longitude[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=15,cex=0.03) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=15,cex=0.03)
if(i <= plateaupoint){
	points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=15,cex=0.03)	
}
}
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2) #misleading
points(reefdata$Longitude[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=15,cex=0.03) #not great but maybe better?
dev.off()

png(paste0("Realized_CU_Predicted_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot(world, main = paste("Predicted Scenario +", numsteps, "Time Steps"), lwd = 0.000002, col = "gainsboro") #lwd = 0.2
points(reefdata$Longitude[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=20,cex=0.2)
if(i <= plateaupoint){
	points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=20,cex=0.2)	
}
}
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2) #misleading
points(reefdata$Longitude[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

#pacific centred versions of the above
png(paste0("PCentre_Realized_CU_Predicted_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
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

png(paste0("PCentre_Realized_CU_Predicted_Reseed_BRB_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
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

png(paste0("PCentre_Realized_Predicted_allblack.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Predicted Scenario +", numsteps, "Time Steps"),bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()



png(paste0("PCentre_Realized_CU_Predicted_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
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
png(paste0("Realized_CU_Predicted_HistPerReseeded_",numsteps,"steps.png"))
hist(reseededclusters$per_reseeded_pr, main = "Percentages of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

#NEW 7.26.2020
#final connectivity metrics
connmat_PCS_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_PCS <- graph.adjacency(as.matrix(connmat_PCS_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]

#number of clusters
cl_PCS <- clusters(g_PCS,mode="weak") #number of clusters, size of each cluster, membership of each cluster
cl_PCS$no #99
#geom_mean of cluster size
gm_mean(cl_PCS$csize) #8.917658
#largest cluster
max(cl_PCS$csize) #1852
#average source strength
so_st_PCS <- ego_size(graph = g_PCS, order = 40, nodes = V(g_PCS), mode = "out")
mean(so_st_PCS) #114.0211
#average node degree
dg_PCS <- degree(g_PCS) #getting node degree
mean(dg_PCS) #10.22306

abundance_mat <- NA
plateaupoint <- 500


#RESEEDING FROM CATASTROPHE SCENARIO
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

png(paste0("Realized_CU_Catastrophe_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseeded_cr,main="Realized Catastrophe Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs")
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

png(paste0("Realized_CU_Catastrophe_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot(world, main = paste("Catastrophe Scenario +",numsteps,"Time Steps"), lwd = 0.000002, col = "gainsboro") #lwd = 0.2
points(reefdata$Longitude[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=15,cex=0.03) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=15,cex=0.03)
if(i <= plateaupoint){
	points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=15,cex=0.03)	
}
}
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2) #misleading
points(reefdata$Longitude[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=15,cex=0.03) #not great but maybe better?
dev.off()

png(paste0("Realized_CU_Catastrophe_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot(world, main = paste("Catastrophe Scenario +",numsteps,"Time Steps"), lwd = 0.000002, col = "gainsboro") #lwd = 0.2
points(reefdata$Longitude[id.totalnotseeded],reefdata$Latitude[id.totalnotseeded],col="gray67",pch=20,cex=0.2) #cex=0.05
for(i in 2:(numsteps+1)){
points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[ultimateplateau+2],pch=20,cex=0.2)
if(i <= plateaupoint){
	points(reefdata$Longitude[id.abund[[i]]],reefdata$Latitude[id.abund[[i]]],col=heatcols[i],pch=20,cex=0.2)	
}
}
#points(reefdata$Longitude[id.50reefs],reefdata$Latitude[id.50reefs],col="green",pch=20,cex=0.2) #misleading
points(reefdata$Longitude[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()

#pacific centred versions of the above
png(paste0("PCentre_Realized_CU_Catastrophe_Reseed_DiffCol_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
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

png(paste0("PCentre_Realized_CU_Catastrophe_Reseed_BRB_",numsteps,"steps_large.png"),width=20,height=20,units="cm",res=1500) #res=1000
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

png(paste0("PCentre_Realized_Catastrophe_allblack.png"),width=20,height=20,units="cm",res=1500,pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",main = paste("Catastrophe Scenario"),bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[starterreefs],reefdata$Latitude[starterreefs],col="black",pch=20,cex=0.2) #not great but maybe better?
dev.off()


png(paste0("PCentre_Realized_CU_Catastrophe_Reseed_DiffCol_",numsteps,"steps.png"),width=20,height=20,units="cm",res=1500) #res=1000
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
png(paste0("Realized_CU_Catastrophe_HistPerReseeded_",numsteps,"steps.png"))
hist(reseededclusters$per_reseeded_cr, main = "Percentages of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

#NEW 7.26.2020
#final connectivity metrics
connmat_Catastrophe_final <- connmat_reduced[-id.totalnotseeded, -id.totalnotseeded]
g_Catastrophe <- graph.adjacency(as.matrix(connmat_Catastrophe_final), weighted = TRUE) 
#connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
#df[[i]] <- reefdata1[-ord[1:final_reef],]

#number of clusters
cl_Catastrophe <- clusters(g_Catastrophe,mode="weak") #number of clusters, size of each cluster, membership of each cluster
cl_Catastrophe$no #175
#geom_mean of cluster size
gm_mean(cl_Catastrophe$csize) #5.803611
#largest cluster
max(cl_Catastrophe$csize) #1038
#average source strength
so_st_Catastrophe <- ego_size(graph = g_Catastrophe, order = 40, nodes = V(g_Catastrophe), mode = "out")
mean(so_st_Catastrophe) #115.4481
#average node degree
dg_Catastrophe <- degree(g_Catastrophe) #getting node degree
mean(dg_Catastrophe) #9.710095

abundance_mat <- NA

#plotting all of the numreseeeded from the 3 different 'scenarios' on one graph
png(paste0("Realized_CU_ReseedFull_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=((numreseeded_50r/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "green", type = 'l', ylim = c(0,30))
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_pr/12292)*100), col = "brown")
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_cr/12292)*100), col  = "pink")
legend("topleft", c("50 Reefs", "Predicted", "Catastrophe"), col = c("green","brown","pink"), lty = c(1,1,1))
dev.off()

#Histograms of Percent of each cluster Reseeded all together
png(paste0("Realized_CU_Full_HistPerReseeded_",numsteps,"steps.png"))
par(mfrow = c(1,3))
hist(reseededclusters$per_reseeded_50r, main = "50 Reefs", xlab = "Percent Re-seeded", col = "green", ylim = c(0,600))
hist(reseededclusters$per_reseeded_pr, main = "Predicted", xlab = "Percent Re-seeded",col = "brown", ylim = c(0,600))
hist(reseededclusters$per_reseeded_cr, main = "Catastrophe", xlab = "Percent Re-seeded", col = "pink", ylim = c(0,600))
dev.off()

#Histograms of Percent of each cluster Reseeded all together
png(paste0("Realized_CU_Full_HistNumReseeded_",numsteps,"steps.png"))
par(mfrow = c(1,3))
hist(reseededclusters$num_reseeded_50r, main = "50 Reefs", xlab = "Number of Cells Re-seeded", col = "green", ylim = c(0,600), xlim = c(0,4000))
hist(reseededclusters$num_reseeded_pr, main = "Predicted", xlab = "Number of Cells Re-seeded",col = "brown", ylim = c(0,600), xlim = c(0,4000))
hist(reseededclusters$num_reseeded_cr, main = "Catastrophe", xlab = "Number of Cells Re-seeded", col = "pink", ylim = c(0,600), xlim = c(0,4000))
dev.off()

#Histograms of Percent of each cluster Reseeded all together, 0s removed
png(paste0("Realized_CU_Full_HistPerReseeded_no0_",numsteps,"steps.png"))
par(mfrow = c(1,3))
hist(reseededclusters$per_reseeded_50r[which(reseededclusters$per_reseeded_50r > 0)], main = "50 Reefs (no 0s)", xlab = "Percent Re-seeded", col = "green", ylim = c(0,200))
hist(reseededclusters$per_reseeded_pr[which(reseededclusters$per_reseeded_pr > 0)], main = "Predicted (no 0s)", xlab = "Percent Re-seeded",col = "brown", ylim = c(0,200))
hist(reseededclusters$per_reseeded_cr[which(reseededclusters$per_reseeded_cr > 0)], main = "Catastrophe (no 0s)", xlab = "Percent Re-seeded", col = "pink", ylim = c(0,200))
dev.off()

png(paste0("Realized_CU_Full_PerReseeded_",numsteps,"steps.png"))
par(mfrow = c(1,3))
plot(x = reseededclusters$clustersize,y = reseededclusters$per_reseeded_50r, main = "50 Reefs", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "green", ylim = c(0,100), pch=20)
plot(y = reseededclusters$per_reseeded_pr, x = reseededclusters$clustersize, main = "Predicted", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "brown", ylim = c(0,100), pch=20)
plot(y = reseededclusters$per_reseeded_cr, x = reseededclusters$clustersize, main = "Catastrophe", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "pink", ylim = c(0,100), pch=20)
dev.off()

png(paste0("Realized_CU_Full_NumReseeded_",numsteps,"steps.png"))
par(mfrow = c(1,3))
plot(x = reseededclusters$clustersize,y = reseededclusters$num_reseeded_50r, main = "50 Reefs", ylab = "Number of Cells Re-seeded", xlab = "Cluster Size", col = "green", pch=20)
plot(y = reseededclusters$num_reseeded_pr, x = reseededclusters$clustersize, main = "Predicted", ylab = "Number of Cells Re-seeded", xlab = "Cluster Size", col = "brown", pch=20)
plot(y = reseededclusters$num_reseeded_cr, x = reseededclusters$clustersize, main = "Catastrophe", ylab = "Number of Cells Re-seeded", xlab = "Cluster Size", col = "pink", pch=20)
dev.off()

png(paste0("Realized_CU_Full_PerReseeded_no0_",numsteps,"steps.png"))
par(mfrow = c(1,3))
plot(x = reseededclusters$clustersize[which(reseededclusters$per_reseeded_50r > 0)],y = reseededclusters$per_reseeded_50r[which(reseededclusters$per_reseeded_50r > 0)], main = "50 Reefs (no 0s)", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "green", ylim = c(0,100), pch=20)
plot(y = reseededclusters$per_reseeded_pr[which(reseededclusters$per_reseeded_pr > 0)], x = reseededclusters$clustersize[which(reseededclusters$per_reseeded_pr > 0)], main = "Predicted (no 0s)", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "brown", ylim = c(0,100), pch=20)
plot(y = reseededclusters$per_reseeded_cr[which(reseededclusters$per_reseeded_cr > 0)], x = reseededclusters$clustersize[which(reseededclusters$per_reseeded_cr > 0)], main = "Catastrophe (no 0s)", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "pink", ylim = c(0,100), pch=20)
dev.off()

png(paste0("Realized_CU_Full_PerReseeded_no0_nobigcluster",numsteps,"steps.png"))
par(mfrow = c(1,3))
plot(x = reseededclusters$clustersize[which(reseededclusters$per_reseeded_50r > 0 & reseededclusters$clustersize < 1000)],y = reseededclusters$per_reseeded_50r[which(reseededclusters$per_reseeded_50r > 0 & reseededclusters$clustersize < 1000)], main = "50 Reefs (no 0s, no big cluster)", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "green", ylim = c(0,100), pch=20)
plot(y = reseededclusters$per_reseeded_pr[which(reseededclusters$per_reseeded_pr > 0 & reseededclusters$clustersize < 1000)], x = reseededclusters$clustersize[which(reseededclusters$per_reseeded_pr > 0 & reseededclusters$clustersize < 1000)], main = "Predicted (no 0s, no big cluster)", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "brown", ylim = c(0,100), pch=20)
plot(y = reseededclusters$per_reseeded_cr[which(reseededclusters$per_reseeded_cr > 0 & reseededclusters$clustersize < 1000)], x = reseededclusters$clustersize[which(reseededclusters$per_reseeded_cr > 0 & reseededclusters$clustersize < 1000)], main = "Catastrophe (no 0s, no big cluster)", ylab = "Percent Re-seeded", xlab = "Cluster Size", col = "pink", ylim = c(0,100), pch=20)
dev.off()

#RESEEDING FROM RANDOM SCENARIO
set.seed(2)
abundmatlist <- starterreefslist <- numreseededrrlist <- id.totalnotseeded <- perreseededrr <- list()
numsteps <- 50
num_remove <- 434
num_iter <- 20
random_iter <- 1000
final_reef <- (num_iter*num_remove)
plateaupoint <- 500
plateaupoint_rr <- rep(NA,random_iter)
reefdata$num <- seq(1,12292,by=1)
for(j in 1:random_iter){
ord <- sample(seq(1,length(reefdata$score),by=1))
starterreefs <- reefdata$num[-ord[1:final_reef]]
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = (numsteps+1)) #50 time steps
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

abundmatlist[[j]] <- abundance_mat
starterreefslist[[j]] <- starterreefs



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
plateaupoint_rr[j] <- plateaupoint
numreseededrrlist[[j]] <- numreseeded_rr

id.totalnotseeded[[j]] <- which(abundance_mat[,(numsteps+1)] == 0)
}

meanplateauptrr <- mean(plateaupoint_rr)

numreseededrr_mean <- rep(NA,(numsteps+1))
numreseededrr_sd <- rep(NA,(numsteps+1))
set <- rep(NA,random_iter)
for(i in 1:(numsteps+1)){
	for(j in 1:random_iter){
	set[j] <- numreseededrrlist[[j]][i]
	}
	numreseededrr_mean[i] <- mean(set)
	numreseededrr_sd[i] <- sd(set)
	set <- rep(NA,random_iter)	
}


png(paste0("Realized_CU_Random_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseededrr_mean,main="Realized Random Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs", type = "n")
for(i in 1:random_iter){
	lines(x = seq(1,(numsteps+1),by=1),y = numreseededrrlist[[i]],col = alpha("blue", 0.2))
}
points(x=seq(1,(numsteps+1),by=1),y=numreseededrr_mean, pch=20)
dev.off()

#NOT SURE IF THIS MAKES SENSE TO PLOT...NEED TO THINK ABOUT
#Histogram of percentage of each of the initial clusters that have been rebuilt
#load in original cluster designations
#reefdata2 <- reefdata[,c(2,3,4)]
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
#range(clusters(g_orig, mode = "weak")$membership)[2] #604 clusters

#reefdata is in the same order as abundance_mat (see: use + construction of id.reefs)
#so can add in another column to reefdata2 with abundance y/n re: re-seeding
#reefdata2 <- reefdata[,c(2,3,4)]

for(j in 1:random_iter){
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

perreseededrr_mean <- rep(NA,(numsteps+1))
perreseededrr_sd <- rep(NA,(numsteps+1))
set <- rep(NA,random_iter)
for(i in 1:(numsteps+1)){
	for(j in 1:random_iter){
	set[j] <- perreseededrr[[j]][i]
	}
	perreseededrr_mean[i] <- mean(set)
	perreseededrr_sd[i] <- sd(set)
	set <- rep(NA,random_iter)	
}

png(paste0("Realized_CU_MeanRandom_HistPerReseeded_",numsteps,"steps.png"))
hist(perreseededrr_mean, main = "Mean Percentage of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
dev.off()

#NEW 7.26.2020
#final connectivity metrics
connmat_Random_final <- g_Random <- cl_Random <- so_st_Random <- dg_Random <- list()
meanclust <- maxclust <- avgsost <- avgdeg <- numclust <- rep(NA,1000)
for(i in 1:random_iter){
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
mean(numclust) #231.937
mean(meanclust) #6.036148
mean(maxclust) #3042.69
mean(avgsost) #105.0341
mean(avgdeg, na.rm = T) #NaN (UGH)

abundance_mat <- NA
plateaupoint <- 500

#plotting all of the numreseeeded from the 4 different 'scenarios' on one graph
png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios.png"))
plot(x=seq(1,(numsteps+1),by=1),y=((numreseeded_50r/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "green", type = 'l', ylim = c(0,60))
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_pr/12292)*100), col = "brown")
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_cr/12292)*100), col  = "pink")
lines(x=seq(1,(numsteps+1),by=1),y=((numreseededrr_mean/12292)*100), col  = "red",lty =3)
legend("bottomright", c("50 Reefs", "Predicted", "Catastrophe","Random"), col = c("green","brown","pink","red"), lty = c(1,1,1,3))
dev.off()

#plotting all of the numreseeeded from the 4 different 'scenarios' on one graph - PERCENTAGE OF RE-SEEDED NOT INCLUDING STARTER REEFS
numpossiblereseeded50r <- 12292-length(id.50reefs)
numpossiblereseeded <- 12292-(12292 - final_reef) #12292 - final_reef = number started with in the predicted and catastrophe situations
png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_fixedpercentages.png"))
plot(x=seq(1,(numsteps+1),by=1),y=((numreseeded_50r/numpossiblereseeded50r)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "green", type = 'l', ylim = c(0,80))
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_pr/numpossiblereseeded)*100), col = "brown")
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_cr/numpossiblereseeded)*100), col  = "pink")
lines(x=seq(1,(numsteps+1),by=1),y=((numreseededrr_mean/numpossiblereseeded)*100), col  = "red",lty =3)
legend("bottomright", c("50 Reefs", "Predicted", "Catastrophe","Random"), col = c("green","brown","pink","red"), lty = c(1,1,1,3))
dev.off()

#for appraisal
#plotting all of the numreseeeded from the 4 different 'scenarios' on one graph - total seeded as a PERCENTTAGE OF TOTAL REEFS
#12292 - final_reef = number started with in the predicted and catastrophe situations
initialreefs = 12292 - final_reef
png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_totalpercentages.png"))
plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_50r+length(id.50reefs))/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "green", type = 'l', ylim = c(20,100), lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100), col = "blue",lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_cr + initialreefs)/12292)*100), col  = "orange", lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseededrr_mean + initialreefs)/12292)*100), col  = "red",lty =3, lwd = 3)
legend("bottomright", c("50 Reefs", "PCS Hyp", "Cautionary Hyp","Random Hyp"), col = c("green","blue","orange","red"), lty = c(1,1,1,3), lwd = c(3,3,3,3))
dev.off()

#NEW 7.27.2020 - adding thin lines around the random scenario one
#plotting all of the numreseeeded from the 4 different 'scenarios' on one graph - total seeded as a PERCENTTAGE OF TOTAL REEFS
#12292 - final_reef = number started with in the predicted and catastrophe situations
initialreefs = 12292 - final_reef
png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_totalpercentages_allrandom.png"))
plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_50r+length(id.50reefs))/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "green", type = 'l', ylim = c(20,100), lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100), col = "blue",lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_cr + initialreefs)/12292)*100), col  = "orange", lwd = 3)
for(i in 1:random_iter){
  lines(x = seq(1,(numsteps+1),by=1),y = (((numreseededrrlist[[i]] + initialreefs)/12292)*100),col = alpha("pink", 0.2),lwd=1)
}
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseededrr_mean + initialreefs)/12292)*100), col  = "red",lty =3, lwd = 3)
legend("bottomright", c("50 Reefs", "PCS Hyp", "Worst Hyp","Mean(Random Hyp)", "Random Hyp"), col = c("green","blue","orange","red","pink"), lty = c(1,1,1,3,1), lwd = c(3,3,3,3,1))
dev.off()

#no 50 reefs for GCC poster
png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_totalpercentages_no50r.png"))
plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100), col = "blue",lwd = 3,main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", type = 'l', ylim = c(20,100))
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_cr + initialreefs)/12292)*100), col  = "orange", lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseededrr_mean + initialreefs)/12292)*100), col  = "red",lty =3, lwd = 3)
legend("bottomright", c("PCS Hyp", "Cautionary Hyp","Random Hyp"), col = c("blue","orange","red"), lty = c(1,1,3), lwd = c(3,3,3))
dev.off()

#not including the 50 reefs scenario

png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_fixedpercentages_no50r.png"))
plot(x=seq(1,(numsteps+1),by=1),y=((numreseeded_pr/numpossiblereseeded)*100), col = "green",main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", type = 'l', ylim = c(0,80))
lines(x=seq(1,(numsteps+1),by=1),y=((numreseeded_cr/numpossiblereseeded)*100), col  = "brown")
lines(x=seq(1,(numsteps+1),by=1),y=((numreseededrr_mean/numpossiblereseeded)*100), col  = "red",lty =3)
legend("bottomright", c("Predicted", "Catastrophe","Random"), col = c("green","brown","red"), lty = c(1,1,3))
dev.off()


