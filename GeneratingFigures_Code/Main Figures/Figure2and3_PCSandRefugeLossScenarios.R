library(rgdal)
library(Matrix)
library(igraph)
library(fields)
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

setwd("~/Github/GlobalCoralConnectivity_SpatialRescue_Project/")

#load connectivity matrix
load('OriginalDataConversion/ConvertingWoodetal2014Data/connmat_reduced.RData')

#thresholding by 3.3% probability, note: connmat_reduced@x turns the matrix into a vector
connmat_reduced@x[connmat_reduced@x <= 0.033] <- 0

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


##Calculating where to stop the destruction
#total SA: 12292*324km^2 = 3,982,608km
#12292/434 = ~28
#SA of BCU reefs (50 reefs)
length(reefdata$BCU_ID[reefdata$BCU_ID > 0]) #3334, 3334*(18*18) = 1,080,216km^2
#3334/434 = ~8
#so if want to stop when the reef SA = 50 reefs SA, should stop at 20 iterations 
num_iter <- 20

g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)

gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###PCS Scenario
#remove from lowest (worst) to highest (best) score, in batches of num_remove (434)
#calculating the metrics as we go
ord <- order(reefdata$score) #lowest to highest (decreasing = FALSE is default)
fr <- so_st <- cl <- dg <- df <- list() #initialize blank lists, each iteration adding in a new dataframe to each list
reefdata1 <- reefdata[,c(2,3,4,5,7,8,9,22)] #each time, just going to remove from the complete dataframe
#i goes to 20 - won't get down to 0 reefs, but that's what we want
num_remove <- 434 #remove 434 the first time, 434 the next run (total of 868 removed), etc
for(i in 1 : num_iter) {
    final_reef <- (i*num_remove)  
    connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
    df[[i]] <- reefdata1[-ord[1:final_reef],]
    g <- graph.adjacency(as.matrix(connmat_i), weighted = TRUE)
    cl[[i]] <- clusters(g,mode="weak") #number of networks, size of each network, membership of each network
    dg[[i]] <- degree(g) #getting node degree
    so_st[[i]] <- ego_size(graph = g, order = 40, nodes = V(g), mode = "out")
    fr[[i]] <- length(reefdata1$BCU_ID[reefdata1$BCU_ID > 0][-ord[1:final_reef]])
}

#Plot how average source strength changes across removal 'time'
nrf_remov <- avg_sost <- rep(NA,(num_iter+1))
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avg_sost[1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    avg_sost[i] <- mean(so_st[[i-1]])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/Source_strength_PCS.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avg_sost, xlab = "% of reefs remaining", ylab = "Avg Source Strength", type = 'l', xlim = rev(range(nrf_remaining)),main = "3.3% screened", ylim =  c(0,100))
dev.off()

#checking source strength distribution at the last step
range(so_st[[20]]) #1 to 190


#Plot how number of networks changes across removal 'time'
nrf_remov <- nclo <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
nclo[1] <- clusters(g_orig, mode = "weak")$no
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    nclo[i] <- cl[[i-1]]$no
}
#nrf_remaining <- 12292 - nrf_remov
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/Num_networks_PCS.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = nclo, xlab = "% of reefs remaining", ylab = "Number of networks", type = 'l', main = "3.3% screened", xlim = rev(range(nrf_remaining)), ylim = c(150,1450)) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

#Plot how the geometric mean size of the networks changes across removal 'time'
nrf_remov <- avgsize_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avgsize_cluster[1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    avgsize_cluster[i] <- gm_mean(cl[[i-1]]$csize)
    print(avgsize_cluster[i])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/GMeanSize_networks_PCS.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avgsize_cluster, xlab = "% of reefs remaining", ylab = "Geom_Mean Network Size", main = "Geom_Mean Network Size, 3.3% screened", xlim = rev(range(nrf_remaining)), type = 'l', ylim = c(0,5)) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()


#Plot how the average node degree changes across removal 'time'
nrf_remov <- avg_degree <- rep(NA,(num_iter+1))
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avg_degree[1] <- mean(degree(g_orig))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    avg_degree[i] <- mean(dg[[i-1]])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/Node_degree_PCS.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", type = 'l', xlim = rev(range(nrf_remaining)),main = "3.3% screened", ylim = c(3,10))
dev.off()


###Figure 3a
world <- readOGR(dsn = "GeneratingFigures_Code/AdditionalFiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")
nrf_remov <-  rep(NA,(num_iter+1))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)

set.seed(3) 
cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #the original has the most clusters, useful to keep the colours of the clusters consistent (to the degree that this does that)
#colouring some of them specific colours
cols[166] <- "purple"
cols[30] <- "pink"
cols[1] <- "blue"
cols[53] <- "green"
cols[54] <- "red"
cols[19] <- "black"

#final maps - removed reefs coloured light blue
library(scales)
i=20
df[[i]]$cl.m <- cl[[i]]$membership
df[[i]]$cl.c <- cols[df[[i]]$cl.m]
png(file=paste0("GeneratingFigures_Code/Main Figures/Figure3a.png"),width=10,height=10,units="cm",res=1000) 
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected, reefdata$Latitude, col = "#eef3ff", pch = 20, cex = 0.1) #alpha("light blue",0.2) #0.05 didn't work for the random map
points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
dev.off()

###Refuge-Loss Scenario 
#remove from highest (best) to lowest (worst) score, in batches of num_remove (434)
#calculating the metrics as we go
ord <- order(reefdata$score, decreasing = TRUE) #refuge-loss
fr <- so_st <- cl <- dg <- df <- list() #initialize blank lists, each iteration adding in a new dataframe to each list
reefdata1 <- reefdata[,c(2,3,4,5,7,8,9,22)] #each time, just going to remove from the complete dataframe
#i goes to 20 - won't get down to 0 reefs, but that's what we want
num_remove <- 434 #remove 434 the first time, 434 the next run (total of 868 removed), etc
for(i in 1 : num_iter) {
    final_reef <- (i*num_remove)  
    connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
    df[[i]] <- reefdata1[-ord[1:final_reef],]
    g <- graph.adjacency(as.matrix(connmat_i), weighted = TRUE)
    cl[[i]] <- clusters(g,mode="weak") #number of networks, size of each network, membership of each network
    dg[[i]] <- degree(g) #getting node degree
    so_st[[i]] <- ego_size(graph = g, order = 40, nodes = V(g), mode = "out")
    fr[[i]] <- length(reefdata1$BCU_ID[reefdata1$BCU_ID > 0][-ord[1:final_reef]])
}

#Plot how average source strength changes across removal 'time'
nrf_remov <- avg_sost <- rep(NA,(num_iter+1))
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avg_sost[1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    avg_sost[i] <- mean(so_st[[i-1]])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/Source_strength_RefugeLoss.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avg_sost, xlab = "% of reefs remaining", ylab = "Avg Source Strength", type = 'l', xlim = rev(range(nrf_remaining)),main = "3.3% screened", ylim =  c(0,100))
dev.off()

#checking source strength distribution at the last step
range(so_st[[20]]) #1 to 190


#Plot how number of networks changes across removal 'time'
nrf_remov <- nclo <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
nclo[1] <- clusters(g_orig, mode = "weak")$no
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    nclo[i] <- cl[[i-1]]$no
}
#nrf_remaining <- 12292 - nrf_remov
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/Num_networks_RefugeLoss.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = nclo, xlab = "% of reefs remaining", ylab = "Number of networks", type = 'l', main = "3.3% screened", xlim = rev(range(nrf_remaining)), ylim = c(150,1450)) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

#Plot how the geometric mean size of the networks changes across removal 'time'
nrf_remov <- avgsize_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avgsize_cluster[1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    avgsize_cluster[i] <- gm_mean(cl[[i-1]]$csize)
    print(avgsize_cluster[i])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/GMeanSize_networks_RefugeLoss.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avgsize_cluster, xlab = "% of reefs remaining", ylab = "Geom_Mean Network Size", main = "Geom_Mean Network Size, 3.3% screened", xlim = rev(range(nrf_remaining)), type = 'l', ylim = c(0,5)) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()


#Plot how the average node degree changes across removal 'time'
nrf_remov <- avg_degree <- rep(NA,(num_iter+1))
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avg_degree[1] <- mean(degree(g_orig))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
    avg_degree[i] <- mean(dg[[i-1]])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="GeneratingFigures_Code/Main Figures/Figure2/Node_degree_RefugeLoss.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", type = 'l', xlim = rev(range(nrf_remaining)),main = "3.3% screened", ylim = c(3,10))
dev.off()


###Figure 3c
world <- readOGR(dsn = "GeneratingFigures_Code/AdditionalFiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")
nrf_remov <-  rep(NA,(num_iter+1))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)

set.seed(3) 
cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #the original has the most clusters, useful to keep the colours of the clusters consistent (to the degree that this does that)
#colouring some of them specific colours
cols[166] <- "purple"
cols[30] <- "pink"
cols[1] <- "blue"
cols[53] <- "green"
cols[54] <- "red"
cols[19] <- "black"

#final maps - removed reefs coloured light blue
library(scales)
i=20
df[[i]]$cl.m <- cl[[i]]$membership
df[[i]]$cl.c <- cols[df[[i]]$cl.m]
png(file=paste0("GeneratingFigures_Code/Main Figures/Figure3c.png"),width=10,height=10,units="cm",res=1000) 
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected, reefdata$Latitude, col = "#eef3ff", pch = 20, cex = 0.1) #alpha("light blue",0.2) #0.05 didn't work for the random map
points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
dev.off()