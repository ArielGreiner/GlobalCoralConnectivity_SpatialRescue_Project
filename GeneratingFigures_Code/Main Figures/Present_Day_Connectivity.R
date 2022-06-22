library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(maps)



#load connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA


#thresholding by 3.3% probability (see 'Comparing_Options' file, note: connmat_reduced@x turns the matrix into a vector)
connmat_reduced@x[connmat_reduced@x <= 0.033] <- 0

#load in the centroids, NOTE: NEED TO FIX THE LONGITUDE IN THESE ONES (>180 ones need to have -360 subtracted)
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/centroids.RData')
#head(centroids) 
reefdata_full <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/RecalculatingResultswithNewScores_2.2020/scaled_scores_latlong_reproj_SJ.csv")
head(reefdata_full)
dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)]
head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")

#calculate coefficient of variation of scores within each reef cell (sd_score/avg_score in each cell)
#calc_CV <- reefdata_full$SD_rescore/reefdata_full$Avg_rescor
#median(calc_CV[!is.nan(calc_CV)]) # = 0
#calculate median avg score within a reef cell
median(reefdata_full$Avg_rescor)

#remove the rows corresponding to the grid cells that have no reef SA
centroids_reduced <- centroids[-id.noReef,]
dim(centroids_reduced) #12292 
#save(centroids_reduced, file = "centroids_reduced.RData")

#what is the distribution of the scores?
#hist(reefdata$score) #centred on 0, -1 to 0.2

#map the 50 reefs in their clusters #4.12.2021 UPDATE: remove this section??
#reefdata1$cl.m <- clusters(g_orig, mode = "weak")$membership #removed 4.7.2022
#cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #defined above instead
#reefdata1$cl.c <- cols[reefdata1$cl.m * 10]
#reefdata1$cl.c <- cols[reefdata1$cl.m] #removed 4.7.2022
#plot(world, main = paste(nrf_remaining[1],"%Reefs Remaining"))
#points(reefdata1$Longitude,reefdata1$Latitude,col=reefdata1$cl.m,pch=20,cex=0.2)

world <- readOGR(dsn = "~/Dropbox/University of Toronto/Research Related/Code from Marco 6.2018 mediterranean larval connectivity/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

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
#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360
centroids_reduced$Longitude[centroids_reduced$Longitude > 180] <- centroids_reduced$Longitude[centroids_reduced$Longitude > 180] - 360

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre
#reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] <- reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] - 360

#warning, running the setwd line below is sometimes annoying later
setwd("~/Dropbox/University of Toronto/Research Related/WorldwideConnectivity_FinalPlots/")
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#calculating how many clusters are 1 reef cell large
#length(which(clusters(g_orig, mode = "weak")$csize == 1))
#largest 6 cluster sizes: 5494 + 609 + 508 + 461 + 379+ 294 => 7745
#+274+247+123+186+164+200+109+274 => 9322
set.seed(3) #added 4.13.2021, 3 seems better than 2
#cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #too much blue #the original has the most clusters, useful to keep the colours of the clusters consistent (to the degree that this does that)
cols <- sample(rainbow(clusters(g_orig, mode = "weak")$no))
#colouring some of them specific colours
cols[166] <- "purple"
cols[30] <- "pink"
cols[1] <- "blue"
cols[53] <- "green"
cols[54] <- "red"
cols[19] <- "black"
#'newcolours' (4.13.2021)
cols[540] <- "brown"
cols[428] <- "Aquamarine"
cols[40] <- "SpringGreen"
cols[454] <- "orange"
cols[453] <- "violetred2"
cols[326] <- "salmon"
cols[60] <- "darkgoldenrod1"
#from looking at the map under set.seed(3)
cols[156] <- "chocolate4"
cols[170] <- "darkgreen"
cols[138] <- "mediumpurple2"
cols[162] <- "indianred2"
#maybe 438

#map the clusters when no reefs removed
reefdata$cl.m <- clusters(g_orig, mode = "weak")$membership
#cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #defined above instead
#reefdata1$cl.c <- cols[reefdata1$cl.m * 10]
reefdata$cl.c <- cols[reefdata$cl.m]

#order(clusters(g_orig, mode = "weak")$csize) smallest clusters listed first

#looking one cluster at a time
plot.map("world", center=newcentre, col="light grey",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[reefdata$cl.m == 30],reefdata$Latitude[reefdata$cl.m == 30],pch=20,cex=0.2)

plot.map("world", center=newcentre, col="light grey",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.m,pch=20,cex=0.2)
 
png("PCentre_3.3_Allreefs_small_othercolours.png",width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=15,cex=0.03) #cex=0.05
dev.off()

png(paste0("PCentre_3.3_Allreefs_othercolours.png"),width=20,height=20,units="cm",res=1500, pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=20,cex=0.2) #cex=0.05
dev.off()

#plot the 50 reefs in their clusters - for appraisal document, need to figure out how made the other parts of fig3...
id.50reefs <- which(!is.na(reefdata$BCU_Name))
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[id.50reefs],reefdata$Latitude[id.50reefs],col=reefdata$cl.c,pch=20,cex=0.2) #cex=0.05

#Out-neighbourhood sizes map
#7.26.2020 - recovered from 'OutDegree_WhatOrder'
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
cols_source <- tim.colors(710)
reefdata$outdg_col <- cols_source[ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out")]
png(paste0("PCentre_GlobalSourceStrengthMap_timcolors.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$outdg_col,pch=20,cex=0.2)
#legend("topleft", c("Neighbourhood Size = 1", "Neighbourhood Size = 350", "Neighbourhood Size = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(20,20,20))
dev.off()

#finding largest source strength values
#order(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"), decreasing = "TRUE")
#9046  9266  9230  9267  9511  9047  9231  9475  9198  9442  9197  9441  9003  9005  9048  9156  9159  9199
#ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out")[2422] #539,514, 376, 311
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected[2422],reefdata$Latitude[2422],col=reefdata$outdg_col[2422],pch=20,cex=0.2)


png(paste0("PCentre_GlobalSourceStrengthMap_timcolors_smallsquares.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$outdg_col,pch=15,cex=0.03)
#legend("topleft", c("Neighbourhood Size = 1", "Neighbourhood Size = 350", "Neighbourhood Size = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(20,20,20))
dev.off()

##4.7.2022: plotting size of network and node degree histogram
#size of network
hist(log(networksizes,10)) #this makes something similar to what is in the manuscript already
networksizes <- clusters(g_orig, mode = "weak")$csize
png(file="LogNetworkSizeDist_4.7.2022.png",width=10,height=10,units="cm",res=1000,pointsize=10)
hist(log(networksizes), breaks = 70, xlim = c(0,10),ylab = "Frequency", xlab = "Log(Network Size)", main = "")
dev.off()
#old
#png(file="InitiallogClusterSizeDist_3.3.png",width=10,height=10,units="cm",res=1000,pointsize=4)
#hist(log(clusters(g_orig, mode = "weak")$csize, base = 10), xlab = "log(Cluster Size)", breaks = 70, ylab = "Frequency", xlim = c(0,4), main = "Size of Clusters before Removal, 3.3%")
#dev.off()

#node degree
png(file="NodeDegreeDist_4.7.2022.png",width=10,height=10,units="cm",res=1000,pointsize=10)
hist(degree(g_orig), xlab = "Node Degree", ylab = "Frequency", xlim = c(0,45), breaks = 30, main = " ")
dev.off()

#old
png(file="InitialNodeDegreeDist_3.3.png",width=10,height=10,units="cm",res=1000,pointsize=4)
hist(degree(g_orig), xlab = "Node Degree", ylab = "Frequency", xlim = c(0,45), breaks = 30, main = "Node Degree Distribution before Removal, 3.3%")
dev.off()

 
#median and interquartile range of out-neighbourhood size
ONS_presentday <- ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out")
median(ONS_presentday) #47
quantile(ONS_presentday, c(0.25,0.75)) #12, 124

#median and interquartile range of node degree
degree_presentday <- degree(g_orig)
median(degree_presentday) #9
quantile(degree_presentday, c(0.25,0.75)) #6, 12

#median + IQR of standard deviation of scores per reef cell
median(reefdata_full$SD_rescore) #0.002544122
quantile(reefdata_full$SD_rescore, c(0.25,0.75)) #0.0001080441, 0.0051180226 


#FOR CSEE2019 Presentation
png("PCentre_3.3_Allreefsred.png",width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,bg = "dark red",col="black",pch=21,cex=0.3,lwd = 0.1) #cex = 0.08,0.15, lwd=0.02,0.05
dev.off()

png("PCentre_3.3_Allreefsgreen.png",width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,bg = "forestgreen",col="black",pch=21,cex=0.3,lwd = 0.1) #cex = 0.08,0.15, lwd=0.02,0.05
dev.off()

#pdf("PCentre_3.3_Allreefsblack.pdf") #res=1000
png("PCentre_3.3_Allreefsblack.png",width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,bg = "black",col="black",pch=21,cex=0.3,lwd = 0.1) #cex = 0.08,0.15, lwd=0.02,0.05
dev.off()

#colouring by score
scorecols <- tim.colors(range(floor(rank(reefdata$score)))[2])
reefdata$scorecols <- scorecols[floor(rank(reefdata$score))]
png("PCentre_3.3_Allreefs_scorecolour.png",width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$scorecols ,pch=20,cex=0.1) #cex = 0.08,0.15, lwd=0.02,0.05
dev.off()

