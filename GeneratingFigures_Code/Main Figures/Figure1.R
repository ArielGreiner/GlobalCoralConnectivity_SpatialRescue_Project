library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(maps)

setwd("~/Github/GlobalCoralConnectivity_SpatialRescue_Project/")

#load connectivity matrix
load('OriginalDataConversion/ConvertingWoodetal2014Data/connmat_reduced.RData')
#dim(connmat_reduced) #12292x12292 with >0 reef SA


#thresholding by 3.3% probability, note: connmat_reduced@x turns the matrix into a vector
connmat_reduced@x[connmat_reduced@x <= 0.033] <- 0

#load in scores and locations of reef cells
reefdata_full <- read.csv(file = "OriginalDataConversion/GeneratingScoresfrom50ReefsScores/scaled_scores_latlong_reproj_SJ.csv")
#head(reefdata_full)
#dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)] #removing rows that aren't needed + re-organizing row orders
#head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")

world <- readOGR(dsn = "GeneratingFigures_Code/AdditionalFiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

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

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre
reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre


#turn connectivity matrix into a graph object
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)

###Figure 1a
set.seed(3) 
cols <- sample(rainbow(clusters(g_orig, mode = "weak")$no))
#colouring some of the networks specific colours
cols[166] <- "purple"
cols[30] <- "pink"
cols[1] <- "blue"
cols[53] <- "green"
cols[54] <- "red"
cols[19] <- "black"
cols[540] <- "brown"
cols[428] <- "Aquamarine"
cols[40] <- "SpringGreen"
cols[454] <- "orange"
cols[453] <- "violetred2"
cols[326] <- "salmon"
cols[60] <- "darkgoldenrod1"
cols[156] <- "chocolate4"
cols[170] <- "darkgreen"
cols[138] <- "mediumpurple2"
cols[162] <- "indianred2"

#map the networks
reefdata$cl.m <- clusters(g_orig, mode = "weak")$membership
reefdata$cl.c <- cols[reefdata$cl.m]

png(paste0("GeneratingFigures_Code/Main Figures/Figure1a.png"),width=20,height=20,units="cm",res=1500, pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=20,cex=0.2) #cex=0.05
dev.off()

#EDIT 11.27.2023: Added a higher res version of the figure, as people keep asking for it
png(paste0("GeneratingFigures_Code/Main Figures/Figure1a_highres.png"),width=20,height=20,units="cm",res=3000, pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=20,cex=0.2) #cex=0.05
dev.off()

###Figure 1b: Source Strength map
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
cols_source <- tim.colors(710)
reefdata$outdg_col <- cols_source[ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out")]
png("GeneratingFigures_Code/Main Figures/Figure1b.png",width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$outdg_col,pch=20,cex=0.2)
dev.off()

###Figure 1c and 1d
#size of network
networksizes <- clusters(g_orig, mode = "weak")$csize
png(file="GeneratingFigures_Code/Main Figures/Figure1c.png",width=10,height=10,units="cm",res=1000,pointsize=10)
hist(log(networksizes), breaks = 70, xlim = c(0,10),ylab = "Frequency", xlab = "Log(Network Size)", main = "")
dev.off()

#node degree
png(file="GeneratingFigures_Code/Main Figures/Figure1d.png",width=10,height=10,units="cm",res=1000,pointsize=10)
hist(degree(g_orig), xlab = "Node Degree", ylab = "Frequency", xlim = c(0,45), breaks = 30, main = " ")
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

#calculate median avg score within a reef cell
median(reefdata_full$Avg_rescor) #0.003442668

