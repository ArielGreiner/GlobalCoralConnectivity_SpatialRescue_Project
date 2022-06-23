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


#load connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA


#thresholding by 3.3% probability (see 'Comparing_Options' file, note: connmat_reduced@x turns the matrix into a vector)
connmat_reduced@x[connmat_reduced@x <= 0.033] <- 0
#check - okay it's actually fine
which(connmat_reduced[which(connmat_reduced <= 0.033)] > 0)

#method 2
connmat_thresholder <- connmat_reduced
#test <- matrix(c(1,3,2,4), nrow = 2, ncol = 2)
#test[test <= 2] <- 0
#     [,1] [,2]
#[1,]    0    0
#[2,]    3    4
connmat_thresholder[connmat_thresholder <= 0.033] <- 0
connmat_reduced <- connmat_thresholder
#check
connmat_reduced[which(connmat_reduced <= 0.033)] #same weirdness (no actually it's fine)

#load in the centroids, NOTE: NEED TO FIX THE LONGITUDE IN THESE ONES (>180 ones need to have -360 subtracted)
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/centroids.RData')
#head(centroids) 

#reefdata is what needs to be replaced, load in new .csv file
#note: the scores from Hawthorne were in a different order than those from Marco but once they got spatially joined into 12292 cells, it's fine
reefdata_full <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/RecalculatingResultswithNewScores_2.2020/scaled_scores_latlong_reproj_SJ.csv")
head(reefdata_full)
dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)]
head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")

load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/cells12292.RData') #this includes 'id.noReef' and 'b' (shapefile)
reefdata_old <- b@data
head(reefdata_old)


#remove the rows corresponding to the grid cells that have no reef SA
centroids_reduced <- centroids[-id.noReef,]
dim(centroids_reduced) #12292 
#save(centroids_reduced, file = "centroids_reduced.RData")

#what is the distribution of the scores?
#hist(reefdata$score) #centred on 0, -1 to 0.2

#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360
centroids_reduced$Longitude[centroids_reduced$Longitude > 180] <- centroids_reduced$Longitude[centroids_reduced$Longitude > 180] - 360

reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre
#reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] <- reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] - 360


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

#COPY ALL OF THIS AND MAKE ANOTHER COPY BELOW FOR THE OTHER SCENARIO

#remove from lowest (worst) to highest (best) score, in batches of num_remove (434)
#calculating the metrics as we go
ord <- order(reefdata$score) #lowest to highest (decreasing = FALSE is default)
#ord <- order(reefdata$score, decreasing = TRUE) #refuge-loss
fr <- so_st <- cl <- dg <- df <- list() #initialize blank lists, each iteration adding in a new dataframe to each list
reefdata1 <- reefdata[,c(2,3,4,5,7,8,9,22)] #each time, just going to remove from the complete dataframe
#i goes to 20 - won't get down to 0 reefs, but that's what we want
num_remove <- 434 #remove 434 the first time, 434 the next run (total of 868 removed), etc
for(i in 1 : num_iter) {
    final_reef <- (i*num_remove)  
    connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
    df[[i]] <- reefdata1[-ord[1:final_reef],]
    g <- graph.adjacency(as.matrix(connmat_i), weighted = TRUE)
    cl[[i]] <- clusters(g,mode="weak") #number of clusters, size of each cluster, membership of each cluster
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
png(file="Source_strength_b_last.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avg_sost, xlab = "% of reefs remaining", ylab = "Avg Source Strength", type = 'l', xlim = rev(range(nrf_remaining)),main = "3.3% screened", ylim =  c(0,100))
dev.off()

#Global Distribution of Source Strength across removal 'time'
maxval <- range(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))[2]
cols_source <- tim.colors(maxval+1)
for(i in 1:num_iter){
df[[i]]$sost_col <- cols_source[so_st[[i]]]
png(file=paste0("GlobalSourceStrength_",round(nrf_remaining[i+1],2),"PercentReefsRemaining_b_last.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot(world, main = paste("Global Distribution of Source Strength",round(nrf_remaining[i+1],2),"% Reefs Remaining"),lwd = 0.000002, col = "light grey")
points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$sost_col,pch=20,cex=0.2)
legend("topleft", c("Out-degree = 1", "Out-degree = 350", "Out-degree = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(20,20,20))
dev.off()

png(file=paste0("GlobalSourceStrength_",round(nrf_remaining[i+1],2),"PercentReefsRemaining_b_last_small.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot(world, main = paste("Global Distribution of Source Strength",round(nrf_remaining[i+1],2),"% Reefs Remaining"),lwd = 0.000002, col = "light grey")
points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$sost_col,pch=15,cex=0.03)
legend("topleft", c("Out-degree = 1", "Out-degree = 350", "Out-degree = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(15,15,15))
dev.off()
}
dev.off()

#checking source strength distribution at the last step
range(so_st[[20]]) #1 to 190
#hist(so_st[[20]])
#for(i in 1:20){
#print(median(so_st[[i]])) 
#}
#quantile(so_st[[20]], c(0.25,0.5,0.75)) #6  16  42 
#boxplot(so_st[[20]])

#Plot how number of clusters changes across removal 'time'
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
png(file="Num_clusters_b_last.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = nclo, xlab = "% of reefs remaining", ylab = "Number of clusters", type = 'l', main = "3.3% screened", xlim = rev(range(nrf_remaining)), ylim = c(150,1450)) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()
#png(file="Num_clusters_b_last.png",width=12,height=10,units="cm",res=300)
#plot(x = nrf_remaining,y = nclo, xlab = "Number of reefs remaining", ylab = "#Clusters", main = "Number of Clusters", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
#dev.off()

#NEW on 7.26.2020 - geometric_mean version
#NEW on 7.15.2020 - mean version
#Plot how the average size of the clusters changes across removal 'time'
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
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
png(file="GMeanSize_clusters_b_first.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avgsize_cluster, xlab = "% of reefs remaining", ylab = "Geom_Mean Cluster Size", main = "Geom_Mean Cluster Size, 3.3% screened", xlim = rev(range(nrf_remaining)), type = 'l', ylim = c(0,5)) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
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
png(file="Node_degree_b_last.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", type = 'l', xlim = rev(range(nrf_remaining)),main = "3.3% screened", ylim = c(3,10))
dev.off()


# Map the clusters
world <- readOGR(dsn = "~/Dropbox/University of Toronto/Research Related/Code from Marco 6.2018 mediterranean larval connectivity/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")
nrf_remov <-  rep(NA,(num_iter+1))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
    nrf_remov[i] <- ((i-1)*num_remove)
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#how many clusters are there, max
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
clusters(g_orig, mode = "weak")$no 
for(i in 1:num_iter){
print(cl[[i]]$no) #626 -> ...
}

set.seed(3) #added 4.22.2021
cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #the original has the most clusters, useful to keep the colours of the clusters consistent (to the degree that this does that)
#colouring some of them specific colours
cols[166] <- "purple"
cols[30] <- "pink"
cols[1] <- "blue"
cols[53] <- "green"
cols[54] <- "red"
cols[19] <- "black"

#map the clusters when no reefs removed
reefdata1$cl.m <- clusters(g_orig, mode = "weak")$membership
#cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #defined above instead
#reefdata1$cl.c <- cols[reefdata1$cl.m * 10]
reefdata1$cl.c <- cols[reefdata1$cl.m]
png(file=paste0(nrf_remaining[1],"PercentReefs_Remaining_3.3screened.png"),width=10,height=10,units="cm",res=1000,pointsize=4)
plot(world, main = paste(nrf_remaining[1],"%Reefs Remaining, 3.3%"))
#points(reefdata1$Longitude,reefdata1$Latitude,col=reefdata1$cl.c,pch=20,cex=0.2)
points(reefdata1$Longitude,reefdata1$Latitude,col=reefdata1$cl.m,pch=20,cex=0.2)
dev.off()

png(file="AllReefs_Remaining_3.3screened.png",width=10,height=10,units="cm",res=1000, pointsize=4)
plot(world, main = paste(nrf_remaining[1],"PercentReefs Remaining, 3.3%"), lwd = 0.000002, col = "light grey")
#points(reefdata1$Longitude,reefdata1$Latitude,col=reefdata1$cl.c,pch=20,cex=0.2)
points(reefdata1$Longitude,reefdata1$Latitude,col=reefdata1$cl.c,pch=20,cex=0.2)
dev.off()

png(file="All_Reefs_Remaining_3.3screened_small.png",width=20,height=20,units="cm",res=1500)
plot(world, main = paste(nrf_remaining[1],"PercentReefs Remaining, 3.3%"), lwd = 0.000002, col = "light grey")
points(reefdata1$Longitude,reefdata1$Latitude,col=reefdata1$cl.c,pch=15,cex=0.03)
dev.off()

plot(world, main = "6thbiggestcluster_AllReefs Remaining, 3.3%", lwd = 0.000002, col = "light grey")
points(reefdata1$Longitude[reefdata1$cl.m == 1],reefdata1$Latitude[reefdata1$cl.m == 1],col=reefdata1$cl.c[reefdata1$cl.m == 1],pch=20,cex=0.2)


#add a column to the dataframe IDing which grid cells are in which clusters + then assigning colours (now w/n forloop)
for(i in 1:num_iter){
df[[i]]$cl.m <- cl[[i]]$membership
#cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#cols <- sample(tim.colors(cl[[i]]$no)) #defined above instead
#df[[i]]$cl.c <- cols[df[[i]]$cl.m * 10]
df[[i]]$cl.c <- cols[df[[i]]$cl.m]
png(file=paste0(nrf_remaining[i+1],"%Reefs_Remaining_3.3screened.png"),width=10,height=10,units="cm",res=1000,pointsize=4)
plot(world, main = paste(nrf_remaining[i+1],"Reefs Remaining"))
points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
dev.off()
}
dev.off()

#PCentred maps version of the above
for(i in 1:num_iter){
df[[i]]$cl.m <- cl[[i]]$membership
#cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#cols <- sample(tim.colors(cl[[i]]$no)) #defined above instead
#df[[i]]$cl.c <- cols[df[[i]]$cl.m * 10]
df[[i]]$cl.c <- cols[df[[i]]$cl.m]
png(file=paste0(nrf_remaining[i+1],"PercentReefs_Remaining_3.3screened.png"),width=10,height=10,units="cm",res=1000) #NEED TO RERUN bc got rid of pointsize since last saved
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
   fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
dev.off()
}
dev.off()

#final maps - removed reefs coloured grey
library(scales)
i=20
df[[i]]$cl.m <- cl[[i]]$membership
df[[i]]$cl.c <- cols[df[[i]]$cl.m]
png(file=paste0("Final_Predicted_PercentReefs_Remaining_3.3screened.png"),width=10,height=10,units="cm",res=1000) 
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected, reefdata$Latitude, col = "#eef3ff", pch = 20, cex = 0.1) #alpha("light blue",0.2) #0.05 didn't work for the random map
points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
dev.off()


