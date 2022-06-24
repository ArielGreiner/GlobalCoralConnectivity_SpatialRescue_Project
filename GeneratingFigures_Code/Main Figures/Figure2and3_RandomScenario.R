library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(ggplot2)
library(maps)

###NOTE: The files needed to run this code were generated using the code in the 'ComputeCanadaCode' and are too large to be uploaded to Github, so this code will not run as is and has been uploaded simply to show how the plots were generated from the data generated from the ComputeCanadaCode files.

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


##Calculating where to stop the destruction
#total SA: 12292*324km^2 = 3,982,608km
#12292/434 = ~28
#SA of BCU reefs (50 reefs)
#length(reefdata$BCU_ID[reefdata$BCU_ID > 0]) #3334, 3334*(18*18) = 1,080,216km^2
#3334/434 = ~8
#so if want to stop when the reef SA = 50 reefs SA, should stop at 20 iterations 
num_iter <- 20
random_iter <- 1000
num_remove <- 434

reefdata$num <- NULL


#Amalgamating Datasets

#cl
#i=1,...,200
#load("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_1.RData")
cl_num_full <- list()
for(i in 1:200){
  load("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_i.RData")
  cl_num_full[[i]] <- cl
  
}


g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
nrf_remov <- sd_sost <- m_sost <- s_sost <- rep(NA,(num_iter+1))
m_sost[1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_sost[i] <- 0
  s_sost[i] <- 0
  for(k in 1:200){
  load(paste0("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/so_st/RandomScenario_3.3_so_st_1kreps_",k,".RData"))
  print(paste("i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate the mean
    m_sost[i] <- m_sost[i] + (mean(so_st[[((i-1)+(num_iter*(j-1)))]])/random_iter)
  }
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_sost[i] <- s_sost[i] + ((mean(so_st[[((i-1)+(num_iter*(j-1)))]]) - m_sost[i])^2)
  }
  }
  #then calc sqrt(summation/(n-1))
  sd_sost[i] <- sqrt(s_sost[i]/(random_iter-1))
}
#save(sd_sost,m_sost, file = "mean_sd_sost.RData")

nrf_remaining <- (12292 - nrf_remov)/12292 * 100

sost <- list()
for(k in 1:200){
  load(paste0("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/so_st/RandomScenario_3.3_so_st_1kreps_",k,".RData"))
for(i in 1:random_iter){
  sost[[i]] <- rep(NA,(num_iter+1))
  sost[[i]][1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
  print(paste("i = ",i, "k = ",k))
  for(j in 2 : (num_iter+1)) {
    sost[[i]][j] <-mean(so_st[[((j-1)+(num_iter*(i-1)))]])
  }
}
}
#save(sost, file = "sost_raw.RData")

#this generates Fig2k
png(file="GeneratingFigures_Code/Main Figures/Figure2/SourceStrength_Random.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_sost, pch = 20, xlab = "% of reefs remaining", ylab = "Source strength", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,100))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = sost[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_sost, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

#Plot how number of networks changes across removal 'time'
nrf_remov <- m_nclo <- sd_nclo <- s_nclo <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_nclo[1] <- clusters(g_orig, mode = "weak")$no
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_nclo[i] <- 0
  for(k in 1:200){
    load(paste0("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate the mean
    m_nclo[i] <- m_nclo[i] + (cl[[((i-1)+(num_iter*(j-1)))]]$no)/random_iter
  }
  s_nclo[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_nclo[i] <- s_nclo[i] + (((cl[[((i-1)+(num_iter*(j-1)))]]$no) - m_nclo[i])^2)
  }
  }
  #then calc sqrt(summation/(n-1))
  sd_nclo[i] <- sqrt(s_nclo[i]/(random_iter-1))
}
#save(sd_nclo,m_nclo, file = "mean_sd_nclo.RData")

nrf_remaining <- (12292 - nrf_remov)/12292 * 100

nclo <- list()
for(k in 1:200){
  load(paste0("/Volumes/BackupPlus/worldwideconn_random_morereps/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
  print(paste("i = ",i, "k = ",k))
for(i in 1:random_iter){
  nclo[[i]] <- rep(NA,(num_iter+1))
  nclo[[i]][1] <- clusters(g_orig, mode = "weak")$no
  for(j in 2 : (num_iter+1)) {
    nclo[[i]][j] <- cl[[((j-1)+(num_iter*(i-1)))]]$no
  }
}
}
#save(nclo, file = "nclo_raw.RData")

##this generates figure 2e
png(file="GeneratingFigures_Code/Main Figures/Figure2/NumNetworks_Random.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nclo, pch = 20, xlab = "% of reefs remaining", ylab = "Number of clusters", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(150,1450))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = nclo[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_nclo, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

#geometric mean network size
nrf_remov <- s_avgsize_cluster <- m_avgsize_cluster <- sd_avgsize_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_avgsize_cluster[1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_avgsize_cluster[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_avgsize_cluster[i] <- m_avgsize_cluster[i] + (gm_mean(cl[[((i-1)+(num_iter*(j-1)))]]$csize))/random_iter
  }
  s_avgsize_cluster[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_avgsize_cluster[i] <- s_avgsize_cluster[i] + (((gm_mean(cl[[((i-1)+(num_iter*(j-1)))]]$csize)) - m_avgsize_cluster[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_avgsize_cluster[i] <- sqrt(s_avgsize_cluster[i]/(random_iter-1))
}


nrf_remaining <- (12292 - nrf_remov)/12292 * 100

avgsize_cluster <- list()
for(i in 1:random_iter){
  avgsize_cluster[[i]] <- rep(NA,(num_iter+1))
  avgsize_cluster[[i]][1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
  for(j in 2 : (num_iter+1)) {
    avgsize_cluster[[i]][j] <- (gm_mean(cl[[((j-1)+(num_iter*(i-1)))]]$csize))
  }
}

##this makes figure 2h
png(file="GeneratingFigures_Code/Main Figures/Figure2/GMeanSize_networks_Random.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avgsize_cluster, pch = 20, xlab = "% of reefs remaining", ylab = "Avg Cluster Size", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,5))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = avgsize_cluster[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_avgsize_cluster, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()



#Plot how the average node degree changes across removal 'time'
nrf_remov <- s_avg_degree <- sd_avg_degree <- m_avg_degree <- rep(NA,(num_iter+1))
#g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_avg_degree[1] <- mean(degree(g_orig))
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_avg_degree[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_avg_degree[i] <- m_avg_degree[i] + (mean(dg[[((i-1)+(num_iter*(j-1)))]]))/random_iter
  }
  s_avg_degree[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_avg_degree[i] <- s_avg_degree[i] + (((mean(dg[[((i-1)+(num_iter*(j-1)))]])) - m_avg_degree[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_avg_degree[i] <- sqrt(s_avg_degree[i]/(random_iter-1))
}

nrf_remaining <- (12292 - nrf_remov)/12292 * 100


avg_degree <- list()
for(i in 1:random_iter){
  avg_degree[[i]] <- rep(NA,(num_iter+1))
  avg_degree[[i]][1] <- mean(degree(g_orig))
  for(j in 2 : (num_iter+1)) {
    avg_degree[[i]][j] <- (mean(dg[[((j-1)+(num_iter*(i-1)))]]))
  }
}

##this makes figure 2b
png(file="GeneratingFigures_Code/Main Figures/Figure2/Node_degree_Random.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avg_degree, pch = 20, xlab = "% of reefs remaining", ylab = "Avg Node Degree", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(3,9))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = avg_degree[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_avg_degree, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

###Figure 3b is made analogously to Figure 3a and 3c, just taking the output of one of the random replicates from the compute canada run. Refer to Figure2and3_PCSandRefugeLossScenarios.R for that code.

