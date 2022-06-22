#Intro Stuff
#library(rgdal) - hopefully don't need this
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(ggplot2)
print("loaded libraries")

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
print("loaded and edited reefdata")

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

#general things 
num_iter <- 20
random_iter <- 1000
num_remove <- 434
reefdata$num <- NULL


#NEW plotting
print("source strength section...")
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
nrf_remov <- sd_sost <- m_sost <- s_sost <- rep(NA,(num_iter+1))
m_sost[1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
nrf_remov[1] <- 0
print("entering first source strength for-loop")
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_sost[i] <- 0
  s_sost[i] <- 0
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/so_st/RandomScenario_3.3_so_st_1kreps_",k,".RData"))
    print(paste("mean loop, i = ",i, "k = ",k))
    for(j in 1:random_iter){ #calculate the mean
      m_sost[i] <- m_sost[i] + (mean(so_st[[((i-1)+(num_iter*(j-1)))]])/(random_iter*200))
    }
  }
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/so_st/RandomScenario_3.3_so_st_1kreps_",k,".RData"))
    print(paste("sd loop, i = ",i, "k = ",k))
    for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
      s_sost[i] <- s_sost[i] + ((mean(so_st[[((i-1)+(num_iter*(j-1)))]]) - m_sost[i])^2)
    }
  }
  
  #then calc sqrt(summation/(n-1))
  sd_sost[i] <- sqrt(s_sost[i]/((random_iter*200)-1))
}
save(sd_sost,m_sost, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/mean_sd_sost.RData")
print("done first source strength for loop")
print(paste0("sd_sost[1] = ", sd_sost[1],"m_sost[1] = ",m_sost[1]))

nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Source_strength_random_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_sost, pch = 20, xlab = "% of reefs remaining", ylab = "Source Strength", main = "3.3% screened", xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

print("entering 2nd source strength for loop")
sost <- list()
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/so_st/RandomScenario_3.3_so_st_1kreps_",k,".RData"))
  for(i in 1:random_iter){
    sost[[(k-1)*1000 + i]] <- rep(NA,(num_iter+1))
    sost[[(k-1)*1000 + i]][1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
    print(paste("i = ",i, "k = ",k))
    for(j in 2 : (num_iter+1)) {
      sost[[(k-1)*1000 + i]][j] <-mean(so_st[[((j-1)+(num_iter*(i-1)))]])
    }
  }
}
save(sost, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/sost_raw.RData")
print("done 2nd source strength for loop")
print(paste0("sost[[1]][2] = ", sost[[1]][2]))

png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Source_strength_random_lines_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_sost, pch = 20, xlab = "% of reefs remaining", ylab = "Source strength", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,100))
for(i in 1:(random_iter*200)){
  lines(x = nrf_remaining,y = sost[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_sost, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()


plotdata <- data.frame(x= nrf_remaining, y = m_sost, lower = m_sost - sd_sost, upper =  m_sost + sd_sost)

png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Source_strength_random_eb_200k.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  #geom_line(aes(y=plotdata[,5],x=x),color="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Source Strength")+
  ggtitle("3.3% screened")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

dev.off()
print("done plotting source strength stuff")

#TO DO, MAYBE: ^ can make maps of this by taking the average for a particular node at the different time points across replicates? not sure how informative that would be though

#Plot how number of clusters changes across removal 'time'
print("entering cluster-zone")
nrf_remov <- m_nclo <- sd_nclo <- s_nclo <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_nclo[1] <- clusters(g_orig, mode = "weak")$no
nrf_remov[1] <- 0
print("entering 1st clusters for loop")
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_nclo[i] <- 0
  s_nclo[i] <- 0
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("mean loop, i = ",i, "k = ",k))
    for(j in 1:random_iter){ #calculate the mean
      m_nclo[i] <- m_nclo[i] + (cl[[((i-1)+(num_iter*(j-1)))]]$no)/(random_iter*200)
    }
  }
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("sd loop, i = ",i, "k = ",k))
    for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
      s_nclo[i] <- s_nclo[i] + (((cl[[((i-1)+(num_iter*(j-1)))]]$no) - m_nclo[i])^2)
    }
  }
  #then calc sqrt(summation/(n-1))
  sd_nclo[i] <- sqrt(s_nclo[i]/((random_iter*200)-1))
}
save(sd_nclo,m_nclo, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/mean_sd_nclo.RData")
print("finished 1st cluster for loop")
print(paste0("sd_nclo[1] =",sd_nclo[1],"m_nclo[1]=",m_nclo[1]))

nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Num_clusters_random_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nclo, pch = 20, xlab = "% of reefs remaining", ylab = "Number of clusters", main = "3.3% screened", xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()
#png(file="Num_clusters_b_last.png",width=12,height=10,units="cm",res=300)
#plot(x = nrf_remaining,y = nclo, xlab = "Number of reefs remaining", ylab = "#Clusters", main = "Number of Clusters", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
#dev.off()

nclo <- list()
print("entering 2nd cluster for loop")
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
  print(paste("i = ",i, "k = ",k))
  for(i in 1:random_iter){
    nclo[[(k-1)*1000 + i]] <- rep(NA,(num_iter+1))
    nclo[[(k-1)*1000 + i]][1] <- clusters(g_orig, mode = "weak")$no
    for(j in 2 : (num_iter+1)) {
      nclo[[(k-1)*1000 + i]][j] <- cl[[((j-1)+(num_iter*(i-1)))]]$no
    }
  }
}
save(nclo, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/nclo_raw.RData")
print("finished 2nd cluster for loop")
print(paste0("nclo[[1]][2] =",nclo[[1]][2]))

png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Num_clusters_random_lines_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nclo, pch = 20, xlab = "% of reefs remaining", ylab = "Number of clusters", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(150,1450))
for(i in 1:(random_iter*200)){
  lines(x = nrf_remaining,y = nclo[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_nclo, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_nclo, lower = m_nclo - sd_nclo, upper =  m_nclo + sd_nclo)
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Num_clusters_random_eb_200k.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Number of clusters")+
  ggtitle("3.3% screened")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()

#Plot how the size of the largest cluster changes across removal 'time'
nrf_remov <- s_largest_cluster <- m_largest_cluster <- sd_largest_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_largest_cluster[1] <- max(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
print("starting largest cluster for loop")
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_largest_cluster[i] <- 0
  s_largest_cluster[i] <- 0
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("mean loop, i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate the mean
    m_largest_cluster[i] <- m_largest_cluster[i] + (max(cl[[((i-1)+(num_iter*(j-1)))]]$csize))/(random_iter*200)
  }
  }
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("sd loop, i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_largest_cluster[i] <- s_largest_cluster[i] + (((max(cl[[((i-1)+(num_iter*(j-1)))]]$csize)) - m_largest_cluster[i])^2)
  }
  }
  #then calc sqrt(summation/(n-1))
  sd_largest_cluster[i] <- sqrt(s_largest_cluster[i]/((random_iter*200)-1))
}

save(sd_largest_cluster,m_largest_cluster, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/mean_sd_largest_cluster.RData")
print("ending largest cluster for loop")
print(paste0("sd_largest_cluster[1] =",sd_largest_cluster[1],"m_largest_cluster[1]=",m_largest_cluster[1]))

nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Size_clusters_random_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_largest_cluster, xlab = "% of reefs remaining", ylab = "Size of Largest Cluster", main = "Size of Largest Cluster, 3.3% screened", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

largest_cluster <- list()
print("starting largest cluster for loop2")
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
  print(paste("i = ",i, "k = ",k))
for(i in 1:random_iter){
  largest_cluster[[(k-1)*1000 + i]] <- rep(NA,(num_iter+1))
  largest_cluster[[(k-1)*1000 + i]][1] <- max(clusters(g_orig, mode = "weak")$csize)
  for(j in 2 : (num_iter+1)) {
    largest_cluster[[(k-1)*1000 + i]][j] <- (max(cl[[((j-1)+(num_iter*(i-1)))]]$csize))
  }
}
}
save(largest_cluster, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/largest_cluster_raw.RData")
print("ending largest cluster for loop 2")
print(paste0("largest_cluster[[1]][2] =",largest_cluster[[1]][2]))


png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Size_clusters_random_lines_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_largest_cluster, pch = 20, xlab = "% of reefs remaining", ylab = "Size of Largest Cluster", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,6000))
for(i in 1:(random_iter*200)){
  lines(x = nrf_remaining,y = largest_cluster[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_largest_cluster, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_largest_cluster, lower = m_largest_cluster - sd_largest_cluster, upper =  m_largest_cluster + sd_largest_cluster)
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Size_clusters_random_eb_200k.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Cluster Size")+
  ggtitle("3.3% screened")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()

#NEW on 7.26.2020
#change to geometric mean bc Marie-Josee says that that is better
#NEW on 7.15.2020
#Plot how the average size of the clusters changes across removal 'time'
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

nrf_remov <- s_avgsize_cluster <- m_avgsize_cluster <- sd_avgsize_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_avgsize_cluster[1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
print("starting avgsize cluster for loop")
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_avgsize_cluster[i] <- 0
  s_avgsize_cluster[i] <- 0
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("mean loop, i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate the mean
    m_avgsize_cluster[i] <- m_avgsize_cluster[i] + (gm_mean(cl[[((i-1)+(num_iter*(j-1)))]]$csize))/(random_iter*200)
  }
  }
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
    print(paste("sd loop, i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_avgsize_cluster[i] <- s_avgsize_cluster[i] + (((gm_mean(cl[[((i-1)+(num_iter*(j-1)))]]$csize)) - m_avgsize_cluster[i])^2)
  }
  }
  #then calc sqrt(summation/(n-1))
  sd_avgsize_cluster[i] <- sqrt(s_avgsize_cluster[i]/((random_iter*200)-1))
}
save(sd_avgsize_cluster,m_avgsize_cluster, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/mean_sd_avgsize_cluster.RData")
print("ending avg size cluster for loop")
print(paste0("sd_avgsize_cluster[1] =",sd_avgsize_cluster[1],"m_avgsize_cluster[1]=",m_avgsize_cluster[1]))


nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/GMean_clusters_random_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avgsize_cluster, xlab = "% of reefs remaining", ylab = "Avg Cluster Size", main = "Avg Cluster Size, 3.3% screened", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

avgsize_cluster <- list()
print("starting avgsize cluster forloop2")
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/cl/RandomScenario_3.3_cl_1kreps_",k,".RData"))
  print(paste("i = ",i, "k = ",k))
for(i in 1:random_iter){
  avgsize_cluster[[(k-1)*1000 + i]] <- rep(NA,(num_iter+1))
  avgsize_cluster[[(k-1)*1000 + i]][1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
  for(j in 2 : (num_iter+1)) {
    avgsize_cluster[[(k-1)*1000 + i]][j] <- (gm_mean(cl[[((j-1)+(num_iter*(i-1)))]]$csize))
  }
}
}
save(avgsize_cluster, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/largest_cluster_raw.RData")
print("ending avg size cluster for loop 2")
print(paste0("avgsize_cluster[[1]][2] =",avgsize_cluster[[1]][2]))


png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/GMean_clusters_random_lines_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avgsize_cluster, pch = 20, xlab = "% of reefs remaining", ylab = "Avg Cluster Size", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,5))
for(i in 1:(random_iter*200)){
  lines(x = nrf_remaining,y = avgsize_cluster[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_avgsize_cluster, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_avgsize_cluster, lower = m_avgsize_cluster - sd_avgsize_cluster, upper =  m_avgsize_cluster + sd_avgsize_cluster)
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/GMean_clusters_random_eb_200k.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Cluster Size")+
  ggtitle("3.3% screened")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()#Plot how the avg size of the clusters changes across removal 'time'

print("average degree time")
#Plot how the average node degree changes across removal 'time'
nrf_remov <- s_avg_degree <- sd_avg_degree <- m_avg_degree <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_avg_degree[1] <- mean(degree(g_orig))
nrf_remov[1] <- 0
print("starting avg degree for loop 1")
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_avg_degree[i] <- 0
  s_avg_degree[i] <- 0
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/dg/RandomScenario_3.3_dg_1kreps_",k,".RData"))
    print(paste("mean loop, i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate the mean
    m_avg_degree[i] <- m_avg_degree[i] + (mean(dg[[((i-1)+(num_iter*(j-1)))]]))/(random_iter*200)
  }
  }
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/dg/RandomScenario_3.3_dg_1kreps_",k,".RData"))
    print(paste("sd loop, i = ",i, "k = ",k))
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_avg_degree[i] <- s_avg_degree[i] + (((mean(dg[[((i-1)+(num_iter*(j-1)))]])) - m_avg_degree[i])^2)
  }
  }
  #then calc sqrt(summation/(n-1))
  sd_avg_degree[i] <- sqrt(s_avg_degree[i]/((random_iter*200)-1))
}
save(sd_avg_degree,m_avg_degree, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/mean_sd_avg_degree.RData")
print("ending avg degree for loop 1")
print(paste0("sd_avg_degree[1] =",sd_avg_degree[1],"m_avg_degree[1]=",m_avg_degree[1]))


nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Node_degree_random_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,pch = 20, y = m_avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", xlim = rev(range(nrf_remaining)),main = "3.3% screened")
dev.off()

print("starting avg degree for loop 2")
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/dg/RandomScenario_3.3_dg_1kreps_",k,".RData"))
avg_degree <- list()
for(i in 1:random_iter){
  print(paste("i = ",i, "k = ",k))
  avg_degree[[(k-1)*1000 + i]] <- rep(NA,(num_iter+1))
  avg_degree[[(k-1)*1000 + i]][1] <- mean(degree(g_orig))
  for(j in 2 : (num_iter+1)) {
    avg_degree[[(k-1)*1000 + i]][j] <- (mean(dg[[((j-1)+(num_iter*(i-1)))]]))
  }
}
}
save(avg_degree, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/avg_degree_raw.RData")
print("ending avg degree for loop 2")
print(paste0("avg_degree[[1]][2] =",avg_degree[[1]][2]))

png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Node_degree_random_lines_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avg_degree, pch = 20, xlab = "% of reefs remaining", ylab = "Avg Node Degree", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n", ylim = c(3,9))
for(i in 1:(random_iter*200)){
  lines(x = nrf_remaining,y = avg_degree[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_avg_degree, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_avg_degree, lower = m_avg_degree - sd_avg_degree, upper =  m_avg_degree + sd_avg_degree)
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Node_degree_random_eb_200k.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Avg Node Degree")+
  ggtitle("3.3% screened")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()

#Plot how 50 reefs are lost across removal 'time'
print("starting 50 reefs zone")
nrf_remov <- sd_nfr <- m_nfr <- s_nfr <- rep(NA,(num_iter+1))
m_nfr[1] <- length(reefdata$BCU_ID[reefdata$BCU_ID > 0])
nrf_remov[1] <- 0
print("starting 50 reefs for loop 1")
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_nfr[i] <- 0
  s_nfr[i] <- 0
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/df/RandomScenario_3.3_df_1kreps_",k,".RData"))
    print(paste("mean loop, i = ",i, "k = ",k))
    for(j in 1:random_iter){ #calculate the mean
      m_nfr[i] <- m_nfr[i] + (length(df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID[df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID > 0]))/(random_iter*200)
    }
  }
  for(k in 1:200){
    load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/df/RandomScenario_3.3_df_1kreps_",k,".RData"))
    print(paste("sd loop, i = ",i, "k = ",k))
    for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
      s_nfr[i] <- s_nfr[i] + ((length(df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID[df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID > 0]) - m_nfr[i])^2)
    }
  }
  #then calc sqrt(summation/(n-1))
  sd_nfr[i] <- sqrt(s_nfr[i]/((random_iter*200)-1))
}
save(sd_nfr,m_nfr, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/mean_sd_nfr.RData")
print("ending 50reefs for loop 1")
print(paste0("sd_nfr[1] =",sd_nfr[1],"m_nfr[1]=",m_nfr[1]))

nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Num_50r_random_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nfr, pch = 20, xlab = "% of reefs remaining", ylab = "Number of 50 reefs", main = "3.3% screened", xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

nfr <- list()
print("starting 50reefs for loop 2")
for(k in 1:200){
  load(paste0("/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/df/RandomScenario_3.3_df_1kreps_",k,".RData"))
  for(i in 1:random_iter){
    nfr[[(k-1)*1000 + i]] <- rep(NA,(num_iter+1))
    nfr[[(k-1)*1000 + i]][1] <- length(reefdata$BCU_ID[reefdata$BCU_ID > 0])
    for(j in 2 : (num_iter+1)) {
      nfr[[(k-1)*1000 + i]][j] <- length(df[[((j-1)+(num_iter*(i-1)))]]$BCU_ID[df[[((j-1)+(num_iter*(i-1)))]]$BCU_ID > 0])
    }
  }
}
save(nfr, file = "/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/nfr_raw.RData")
print("ending 50reefs for loop 2")
print(paste0("nfr[[1]][2] =",nfr[[1]][2]))

png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Num_50r_random_lines_200k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nfr, pch = 20, xlab = "% of reefs remaining", ylab = "Number of 50 reefs", main = "3.3% screened", xlim = rev(range(nrf_remaining)), type = "n")
for(i in 1:(random_iter*200)){
  lines(x = nrf_remaining,y = nfr[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_nfr, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()


plotdata <- data.frame(x= nrf_remaining, y = m_nfr, lower = m_nfr - sd_nfr, upper =  m_nfr + sd_nfr)

png(file="/home/agreiner/scratch/worldconn_randomrep/1kreps_redux/outputs/Num_50r_random_eb_200k.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  #geom_line(aes(y=plotdata[,5],x=x),color="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Number of 50 reefs")+
  ggtitle("3.3% screened")+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

dev.off()


