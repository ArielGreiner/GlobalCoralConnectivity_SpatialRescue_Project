library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(ggplot2)
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


#load connectivity matrix
load('OriginalDataConversion/ConvertingWoodetal2014Data/connmat_reduced.RData')

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

world <- readOGR(dsn = "GeneratingFigures_Code/AdditionalFiles/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")


##Calculating where to stop the destruction
#total SA: 12292*324km^2 = 3,982,608km
#12292/434 = ~28
#SA of BCU reefs (50 reefs)
#length(reefdata$BCU_ID[reefdata$BCU_ID > 0]) #3334, 3334*(18*18) = 1,080,216km^2
#3334/434 = ~8
#so if want to stop when the reef SA = 50 reefs SA, should stop at 20 iterations 
num_iter <- 20

final_avg_sost_random_list <- list()
final_nclo_random_list <- list()
final_avgsizecluster_random_list <- list()
final_avgnd_random_list <- list()

plateaupoint_rr_list <- list()
numreseeded_rr_list <- list()
mean_numclust_rr <- rep(NA,27)
mean_gmmean_rr <- rep(NA,27)
mean_avgss_rr <- rep(NA,27)
mean_avgnd_rr <- rep(NA,27)


for(k in 1:27){
  print(paste("k =",k))
  if(k == 1){
    avgdist <- 10
    probthresh <- 19.5
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 2){
    avgdist <- 20
    probthresh <- 9
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 3){
    avgdist <- 30
    probthresh <- 5.3
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 4){
    avgdist <- 40
    probthresh <- 3.5
    connmat_reduced <- connmat_reduced_list[[k]]
  }  
  
  if(k == 5){
    avgdist <- 50
    probthresh <- 2.6
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 6){
    avgdist <- 60
    probthresh <- 2
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 7){
    avgdist <- 70
    probthresh <- 1.6
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 8){
    avgdist <- 80
    probthresh <- 1.3
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 9){
    avgdist <- 90
    probthresh <- 1.07
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 10){
    avgdist <- 100
    probthresh <- 0.89
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 11){
    avgdist <- 110
    probthresh <- 0.7526
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 12){
    avgdist <- 120
    probthresh <- 0.6451
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 13){
    avgdist <- 130
    probthresh <- 0.6451
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 14){
    avgdist <- 140
    probthresh <- 0.4732
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 15){
    avgdist <- 150
    probthresh <- 0.4
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 16){
    avgdist <- 160
    probthresh <- 0.35
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 17){
    avgdist <- 170
    probthresh <- 0.3011
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 18){
    avgdist <- 180
    probthresh <- 0.26
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 19){
    avgdist <- 190
    probthresh <- 0.2259
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 20){
    avgdist <- 200
    probthresh <- 0.2
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 21){
    avgdist <- 210
    probthresh <- 0.1721
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 22){
    avgdist <- 110
    probthresh <- 0.1506
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 23){
    avgdist <- 110
    probthresh <- 0.1299
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 24){
    avgdist <- 240
    probthresh <- 0.1182
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 25){
    avgdist <- 250
    probthresh <- 0.1
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 26){
    avgdist <- 260
    probthresh <- 0.09
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
  if(k == 27){
    avgdist <- 42
    probthresh <- 3.3
    connmat_reduced <- connmat_reduced_list[[k]]
  }
  
#remove from lowest (worst) to highest (best) score, in batches of num_remove (434)
#calculating the metrics as we go
#ord <- order(reefdata$score) #lowest to highest
set.seed(2) #moved inside for loop so that every new k experiences the same random scenarios
random_iter <- 100
final_iterator <- 0 #IDEA: just keep adding to the original things and then take the mean and std. dev. of each set (1, 20+1, 40+1 representing one step)
so_st <- cl <- dg <- df <- list() #initialize blank lists, each iteration adding in a new dataframe to each list
reefdata1 <- reefdata[,c(2,3,4,5,7,8,9,22)] #each time, just going to remove from the complete dataframe

abundmatlist <- starterreefslist <- numreseededrrlist <- id.totalnotseeded <- perreseededrr <- list()
numsteps <- 50
num_remove <- 434
num_iter <- 20
random_iter <- 100
final_reef <- (num_iter*num_remove)
plateaupoint <- 500
plateaupoint_rr <- rep(NA,random_iter)
reefdata$num <- seq(1,12292,by=1)

i=1
for(j in 1:random_iter){
  print(paste("j = ", j))
  ord <- sample(seq(1,length(reefdata$score),by=1)) #highest to lowest
  #i goes to 20 - won't get down to 0 reefs, but that's what we want
  num_remove <- 434 #remove 434 the first time, 434 the next run (total of 868 removed), etc
  print(paste("i = ", i))
  for(i in (final_iterator + 1) : (num_iter+final_iterator)) {
    print(paste("i = ",i, "final_i =", final_iterator))
    final_reef <- ((i-final_iterator)*num_remove)  
    connmat_i <- connmat_reduced[-ord[1:final_reef],-ord[1:final_reef]]
    df[[i]] <- reefdata1[-ord[1:final_reef],]
    g <- graph.adjacency(as.matrix(connmat_i), weighted = TRUE)
    cl[[i]] <- clusters(g,mode="weak") #number of clusters, size of each cluster, membership of each cluster
    dg[[i]] <- degree(g) #getting node degree
    so_st[[i]] <- ego_size(graph = g, order = 40, nodes = V(g), mode = "out")
    if(i == num_iter+final_iterator){
      final_iterator <- num_iter+final_iterator
      i = 0
    }
  }
  #reseeding
    #ord <- sample(seq(1,length(reefdata$score),by=1))
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
i=0  
}
plateaupoint_rr_list[[k]] <- plateaupoint_rr
numreseeded_rr_list[[k]] <- numreseededrrlist
#save(so_st,cl,dg,df, file = paste0("RandomScenario_",probthresh,"_alldata_100reps.RData"))

print("saved data, onto source strength")
#Plot how average source strength changes across removal 'time'
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
nrf_remov <- sd_sost <- m_sost <- s_sost <- rep(NA,(num_iter+1))
m_sost[1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_sost[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_sost[i] <- m_sost[i] + (mean(so_st[[((i-1)+(num_iter*(j-1)))]])/random_iter)
  }
  s_sost[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_sost[i] <- s_sost[i] + ((mean(so_st[[((i-1)+(num_iter*(j-1)))]]) - m_sost[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_sost[i] <- sqrt(s_sost[i]/(random_iter-1))
}

nrf_remaining <- (12292 - nrf_remov)/12292 * 100


final_avg_sost_random_list[[k]] <- m_sost[num_iter+1]

sost <- list()
for(i in 1:random_iter){
  sost[[i]] <- rep(NA,(num_iter+1))
  sost[[i]][1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
  for(j in 2 : (num_iter+1)) {
    sost[[i]][j] <-mean(so_st[[((j-1)+(num_iter*(i-1)))]])
  }
}



print("number of networks time")
#Plot how number of networks changes across removal 'time'
nrf_remov <- m_nclo <- sd_nclo <- s_nclo <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_nclo[1] <- clusters(g_orig, mode = "weak")$no
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_nclo[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_nclo[i] <- m_nclo[i] + (cl[[((i-1)+(num_iter*(j-1)))]]$no)/random_iter
  }
  s_nclo[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_nclo[i] <- s_nclo[i] + (((cl[[((i-1)+(num_iter*(j-1)))]]$no) - m_nclo[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_nclo[i] <- sqrt(s_nclo[i]/(random_iter-1))
}

nrf_remaining <- (12292 - nrf_remov)/12292 * 100

final_nclo_random_list[[k]] <- m_nclo[num_iter+1]

nclo <- list()
for(i in 1:random_iter){
  nclo[[i]] <- rep(NA,(num_iter+1))
  nclo[[i]][1] <- clusters(g_orig, mode = "weak")$no
  for(j in 2 : (num_iter+1)) {
    nclo[[i]][j] <- cl[[((j-1)+(num_iter*(i-1)))]]$no
  }
}




#Plot how the geometric mean size of the network changes across removal 'time'
print("avg network size time")
nrf_remov <- s_avgsize_network <- m_avgsize_network <- sd_avgsize_network <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_avgsize_network[1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_avgsize_network[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_avgsize_network[i] <- m_avgsize_network[i] + (gm_mean(cl[[((i-1)+(num_iter*(j-1)))]]$csize))/random_iter
  }
  s_avgsize_network[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_avgsize_network[i] <- s_avgsize_network[i] + (((gm_mean(cl[[((i-1)+(num_iter*(j-1)))]]$csize)) - m_avgsize_network[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_avgsize_network[i] <- sqrt(s_avgsize_network[i]/(random_iter-1))
}


nrf_remaining <- (12292 - nrf_remov)/12292 * 100


final_avgsizecluster_random_list[[k]] <- m_avgsize_network[num_iter+1]

avgsize_network <- list()
for(i in 1:random_iter){
  avgsize_network[[i]] <- rep(NA,(num_iter+1))
  avgsize_network[[i]][1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
  for(j in 2 : (num_iter+1)) {
    avgsize_network[[i]][j] <- (gm_mean(cl[[((j-1)+(num_iter*(i-1)))]]$csize))
  }
}


print("avg node degree")
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


final_avgnd_random_list[[k]] <- m_avg_degree[num_iter+1]

avg_degree <- list()
for(i in 1:random_iter){
  avg_degree[[i]] <- rep(NA,(num_iter+1))
  avg_degree[[i]][1] <- mean(degree(g_orig))
  for(j in 2 : (num_iter+1)) {
    avg_degree[[i]][j] <- (mean(dg[[((j-1)+(num_iter*(i-1)))]]))
  }
}


print("starting reseeding")
#RESEEDING FROM RANDOM SCENARIO
#note: commented this section out because didn't put these folders into the Github folder, but useful to save periodically when running this code so going to leave this here (commented out) as an example of how best to do this
#setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Random/Reseeding"))
#these were generated above, but to save them into the reseeding folder...
#save(plateaupoint_rr, numreseededrrlist, id.totalnotseeded, file = paste0("RandomScenarioReseeding_",probthresh,"_alldata_100reps.RData"))


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
  
  
  
  
  #final connectivity metrics
  print("beginning final connectivity metrics")
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
    avgdeg[i] <- mean(dg_Random[[i]], na.rm = T) 
  }
  #mean(numclust) #231.937
  #mean(meanclust) #6.036148
  #mean(maxclust) #3042.69
  #mean(avgsost) #105.0341
  #mean(avgdeg, na.rm = T) #NaN (UGH)
  
  mean_numclust_rr[k] <- mean(numclust)
  mean_gmmean_rr[k] <- mean(meanclust)
  mean_largest_clust_rr[k] <- mean(maxclust)
  mean_avgss_rr[k] <- mean(avgsost)
  mean_avgnd_rr[k] <- mean(avgdeg, na.rm = T)
  
  abundance_mat <- NA
  plateaupoint <- 500
  
  
  numreseeded_pr <- numreseeded_pr_list[[k]]
  numreseeded_cr <- numreseeded_cr_list[[k]]


}

save(plateaupoint_rr_list, numreseeded_rr_list, mean_numclust_rr, mean_gmmean_rr, mean_largest_clust_rr, mean_avgss_rr, mean_avgnd_rr, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/FullRandom_multiprob_reseeding.RData")

save(final_avg_sost_random_list, final_nclo_random_list, final_avgsizecluster_random_list, final_avgnd_random_list, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/FullRandom_multiprob_removal.RData")