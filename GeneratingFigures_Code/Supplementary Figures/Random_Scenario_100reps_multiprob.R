library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(ggplot2)
library(maps)

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
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA


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

#load in the centroids, NOTE: NEED TO FIX THE LONGITUDE IN THESE ONES (>180 ones need to have -360 subtracted)
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/centroids.RData')
#head(centroids) 
reefdata_full <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/RecalculatingResultswithNewScores_2.2020/scaled_scores_latlong_reproj_SJ.csv")
head(reefdata_full)
dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)]
head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")


#remove the rows corresponding to the grid cells that have no reef SA
centroids_reduced <- centroids[-id.noReef,]
#dim(centroids_reduced) #12292 
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
#length(reefdata$BCU_ID[reefdata$BCU_ID > 0]) #3334, 3334*(18*18) = 1,080,216km^2
#3334/434 = ~8
#so if want to stop when the reef SA = 50 reefs SA, should stop at 20 iterations 
num_iter <- 20

final_avg_sost_random_list <- list()
final_nclo_random_list <- list()
final_avgsizecluster_random_list <- list()
final_avgnd_random_list <- list()

load("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/AllReseedingresults_norandom.RData") #for the final graph
plateaupoint_rr_list <- list()
numreseeded_rr_list <- list()
mean_numclust_rr <- rep(NA,27)
mean_gmmean_rr <- rep(NA,27)
mean_largest_clust_rr <- rep(NA,27)
mean_avgss_rr <- rep(NA,27)
mean_avgnd_rr <- rep(NA,27)


#apparently stopped at k = 24 but don't see anything beyond k = 21 done?
for(k in 22:27){
  print(paste("k =",k))
  setwd("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject")
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
  
setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Random"))
getwd()
#remove from lowest (worst) to highest (best) score, in batches of num_remove (434)
#calculating the metrics as we go
#ord <- order(reefdata$score) #lowest to highest
set.seed(2) #moved inside for loop so that every new k experiences the same random scenarios
random_iter <- 100
final_iterator <- 0 #IDEA: just keep adding to the original things and then take the mean and std. dev. of each set (1, 20+1, 40+1 representing one step)
fr <- so_st <- cl <- dg <- df <- list() #initialize blank lists, each iteration adding in a new dataframe to each list
reefdata1 <- reefdata[,c(2,3,4,5,7,8,9,22)] #each time, just going to remove from the complete dataframe

setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Random/Reseeding"))

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
    fr[[i]] <- length(reefdata1$BCU_ID[reefdata1$BCU_ID > 0][-ord[1:final_reef]])
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
save(fr,so_st,cl,dg,df, file = paste0("RandomScenario_",probthresh,"_alldata_100reps.RData"))
#load(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Random/RandomScenario_",probthresh,"_alldata_100reps.RData")) #100 reps

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
png(file="Source_strength_random_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_sost, pch = 20, xlab = "% of reefs remaining", ylab = "Source Strength", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

final_avg_sost_random_list[[k]] <- m_sost[num_iter+1]

sost <- list()
for(i in 1:random_iter){
  sost[[i]] <- rep(NA,(num_iter+1))
  sost[[i]][1] <- mean(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))
  for(j in 2 : (num_iter+1)) {
    sost[[i]][j] <-mean(so_st[[((j-1)+(num_iter*(i-1)))]])
  }
}

png(file="Source_strength_random_lines_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_sost, pch = 20, xlab = "% of reefs remaining", ylab = "Source strength", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,100))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = sost[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_sost, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

raw_string <- rep(NA,random_iter)
for(i in 1:random_iter){
  raw_string[i] <- paste0("raw",i)
}
plotdata <- data.frame(x= nrf_remaining, y = m_sost, lower = m_sost - sd_sost, upper =  m_sost + sd_sost)
#the two lines below work but are not needed
#plotdata <- data.frame(x= nrf_remaining, y = m_nfr, lower = m_nfr - sd_nfr, upper =  m_nfr + sd_nfr, raw = nfr)
#names(plotdata)[5:(random_iter+4)] <- as.character(raw_string)

png(file="Source_strength_random_eb_100.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  #geom_line(aes(y=plotdata[,5],x=x),color="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Source Strength")+
  ggtitle(paste0(probthresh,"% screened"))+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

dev.off()
#Plot how 50 reefs are lost across removal 'time'
nrf_remov <- sd_nfr <- m_nfr <- s_nfr <- rep(NA,(num_iter+1))
m_nfr[1] <- length(reefdata$BCU_ID[reefdata$BCU_ID > 0])
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_nfr[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_nfr[i] <- m_nfr[i] + (length(df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID[df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID > 0]))/random_iter
  }
  s_nfr[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_nfr[i] <- s_nfr[i] + ((length(df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID[df[[((i-1)+(num_iter*(j-1)))]]$BCU_ID > 0]) - m_nfr[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_nfr[i] <- sqrt(s_nfr[i]/(random_iter-1))
}

nrf_remaining <- (12292 - nrf_remov)/12292 * 100
png(file="Num_50r_random_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nfr, pch = 20, xlab = "% of reefs remaining", ylab = "Number of 50 reefs", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

nfr <- list()
for(i in 1:random_iter){
  nfr[[i]] <- rep(NA,(num_iter+1))
  nfr[[i]][1] <- length(reefdata$BCU_ID[reefdata$BCU_ID > 0])
  for(j in 2 : (num_iter+1)) {
    nfr[[i]][j] <- length(df[[((j-1)+(num_iter*(i-1)))]]$BCU_ID[df[[((j-1)+(num_iter*(i-1)))]]$BCU_ID > 0])
  }
}

png(file="Num_50r_random_lines_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nfr, pch = 20, xlab = "% of reefs remaining", ylab = "Number of 50 reefs", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining)), type = "n")
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = nfr[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_nfr, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

raw_string <- rep(NA,random_iter)
for(i in 1:random_iter){
  raw_string[i] <- paste0("raw",i)
}
plotdata <- data.frame(x= nrf_remaining, y = m_nfr, lower = m_nfr - sd_nfr, upper =  m_nfr + sd_nfr)
#the two lines below work but are not needed
#plotdata <- data.frame(x= nrf_remaining, y = m_nfr, lower = m_nfr - sd_nfr, upper =  m_nfr + sd_nfr, raw = nfr)
#names(plotdata)[5:(random_iter+4)] <- as.character(raw_string)

png(file="Num_50r_random_eb_100.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  #geom_line(aes(y=plotdata[,5],x=x),color="blue",alpha=0.2)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Number of 50 reefs")+
  ggtitle(paste0(probthresh,"% screened"))+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

dev.off()
print("number of clusters time")
#Plot how number of clusters changes across removal 'time'
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
png(file="Num_clusters_random_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nclo, pch = 20, xlab = "% of reefs remaining", ylab = "Number of clusters", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()
#png(file="Num_clusters_b_last.png",width=12,height=10,units="cm",res=300)
#plot(x = nrf_remaining,y = nclo, xlab = "Number of reefs remaining", ylab = "#Clusters", main = "Number of Clusters", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
#dev.off()

final_nclo_random_list[[k]] <- m_nclo[num_iter+1]

nclo <- list()
for(i in 1:random_iter){
  nclo[[i]] <- rep(NA,(num_iter+1))
  nclo[[i]][1] <- clusters(g_orig, mode = "weak")$no
  for(j in 2 : (num_iter+1)) {
    nclo[[i]][j] <- cl[[((j-1)+(num_iter*(i-1)))]]$no
  }
}

png(file="Num_clusters_random_lines_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_nclo, pch = 20, xlab = "% of reefs remaining", ylab = "Number of clusters", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining)), type = "n", ylim = c(150,1450))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = nclo[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_nclo, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_nclo, lower = m_nclo - sd_nclo, upper =  m_nclo + sd_nclo)
png(file="Num_clusters_random_eb_100.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Number of clusters")+
  ggtitle(paste0(probthresh,"% screened"))+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()

#Plot how the size of the largest cluster changes across removal 'time'
nrf_remov <- s_largest_cluster <- m_largest_cluster <- sd_largest_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
m_largest_cluster[1] <- max(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
for(i in 2:(num_iter+1)){
  nrf_remov[i] <- ((i-1)*num_remove)
  m_largest_cluster[i] <- 0
  for(j in 1:random_iter){ #calculate the mean
    m_largest_cluster[i] <- m_largest_cluster[i] + (max(cl[[((i-1)+(num_iter*(j-1)))]]$csize))/random_iter
  }
  s_largest_cluster[i] <- 0
  for(j in 1:random_iter){ #calculate  the standard deviation, first calc summation((x_i - x_bar)^2)
    s_largest_cluster[i] <- s_largest_cluster[i] + (((max(cl[[((i-1)+(num_iter*(j-1)))]]$csize)) - m_largest_cluster[i])^2)
  }
  #then calc sqrt(summation/(n-1))
  sd_largest_cluster[i] <- sqrt(s_largest_cluster[i]/(random_iter-1))
}


nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
png(file="Size_clusters_random_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_largest_cluster, xlab = "% of reefs remaining", ylab = "Size of Largest Cluster", main = paste0("Size of Largest Cluster, ",probthresh,"% screened"), xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

largest_cluster <- list()
for(i in 1:random_iter){
  largest_cluster[[i]] <- rep(NA,(num_iter+1))
  largest_cluster[[i]][1] <- max(clusters(g_orig, mode = "weak")$csize)
  for(j in 2 : (num_iter+1)) {
    largest_cluster[[i]][j] <- (max(cl[[((j-1)+(num_iter*(i-1)))]]$csize))
  }
}

png(file="Size_clusters_random_lines_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_largest_cluster, pch = 20, xlab = "% of reefs remaining", ylab = "Size of Largest Cluster", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,6000))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = largest_cluster[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_largest_cluster, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_largest_cluster, lower = m_largest_cluster - sd_largest_cluster, upper =  m_largest_cluster + sd_largest_cluster)
png(file="Size_clusters_random_eb_100.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Cluster Size")+
  ggtitle(paste0(probthresh,"% screened"))+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()

#NEW on 7.26.2020
#change to geometric mean bc Marie-Josee says that that is better
#NEW on 7.15.2020
#Plot how the average size of the clusters changes across removal 'time'
print("avg cluster size time")
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
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
png(file="GMean_clusters_random_1k.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avgsize_cluster, xlab = "% of reefs remaining", ylab = "Avg Cluster Size", main = paste0("Avg Cluster Size, ",probthresh,"% screened"), xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
dev.off()

final_avgsizecluster_random_list[[k]] <- m_avgsize_cluster[num_iter+1]

avgsize_cluster <- list()
for(i in 1:random_iter){
  avgsize_cluster[[i]] <- rep(NA,(num_iter+1))
  avgsize_cluster[[i]][1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
  for(j in 2 : (num_iter+1)) {
    avgsize_cluster[[i]][j] <- (gm_mean(cl[[((j-1)+(num_iter*(i-1)))]]$csize))
  }
}

png(file="GMean_clusters_random_lines_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avgsize_cluster, pch = 20, xlab = "% of reefs remaining", ylab = "Avg Cluster Size", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining)), type = "n", ylim = c(0,5))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = avgsize_cluster[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_avgsize_cluster, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_avgsize_cluster, lower = m_avgsize_cluster - sd_avgsize_cluster, upper =  m_avgsize_cluster + sd_avgsize_cluster)
png(file="GMean_clusters_random_eb_100.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Cluster Size")+
  ggtitle(paste0(probthresh,"% screened"))+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()#Plot how the avg size of the clusters changes across removal 'time'

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
png(file="Node_degree_random_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,pch = 20, y = m_avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", xlim = rev(range(nrf_remaining)),main = paste0(probthresh,"% screened"))
dev.off()

final_avgnd_random_list[[k]] <- m_avg_degree[num_iter+1]

avg_degree <- list()
for(i in 1:random_iter){
  avg_degree[[i]] <- rep(NA,(num_iter+1))
  avg_degree[[i]][1] <- mean(degree(g_orig))
  for(j in 2 : (num_iter+1)) {
    avg_degree[[i]][j] <- (mean(dg[[((j-1)+(num_iter*(i-1)))]]))
  }
}

png(file="Node_degree_random_lines_100.png",width=12,height=10,units="cm",res=300)
plot(x = nrf_remaining,y = m_avg_degree, pch = 20, xlab = "% of reefs remaining", ylab = "Avg Node Degree", main = paste0(probthresh,"% screened"), xlim = rev(range(nrf_remaining)), type = "n", ylim = c(3,9))
for(i in 1:random_iter){
  lines(x = nrf_remaining,y = avg_degree[[i]],col = alpha("blue", 0.2), xlim = rev(range(nrf_remaining)))
}
points(x = nrf_remaining, y = m_avg_degree, xlim = rev(range(nrf_remaining)), pch = 20)
dev.off()

plotdata <- data.frame(x= nrf_remaining, y = m_avg_degree, lower = m_avg_degree - sd_avg_degree, upper =  m_avg_degree + sd_avg_degree)
png(file="Node_degree_random_eb_100.png",width=12,height=10,units="cm",res=300)
ggplot(plotdata)+
  geom_point(aes(y=y,x=x),pch=20)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=x), alpha = 0.3)+
  scale_x_reverse(name = "% of reefs remaining")+
  ylab("Avg Node Degree")+
  ggtitle(paste0(probthresh,"% screened"))+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
dev.off()

print("final map")
# Map the final map
world <- readOGR(dsn = "~/Dropbox/University of Toronto/Research Related/Code from Marco 6.2018 mediterranean larval connectivity/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")
nrf_remov <-  rep(NA,(num_iter+1))
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
  nrf_remov[i] <- ((i-1)*num_remove)
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#how many clusters are there, max 
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)


#assigning colours
#present-day connectivity (added 5.18.2021)
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
#set.seed(3) #added 4.13.2021, 3 seems better than 2
cols <- sample(rainbow(clusters(g_orig, mode = "weak")$no)) #might not work for all bc num clusters might be too small idk
#colouring some of them specific colours #this was taken from the 3.3% full graph, might be weird here but whatever
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

#plot the last iteration of the first set
print("final maps section")
#final maps - removed reefs coloured grey
library(scales)
df[[num_iter]]$cl.m <- cl[[num_iter]]$membership
df[[num_iter]]$cl.c <- cols[df[[num_iter]]$cl.m]
png(file=paste0("Final_Predicted_PercentReefs_Remaining_",probthresh,"screened_random.png"),width=10,height=10,units="cm",res=1000)
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected, reefdata$Latitude, col = "#eef3ff", pch = 20, cex = 0.1) #alpha("light blue",0.2) #0.05 didn't work for the random map
points(df[[num_iter]]$Longitude_corrected,df[[num_iter]]$Latitude,col=df[[num_iter]]$cl.c,pch=20,cex=0.2)
dev.off()

load("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/AllReseedingresults_norandom.RData") #for the final graph

print("starting reseeding")
#RESEEDING FROM RANDOM SCENARIO
setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Random/Reseeding"))
#these were generated above, but to save them into the reseeding folder...
save(plateaupoint_rr, numreseededrrlist, id.totalnotseeded, file = paste0("RandomScenarioReseeding_",probthresh,"_alldata_100reps.RData"))


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
  
  
  png(paste0("Realized_CU_Random_Reseed_",numsteps,"steps_mean",avgdist,"km.png"))
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
  
  # for(j in 1:random_iter){
  #   reefdata2$final_abund_rr <- 1
  #   reefdata2$final_abund_rr[id.totalnotseeded[[j]]] <- 0
  #   head(reefdata2$final_abund_rr)
  #   
  #   #reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)
  #   for(i in 1:(clusters(g_orig, mode = "weak")$no)){
  #     reseededclusters$num_reseeded_rr[i] <- sum(reefdata2$final_abund_rr[reefdata2$cl.m == i])
  #     reseededclusters$per_reseeded_rr[i] <- (reseededclusters$num_reseeded_rr[i]/reseededclusters$clustersize[i])*100
  #   }
  #   perreseededrr[[j]] <- reseededclusters$per_reseeded_rr
  # }
  # 
  # perreseededrr_mean <- rep(NA,(numsteps+1))
  # perreseededrr_sd <- rep(NA,(numsteps+1))
  # set <- rep(NA,random_iter)
  # for(i in 1:(numsteps+1)){
  #   for(j in 1:random_iter){
  #     set[j] <- perreseededrr[[j]][i]
  #   }
  #   perreseededrr_mean[i] <- mean(set)
  #   perreseededrr_sd[i] <- sd(set)
  #   set <- rep(NA,random_iter)	
  # }
  # 
  # png(paste0("Realized_CU_MeanRandom_HistPerReseeded_",numsteps,"steps_mean",avgdist,"km.png"))
  # hist(perreseededrr_mean, main = "Mean Percentage of Each Original Cluster that was Re-seeded", xlab = "Percent Re-seeded")
  # dev.off()
  
  #NEW 7.26.2020
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
  #GOAL GRAPH
  #NEW 7.27.2020 - adding thin lines around the random scenario one
  #plotting all of the numreseeeded from the 4 different 'scenarios' on one graph - total seeded as a PERCENTTAGE OF TOTAL REEFS
  #12292 - final_reef = number started with in the predicted and catastrophe situations
  initialreefs = 12292 - final_reef
  png(paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_totalpercentages_allrandom_100reps_no50r_mean",avgdist,"km.png"))
  plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "#1d457f", type = 'l', ylim = c(20,100), lwd = 3) #blue
  lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_cr + initialreefs)/12292)*100), col  = "#f9ad2a", lwd = 3) #darkgoldenrod1
  for(i in 1:(random_iter)){
    lines(x = seq(1,(numsteps+1),by=1),y = (((numreseededrrlist[[i]] + initialreefs)/12292)*100),col = alpha("#e6b6c0", 0.2),lwd=1) ##c7c0d1
  }
  lines(x=seq(1,(numsteps+1),by=1),y=(((numreseededrr_mean + initialreefs)/12292)*100), col  = "#cc5c76",lty =3, lwd = 3) ##7f43d1
  legend("bottomright", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)", "Random Scenario"), col = c("#1d457f","#f9ad2a","#cc5c76","#e6b6c0"), lty = c(1,1,3,1), lwd = c(3,3,3,1))
  dev.off()

}
#bc it didn't finish running and i had to close my computer and stop it...i lost 1-21 of these bc i was dumb and didn't save them first (also whenever i try to stop it...it just quits) so it won't save the whole thing when this saves (just 22-27)
save(plateaupoint_rr_list, numreseeded_rr_list, mean_numclust_rr, mean_gmmean_rr, mean_largest_clust_rr, mean_avgss_rr, mean_avgnd_rr, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/FullRandom_multiprob_reseeding.RData")

save(final_avg_sost_random_list, final_nclo_random_list, final_avgsizecluster_random_list, final_avgnd_random_list, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/FullRandom_multiprob_removal.RData")