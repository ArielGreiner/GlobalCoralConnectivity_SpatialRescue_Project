library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
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


#load connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA

#what do i need to threshold it by for the average dispersal distance to be 10km? 100km?

#matrix of euclidean distances as calculated from the original Sally Wood connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/Distances/Euclid.dist.mat.12292.RData')
#histogram of euclidean distances 
hist(euclid.dist.mat@x,breaks=seq(0,8000,40))

# Check if the two matrices are in the same order (they should because one has been built from the other!)
all(connmat_reduced@i == euclid.dist.mat@i)
all(connmat_reduced@p == euclid.dist.mat@p)
#yes :)

#Method 1: Choose a probability threshold for the connectivity matrix
euclid_dist_threshold <- euclid.dist.mat@x
#probabilities negatively correlated with distances (lower probabilities associated with larger distances)
quantile(connmat_reduced@x, c(0.05,0.25,0.5,0.75,0.9))
#          5%          25%          50%          75%          90% 
#0.0001075269 0.0002150538 0.0006451613 0.0032258065 0.0121505376 
connmat_reduced_threshold <- connmat_reduced@x
mean(euclid_dist_threshold) #529.5538
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0001)]) #529.5538 (makes sense, since this threshold is so low)
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0008)]) #270.4602
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0009)]) #260.1317 #260!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001)]) #251.4931 #250
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00101)]) #251.4931
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00109)]) #243.3224
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0011)]) #243.3224
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00115)]) #243.3224
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00117)]) #243.3224
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00118)]) #243.3224 
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001182)]) #243.3224 #240
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001183)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001184)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001185)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00119)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0012)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00125)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00129)]) #236.2528
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001299)]) #229.9041 #230
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0013)]) #229.9041
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0015)]) #223.9827
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001505)]) #223.9827 
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001506)]) #218.6798 #220
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001507)]) #218.6798
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00151)]) #218.6798
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00152)]) #218.6798
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00155)]) #218.6798
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0016)]) #218.6798
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0017)]) #213.6867
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00171)]) #213.6867
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00172)]) #213.6867
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001721)]) #209.1456 #210
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.001725)]) #209.1456
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00173)]) #209.1456
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00175)]) #209.1456
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0018)]) #209.1456
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.002)]) #200.8497 #200
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0021)]) #197.0314
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00215)]) #197.0314
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0022)]) #193.4183
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00225)]) #193.4183
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.002255)]) #193.4183
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.002257)]) #193.4183
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.002258)]) #193.4183
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.002259)]) #189.8827 #190
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00226)]) #189.8827
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00227)]) #189.8827
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00228)]) #189.8827
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0023)]) #189.8827
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0025)]) #183.4353
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0026)]) #180.48 #180!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0028)]) #174.9512
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.003)]) #172.2349
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00301)]) #172.2349
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.003011)]) #169.5408 #170!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.003015)]) #169.5408
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00302)]) #169.5408
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00305)]) #169.5408
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0031)]) #169.5408
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0033)]) #164.8342
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0035)]) #160.4076 #160!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.004)]) #150.5634 #150!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0045)]) #144.0264
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0047)]) #141.013
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00471)]) #141.013
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.004715)]) #141.013
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00473)]) #141.013
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.004731)]) #141.013
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.004732)]) #139.4852 #140
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.004733)]) #139.4852
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.004735)]) #139.4852
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00474)]) #139.4852
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00475)]) #139.4852
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0048)]) #139.4852
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.005)]) #136.6488
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0055)]) #130.3428 #130!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.006)]) #125.3663
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0061)]) #124.2506
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0063)]) #121.9913
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0064)]) #120.8857
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00645)]) #120.8857
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.006451)]) #120.8857 #120!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.006452)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.006453)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.006455)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00646)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00648)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00649)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0065)]) #119.7784
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.007)]) #114.5726
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0075)]) #110.7877
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00751)]) #110.7877
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00752)]) #110.7877
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.007525)]) #110.7877
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.007526)]) #110.7877 #110
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.007527)]) #109.9403
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.007529)]) #109.9403
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00753)]) #109.9403
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.00755)]) #109.9403
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0076)]) #109.9403
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.008)]) #106.4399
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0085)]) #102.34
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0088)]) #100.8792
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0089)]) #100.1995
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.009)]) #99.49532
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0095)]) #96.16046
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.01)]) #93.3437
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0105)]) #91.11556
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.0107)]) #90.09244 #90!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.011)]) #88.55332
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.012)]) #84.11251
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.013)]) #80.31571 #80!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.014)]) #76.18465
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.015)]) #73.07008
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.016)]) #70.23352 #70!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.018)]) #64.8124
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.02)]) #60.36905 #60!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.025)]) #51.71971
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.026)]) #50.37229 #50!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.027)]) #48.93655
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.028)]) #47.59688
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.03)]) #45.3572
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.031)]) #44.18484
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.032)]) #43.18716
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.033)]) #42.27354
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.034)]) #41.30342
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.035)]) #40.43335 #40!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.037)]) #38.72758
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.04)]) #36.362
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.05)]) #31.31909
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.053)]) #30.0356 #30!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.055)]) #29.29386
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.06)]) #27.73867 
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.07)]) #25
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.08)]) #22.71214
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.09)]) #20.49308 #20!
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.1)]) #18.57658
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.15)]) #13.69761
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.18)]) #11.18774
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.19)]) #10.40349
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.195)]) #10.0542
mean(euclid_dist_threshold[which(connmat_reduced_threshold > 0.2)]) #9.861753

#note: median(euclid_dist_threshold[which(connmat_reduced_threshold > 0.01345)]) = 40 (0.0134 is 46)

probs_touse <- which(connmat_reduced_threshold > 0.0089) #mean dispersal distance = 100km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 0.89%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 27940km
median(euclid_dist_threshold[probs_touse]) #56
length(euclid_dist_threshold[probs_touse]) #197039 connections included

probs_touse <- which(connmat_reduced_threshold > 0.0107) #mean dispersal distance = 90km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 1.07%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 2794km
median(euclid_dist_threshold[probs_touse]) #53
length(euclid_dist_threshold[probs_touse]) #169774 connections included

probs_touse <- which(connmat_reduced_threshold > 0.013) #mean dispersal distance = 80km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 1.3%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 2794km
median(euclid_dist_threshold[probs_touse]) #48
length(euclid_dist_threshold[probs_touse]) #144612 connections included

probs_touse <- which(connmat_reduced_threshold > 0.016) #mean dispersal distance = 70km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 1.6%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 1355km
median(euclid_dist_threshold[probs_touse]) #39
length(euclid_dist_threshold[probs_touse]) #119524 connections included

probs_touse <- which(connmat_reduced_threshold > 0.02) #mean dispersal distance = 60km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 2%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 1355km
median(euclid_dist_threshold[probs_touse]) #36
length(euclid_dist_threshold[probs_touse]) #95458 connections included

probs_touse <- which(connmat_reduced_threshold > 0.026) #mean dispersal distance = 50km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 2.6%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 1337km
median(euclid_dist_threshold[probs_touse]) #25
length(euclid_dist_threshold[probs_touse]) #72492 connections included

probs_touse <- which(connmat_reduced_threshold > 0.035) #mean dispersal distance = 40km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 3.5%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 1200km
median(euclid_dist_threshold[probs_touse]) #24
length(euclid_dist_threshold[probs_touse]) #51236 connections included

probs_touse <- which(connmat_reduced_threshold > 0.053) #mean dispersal distance = 30km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 5.3%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 906km
median(euclid_dist_threshold[probs_touse]) #18
length(euclid_dist_threshold[probs_touse]) #30115 connections included

probs_touse <- which(connmat_reduced_threshold > 0.09) #mean dispersal distance = 20km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 9%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 809km
median(euclid_dist_threshold[probs_touse]) #17
length(euclid_dist_threshold[probs_touse]) #13298 connections included

probs_touse <- which(connmat_reduced_threshold > 0.195) #mean dispersal distance = 10km
hist(euclid_dist_threshold[probs_touse], main = "Thresholded at 19.5%", xlab = "Distances (km)")
range(euclid_dist_threshold[probs_touse]) #0km -> 418km
median(euclid_dist_threshold[probs_touse]) #0
length(euclid_dist_threshold[probs_touse]) #2878 connections included

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

getwd()
# "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject"

final_avg_sost_pcs_list <- final_avg_sost_rl_list <- list()
final_nclo_pcs_list <- final_nclo_rl_list <- list()
final_avgsizecluster_pcs_list <- final_avgsizecluster_rl_list <- list()
final_avgnd_pcs_list <- final_avgnd_rl_list <- list()
for(j in 1:27){ 
setwd("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject")
if(j == 1){
  avgdist <- 10
  probthresh <- 19.5
  connmat_reduced <- connmat_reduced_list[[j]]
}
  
if(j == 2){
  avgdist <- 20
  probthresh <- 9
  connmat_reduced <- connmat_reduced_list[[j]]
}

if(j == 3){
  avgdist <- 30
  probthresh <- 5.3
  connmat_reduced <- connmat_reduced_list[[j]]
}
  
if(j == 4){
  avgdist <- 40
  probthresh <- 3.5
  connmat_reduced <- connmat_reduced_list[[j]]
}  

if(j == 5){
  avgdist <- 50
  probthresh <- 2.6
  connmat_reduced <- connmat_reduced_list[[j]]
}
  
  if(j == 6){
    avgdist <- 60
    probthresh <- 2
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 7){
    avgdist <- 70
    probthresh <- 1.6
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 8){
    avgdist <- 80
    probthresh <- 1.3
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 9){
    avgdist <- 90
    probthresh <- 1.07
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
if(j == 10){
  avgdist <- 100
  probthresh <- 0.89
  connmat_reduced <- connmat_reduced_list[[j]]
}
  
  if(j == 11){
    avgdist <- 110
    probthresh <- 0.7526
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 12){
    avgdist <- 120
    probthresh <- 0.6451
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 13){
    avgdist <- 130
    probthresh <- 0.6451
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 14){
    avgdist <- 140
    probthresh <- 0.4732
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 15){
    avgdist <- 150
    probthresh <- 0.4
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 16){
    avgdist <- 160
    probthresh <- 0.35
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 17){
    avgdist <- 170
    probthresh <- 0.3011
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 18){
    avgdist <- 180
    probthresh <- 0.26
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 19){
    avgdist <- 190
    probthresh <- 0.2259
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 20){
    avgdist <- 200
    probthresh <- 0.2
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 21){
    avgdist <- 210
    probthresh <- 0.1721
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 22){
    avgdist <- 110
    probthresh <- 0.1506
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 23){
    avgdist <- 110
    probthresh <- 0.1299
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 24){
    avgdist <- 240
    probthresh <- 0.1182
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 25){
    avgdist <- 250
    probthresh <- 0.1
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 26){
    avgdist <- 260
    probthresh <- 0.09
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  if(j == 27){
    avgdist <- 42
    probthresh <- 3.3
    connmat_reduced <- connmat_reduced_list[[j]]
  }
  
  
for(k in 1:2){
#k = 1: predicted scenario
#k = 2: refuge-loss scenario
  
if(k == 1){
ord <- order(reefdata$score) #lowest to highest (decreasing = FALSE is default)
setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/PCS"))
getwd()
}
  
if(k == 2){
ord <- order(reefdata$score, decreasing = TRUE) #refuge-loss
setwd(paste0("/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/",avgdist,"km/Refuge_Loss"))
getwd()
}

print(paste("k = ",k,"and j = ",j))

print("starting the removal")
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
print("source strength section")
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
if(k == 1){
  png(file="Source_strength_b_last.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = avg_sost, xlab = "% of reefs remaining", ylab = "Avg Source Strength", type = 'l', xlim = rev(range(nrf_remaining)),main = paste0(probthresh, "% screened"), ylim =  c(0,10000)) #ylim =  c(0,100)
  dev.off()
  final_avg_sost_pcs_list[[j]] <- avg_sost[num_iter+1]
}

if(k == 2){
  png(file="Source_strength_b_first.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = avg_sost, xlab = "% of reefs remaining", ylab = "Avg Source Strength", type = 'l', xlim = rev(range(nrf_remaining)),main = paste0(probthresh, "% screened"), ylim =  c(0,10000)) #ylim =  c(0,100)
  dev.off()
  final_avg_sost_rl_list[[j]] <- avg_sost[num_iter+1]
}

#Global Distribution of Source Strength across removal 'time'
 # maxval <- range(ego_size(graph = g_orig, order = 40, nodes = V(g_orig), mode = "out"))[2]
 # cols_source <- tim.colors(maxval+1)
 # if(k == 1){
 # for(i in 1:num_iter){
 #   df[[i]]$sost_col <- cols_source[so_st[[i]]]
 #   png(file=paste0("GlobalSourceStrength_",round(nrf_remaining[i+1],2),"PercentReefsRemaining_b_last.png"),width=20,height=20,units="cm",res=1500) #res=1000
 #   plot(world, main = paste("Global Distribution of Source Strength",round(nrf_remaining[i+1],2),"% Reefs Remaining"),lwd = 0.000002, col = "light grey")
 #   points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$sost_col,pch=20,cex=0.2)
 #   legend("topleft", c("Out-degree = 1", "Out-degree = 350", "Out-degree = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(20,20,20))
 #   dev.off()
 #   
 #   png(file=paste0("GlobalSourceStrength_",round(nrf_remaining[i+1],2),"PercentReefsRemaining_b_last_small.png"),width=20,height=20,units="cm",res=1500) #res=1000
 #   plot(world, main = paste("Global Distribution of Source Strength",round(nrf_remaining[i+1],2),"% Reefs Remaining"),lwd = 0.000002, col = "light grey")
 #   points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$sost_col,pch=15,cex=0.03)
 #   legend("topleft", c("Out-degree = 1", "Out-degree = 350", "Out-degree = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(15,15,15))
 #   dev.off()
 # }
 # }
 # 
 # if(k == 2){
 #   for(i in 1:num_iter){
 #     df[[i]]$sost_col <- cols_source[so_st[[i]]]
 #     png(file=paste0("GlobalSourceStrength_",round(nrf_remaining[i+1],2),"PercentReefsRemaining_b_first.png"),width=20,height=20,units="cm",res=1500) #res=1000
 #     plot(world, main = paste("Global Distribution of Source Strength",round(nrf_remaining[i+1],2),"% Reefs Remaining"),lwd = 0.000002, col = "light grey")
 #     points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$sost_col,pch=20,cex=0.2)
 #     legend("topleft", c("Out-degree = 1", "Out-degree = 350", "Out-degree = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(20,20,20))
 #     dev.off()
 #     
 #     png(file=paste0("GlobalSourceStrength_",round(nrf_remaining[i+1],2),"PercentReefsRemaining_b_first_small.png"),width=20,height=20,units="cm",res=1500) #res=1000
 #     plot(world, main = paste("Global Distribution of Source Strength",round(nrf_remaining[i+1],2),"% Reefs Remaining"),lwd = 0.000002, col = "light grey")
 #     points(df[[i]]$Longitude,df[[i]]$Latitude,col=df[[i]]$sost_col,pch=15,cex=0.03)
 #     legend("topleft", c("Out-degree = 1", "Out-degree = 350", "Out-degree = 709"), col = c(cols_source[1],cols_source[350],cols_source[709]), pch = c(15,15,15))
 #     dev.off()
 #   }
 # }

#checking source strength distribution at the last step
range(so_st[[20]]) #1 to 190
#hist(so_st[[20]])
#for(i in 1:20){
#print(median(so_st[[i]])) 
#}
#quantile(so_st[[20]], c(0.25,0.5,0.75)) #6  16  42 
#boxplot(so_st[[20]])

print("50 reefs section")
#Plot how 50 reefs are lost across removal 'time'
nrf_remov <- nfr <- rep(NA,(num_iter+1))
nfr[1] <- length(reefdata$BCU_ID[reefdata$BCU_ID > 0])
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
  nrf_remov[i] <- ((i-1)*num_remove)
  nfr[i] <- length(df[[i-1]]$BCU_ID[df[[i-1]]$BCU_ID > 0])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
if(k == 1){
  png(file="Num_50r_b_last.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = nfr, xlab = "% of reefs remaining", ylab = "Number of 50 reefs", type = 'l', main = paste0(probthresh, "% screened"), xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
}

if(k == 2){
  png(file="Num_50r_b_first.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = nfr, xlab = "% of reefs remaining", ylab = "Number of 50 reefs", type = 'l', main = paste0(probthresh, "% screened"), xlim = rev(range(nrf_remaining))) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
}

print("entering number of cluster section")
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
if(k == 1){
  png(file="Num_clusters_b_last.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = nclo, xlab = "% of reefs remaining", ylab = "Number of clusters", type = 'l', main = paste0(probthresh, "% screened"), xlim = rev(range(nrf_remaining)), ylim = c(0,12250)) #ylim = c(150,1450)
  #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
  #png(file="Num_clusters_b_last.png",width=12,height=10,units="cm",res=300)
  #plot(x = nrf_remaining,y = nclo, xlab = "Number of reefs remaining", ylab = "#Clusters", main = "Number of Clusters", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  #dev.off()
  
  final_nclo_pcs_list[[j]] <- nclo[num_iter+1]
}

if(k == 2){
  png(file="Num_clusters_b_first.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = nclo, xlab = "% of reefs remaining", ylab = "Number of clusters", type = 'l', main = paste0(probthresh, "% screened"), xlim = rev(range(nrf_remaining)), ylim = c(0,12250)) #ylim = c(150,1450)
  #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
  #png(file="Num_clusters_b_first.png",width=12,height=10,units="cm",res=300)
  #plot(x = nrf_remaining,y = nclo, xlab = "Number of reefs remaining", ylab = "#Clusters", main = "Number of Clusters", xlim = rev(range(nrf_remaining)), pch = 20) #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  #dev.off()
  final_nclo_rl_list[[j]] <- nclo[num_iter+1]
}

print("Average cluster size section")
#NEW on 7.26.2020 - geometric_mean version
#NEW on 7.15.2020 - mean version
#Plot how the average size of the clusters changes across removal 'time'
nrf_remov <- avgsize_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
avgsize_cluster[1] <- gm_mean(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
print("average cluster size...")
for(i in 2 : (num_iter+1)) {
  nrf_remov[i] <- ((i-1)*num_remove)
  avgsize_cluster[i] <- gm_mean(cl[[i-1]]$csize)
  print(avgsize_cluster[i])
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
if(k == 1){
  png(file="GMeanSize_clusters_b_last.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = avgsize_cluster, xlab = "% of reefs remaining", ylab = "Geom_Mean Cluster Size", main = paste0("Geom_Mean Cluster Size, ", probthresh, "% screened"), xlim = rev(range(nrf_remaining)), type = 'l', ylim = c(0,8)) #ylim = c(0,5)
  #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
  final_avgsizecluster_pcs_list[[j]] <- avgsize_cluster[num_iter+1]
}

if(k == 2){
  png(file="GMeanSize_clusters_b_first.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = avgsize_cluster, xlab = "% of reefs remaining", ylab = "Geom_Mean Cluster Size", main = paste0("Geom_Mean Cluster Size, ", probthresh, "% screened"), xlim = rev(range(nrf_remaining)), type = 'l', ylim = c(0,8)) #ylim = c(0,5)
  #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
  final_avgsizecluster_rl_list[[j]] <- avgsize_cluster[num_iter+1]
}

print("Largest cluster section")
#Plot how the size of the largest cluster changes across removal 'time'
nrf_remov <- largest_cluster <- rep(NA,(num_iter+1))
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
largest_cluster[1] <- max(clusters(g_orig, mode = "weak")$csize)
nrf_remov[1] <- 0
for(i in 2 : (num_iter+1)) {
  nrf_remov[i] <- ((i-1)*num_remove)
  largest_cluster[i] <- max(cl[[i-1]]$csize)
}
nrf_remaining <- (12292 - nrf_remov)/12292 * 100
#plot(x = nrf_remaining,y = largest_cluster, xlab = "Number of reefs remaining", ylab = "Largest cluster", main = "Over 40km screened", xlim = rev(range(nrf_remaining)))
if(k==1){
  png(file="Size_clusters_b_last.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = largest_cluster, xlab = "% of reefs remaining", ylab = "Size of Largest Cluster", main = paste0("Size of Largest Cluster, ", probthresh, "% screened"), xlim = rev(range(nrf_remaining)), type = 'l') #ylim = c(0,6000)
  #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
}
if(k==2){
  png(file="Size_clusters_b_first.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = largest_cluster, xlab = "% of reefs remaining", ylab = "Size of Largest Cluster", main = paste0("Size of Largest Cluster, ", probthresh, "% screened"), xlim = rev(range(nrf_remaining)), type = 'l') #ylim = c(0,6000)
  #had to do this weird xlim thing to get the x-axis to go from largest -> smallest
  dev.off()
}

print("Average node degree section")
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
if(k==1){
  png(file="Node_degree_b_last.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", type = 'l', xlim = rev(range(nrf_remaining)),main = paste0(probthresh, "% screened"), ylim = c(0,35)) #ylim = c(3,10)
  dev.off()
  final_avgnd_pcs_list[[j]] <- avg_degree[num_iter+1]
}

if(k==2){
  png(file="Node_degree_b_first.png",width=12,height=10,units="cm",res=300)
  plot(x = nrf_remaining,y = avg_degree, xlab = "% of reefs remaining", ylab = "Avg Node Degree", type = 'l', xlim = rev(range(nrf_remaining)),main = paste0(probthresh, "% screened"), ylim = c(0,35)) #ylim = c(3,10)
  dev.off()
  final_avgnd_rl_list[[j]] <- avg_degree[num_iter+1]
}

print("map cluster section")
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
print("number of clusters...")
for(i in 1:num_iter){
  print(cl[[i]]$no) #626 -> ...
}

#present-day connectivity (added 5.18.2021)
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
set.seed(3) #added 4.13.2021, 3 seems better than 2
cols <- sample(rainbow(clusters(g_orig, mode = "weak")$no)) #might not work for all bc num clusters might be too small idk
#colouring some of them specific colours #might not work well for all of them but whatever
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
reefdata$cl.c <- cols[reefdata$cl.m]

png(paste0("PCentre_",probthresh,"_Allreefs_small_othercolours.png"),width=20,height=20,units="cm",res=1500) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=15,cex=0.03) #cex=0.05
dev.off()

png(paste0("PCentre_",probthresh,"_Allreefs_othercolours.png"),width=20,height=20,units="cm",res=1500, pointsize=4) #res=1000
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected,reefdata$Latitude,col=reefdata$cl.c,pch=20,cex=0.2) #cex=0.05
dev.off()

#set.seed(3) #added 4.22.2021
#5.17.2021: might need to change this because the number of clusters in g_orig might be too small
#cols <- sample(tim.colors(clusters(g_orig, mode = "weak")$no)) #the original has the most clusters, useful to keep the colours of the clusters consistent (to the degree that this does that)
#colouring some of them specific colours
#cols[166] <- "purple"
#cols[30] <- "pink"
#cols[1] <- "blue"
#cols[53] <- "green"
#cols[54] <- "red"
#cols[19] <- "black"

if(k==1){
#PCentred maps version
for(i in 1:num_iter){
  df[[i]]$cl.m <- cl[[i]]$membership
  #cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  #cols <- sample(tim.colors(cl[[i]]$no)) #defined above instead
  #df[[i]]$cl.c <- cols[df[[i]]$cl.m * 10]
  df[[i]]$cl.c <- cols[df[[i]]$cl.m]
  png(file=paste0(nrf_remaining[i+1],"PercentReefs_Remaining_",probthresh,"screened_b_last.png"),width=10,height=10,units="cm",res=1000) 
  plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
           fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
  points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
  dev.off()
}
#dev.off()

print("final maps section")
#final maps - removed reefs coloured grey
library(scales)
i=20
df[[i]]$cl.m <- cl[[i]]$membership
df[[i]]$cl.c <- cols[df[[i]]$cl.m]
png(file=paste0("Final_Predicted_PercentReefs_Remaining_",probthresh,"screened_b_last.png"),width=10,height=10,units="cm",res=1000)
plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
         fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
points(reefdata$Longitude_corrected, reefdata$Latitude, col = "#eef3ff", pch = 20, cex = 0.1) #alpha("light blue",0.2) #0.05 didn't work for the random map
points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
dev.off()
}

if(k==2){
  #PCentred maps version
  for(i in 1:num_iter){
    df[[i]]$cl.m <- cl[[i]]$membership
    #cols <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    #cols <- sample(tim.colors(cl[[i]]$no)) #defined above instead
    #df[[i]]$cl.c <- cols[df[[i]]$cl.m * 10]
    df[[i]]$cl.c <- cols[df[[i]]$cl.m]
    png(file=paste0(nrf_remaining[i+1],"PercentReefs_Remaining_",probthresh,"screened_b_first.png"),width=10,height=10,units="cm",res=1000) #NEED TO RERUN bc got rid of pointsize since last saved
    plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
             fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
    points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
    dev.off()
  }
  #dev.off()

  print("final maps section")
  #final maps - removed reefs coloured grey
  library(scales)
  i=20
  df[[i]]$cl.m <- cl[[i]]$membership
  df[[i]]$cl.c <- cols[df[[i]]$cl.m]
  png(file=paste0("Final_Predicted_PercentReefs_Remaining_",probthresh,"screened_b_first.png"),width=10,height=10,units="cm",res=1000)
  plot.map("world", center=newcentre, col="gainsboro",bg="white",lwd = 0.000002,
           fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0))
  points(reefdata$Longitude_corrected, reefdata$Latitude, col = "#eef3ff", pch = 20, cex = 0.1) #alpha("light blue",0.2) #0.05 didn't work for the random map
  points(df[[i]]$Longitude_corrected,df[[i]]$Latitude,col=df[[i]]$cl.c,pch=20,cex=0.2)
  dev.off()
}

}
}

save(final_avg_sost_pcs_list, final_avg_sost_rl_list, final_nclo_pcs_list, final_nclo_rl_list, final_avgsizecluster_pcs_list, final_avgsizecluster_rl_list, final_avgnd_pcs_list, final_avgnd_rl_list, file = "/Users/arielgreiner/GitHub/PhDThesisProjects/GlobalCoralConnectivityProject/SensitivityTesting/finalremovallists_norandom.RData")

