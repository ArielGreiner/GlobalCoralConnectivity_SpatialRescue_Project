library(rgdal)

#read 50 reefs data in (from Marco)
fiftyreefsdata <- readOGR("~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/DatafromMarco/reef_cells_allstat/reef_cells_allstat.shp")
#it's a spatial data frame
#str(fiftyreefsdata@data)
#order of scores: thermal stress history metrics (13), thermal stress projection metrics (8*19), recent thermal conditions (4), connectivity (2), cyclones (3)
orig_data <- fiftyreefsdata@data
head(orig_data)

#load('~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/reefs_R_dataframe.rdata')
#head(reefs.df) #this is the data that Hawthorne sent, his IDs are different and it's not a shape file (will have to deal with that later, or wait for him to send me the geospatial data bc ill need to have the shape file eventually? can i insert this into the shape file above? not sure) - I need to use this dataframe and not the one from Marco to re-calculate the scores because the 50reefs people took >47000 unique combinations - each time taking one particular metric from each theme...with each unique combination they summed over all 38 weighting schemes and then summed over all unique combinations
#new data from Hawthorne
reefs.df <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/FixedData_1.20.2020/50Reefs_data/tabular/50Reefs_dataset.csv")
head(reefs.df)

dim(reefs.df)
#write.table(reefs.df,file=paste0("reefsdf.csv"),sep=",")
reefs.df.stdized <- unname(reefs.df) #removes the column names
reefs.df.stdized <- as.matrix(reefs.df.stdized) #these two rows were to avoid the for loops below and be able to just use Row or Col sum but it didn't work for some reason
for(i in 5:178){
	reefs.df.stdized[,i] <- (reefs.df.stdized[,i] - mean(reefs.df.stdized[,i]))/sd(reefs.df.stdized[,i]) #'mean 0 and unit variance standardization'
	if(sd(reefs.df.stdized[,i]) == 0){
		break;
	}
	reefs.df.stdized[,i] <- -reefs.df.stdized[,i] #'sign changed such that higher values are more desirable' #did it so higher values are bad and that didn't work either
}

#order of scores: thermal stress history metrics (13, col 5:17), thermal stress projection metrics (8*19), recent thermal conditions (4), connectivity (2), cyclones (3)
#names(reefs.df)[5:17] #thermal stress history metrics
tsh_met <- as.matrix(reefs.df.stdized[,5:17])
#names(reefs.df)[18:169] #thermal stress projection metrics
tsp_met <- as.matrix(reefs.df.stdized[,18:169])
#names(reefs.df)[170:173] #recent thermal conditions
rtc_met <- as.matrix(reefs.df.stdized[,170:173])
#names(reefs.df)[174:175] #connectivity
conn_met <- as.matrix(-reefs.df.stdized[,174:175]) #because had the right sign to begin with, higher values = better for corals
#names(reefs.df)[176:178] #cyclones
cyc_met <- as.matrix(-reefs.df.stdized[,176:178]) #this one apparently flips too..

#1.28.2020 - calculating scores without connectivity data
tshnum <- 13
tspnum <- 8*19
rtcnum <- 4
#connnum <- 2
cycnum <- 3
#library(matrixStats)
#tsh_num <- rowSums[tsh_met] <- this isn't working for some reason
tsh_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	tsh_num[i] <- sum(tsh_met[i,])
}
tsp_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	tsp_num[i] <- sum(tsp_met[i,])
}
rtc_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	rtc_num[i] <- sum(rtc_met[i,])
}
#conn_num <- c(NA,dim(reefs.df)[1])
#for(i in 1:dim(reefs.df)[1]){
#	conn_num[i] <- sum(conn_met[i,])
#}
cyc_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	cyc_num[i] <- sum(cyc_met[i,])
}
#weights<- read.csv("~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/metric_weights.csv")
#new data:
weights<-read.csv("~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/FixedData_1.20.2020/50Reefs_data/tabular/metric_weights.csv")
reweighted <- weights[,c(1,2,3,5)]
rescore_onescheme_dadm <- rep(NA,dim(reweighted)[1])
rescore_onemet_dadm <- rep(NA, 23712)
#redo_ogsum_onemet <- matrix(data = NA, nrow= 47424, ncol = dim(reefs.df)[1])
rescore_unscaledscore_dadm <- rep(NA, dim(reefs.df.stdized)[1])
for(i in 1:dim(reefs.df)[1]){
print(paste("i = ",i))
	for(j in 1:dim(reweighted)[1]){
rescore_onemet_dadm <- as.numeric(reweighted[j,1]*(tspnum*rtcnum*cycnum)*tsh_num[i] + reweighted[j,2]*(tshnum*rtcnum*cycnum)*tsp_num[i] + reweighted[j,3]*(tshnum*tspnum*cycnum)*rtc_num[i] + reweighted[j,4]*(tshnum*tspnum*rtcnum)*cyc_num[i])
rescore_onescheme_dadm[j] <- sum(rescore_onemet_dadm)
}
rescore_unscaledscore_dadm[i] <- sum(rescore_onescheme_dadm)
}
rescore_scaledscore_dadm <- rescore_unscaledscore_dadm/max(abs(rescore_unscaledscore_dadm))
head(rescore_scaledscore_dadm)
head(reefs.df[,4]) #shouldn't be identical, should be similar
#hist(rescore_scaledscore_dadm) #histograms look fairly similar
save(rescore_scaledscore_dadm,file="scaled_scores_1.28.2020.RData")
scaledscoreswlatlong <- data.frame(longitude = reefs.df[,2],laitude = reefs.df[,3],rescore = rescore_scaledscore_dadm)
write.table(scaledscoreswlatlong,file="scaled_scores_latlong_1.28.2020.csv",sep=",",row.names=FALSE)

oldscoreswlatlong <- data.frame(longitude = reefs.df[,2],laitude = reefs.df[,3],score = reefs.df[,4])
write.table(oldscoreswlatlong,file="old_latlong_2.5.2020.csv",sep=",",row.names=FALSE)




#12.26.2019 - add up all of the numbers in each theme for each reef, multiply each of the themes by the number of values in each of the other themes and then multiply that by w
tshnum <- 13
tspnum <- 8*19
rtcnum <- 4
connnum <- 2
cycnum <- 3
#library(matrixStats)
#tsh_num <- rowSums[tsh_met] <- this isn't working for some reason
tsh_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	tsh_num[i] <- sum(tsh_met[i,])
}
tsp_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	tsp_num[i] <- sum(tsp_met[i,])
}
rtc_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	rtc_num[i] <- sum(rtc_met[i,])
}
conn_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	conn_num[i] <- sum(conn_met[i,])
}
cyc_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	cyc_num[i] <- sum(cyc_met[i,])
}
#weights<- read.csv("~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/metric_weights.csv")
#new data:
weights<-read.csv("~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/FixedData_1.20.2020/50Reefs_data/tabular/metric_weights.csv")
redo_ogsum_onescheme_dadm <- rep(NA,dim(weights)[1])
redo_ogsum_onemet_dadm <- rep(NA, 47424)
#redo_ogsum_onemet <- matrix(data = NA, nrow= 47424, ncol = dim(reefs.df)[1])
redo_ogunscaledscore_dadm <- rep(NA, dim(reefs.df.stdized)[1])
for(i in 1:dim(reefs.df)[1]){
print(paste("i = ",i))
	for(j in 1:dim(weights)[1]){
redo_ogsum_onemet_dadm <- as.numeric(weights[j,1]*(tspnum*rtcnum*connnum*cycnum)*tsh_num[i] + weights[j,2]*(tshnum*rtcnum*connnum*cycnum)*tsp_num[i] + weights[j,3]*(tshnum*tspnum*connnum*cycnum)*rtc_num[i] + weights[j,4]*(tshnum*tspnum*rtcnum*cycnum)*conn_num[i]  + weights[j,5]*(tshnum*tspnum*rtcnum*connnum)*cyc_num[i])
redo_ogsum_onescheme_dadm[j] <- sum(redo_ogsum_onemet_dadm)
}
redo_ogunscaledscore_dadm[i] <- sum(redo_ogsum_onescheme_dadm)
}
redo_ogscaledscore_dadm <- redo_ogunscaledscore_dadm/max(abs(redo_ogunscaledscore_dadm))
head(redo_ogscaledscore_dadm)
head(reefs.df[,4])
sum(round(redo_ogscaledscore_dadm,10) - round(reefs.df[,4],10)) #9.999999e-11, okay so basically identical
save(redo_ogscaledscore_dadm,file="scaled_ogscores_1.20.2020.RData")


#make dataframe with all of the different possible unique combos (should be 47424)
recalc_ogscores_perms <- data.frame(num = 1:47424, tsh_met = rep(c(1:13),each = 8*19*4*2*3), tsp_met = rep(c(1:(8*19)), each = 4*2*3), rtc_met = rep(c(1:4), each = 2*3), conn_met = rep(c(1,2), each = 3), cyc_met = rep(c(1,2,3)))
#recalc_ogscores_perms[1:100,]
scores_perms <- data.frame(num = 1:23712, tsh_met = rep(c(1:13),each = 8*19*4*3), tsp_met = rep(c(1:(8*19)), each = 4*3), rtc_met = rep(c(1:4), each = 3), cyc_met = rep(c(1,2,3)))

#load in weights
weights<- read.csv("~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/metric_weights.csv")
head(weights)
names(weights) <- c("thermalhistoryweights", "thermalprojectionweights", "recentthermalweights", "connectivityweights", "cyclonesweights")
reweighted <- weights
#remove the connectivity column
reweighted <- weights[,c(1,2,3,5)]

#note: final scores were scaled so to see if this formula works - need to recalculate them, then scale them the same way to see if those are the same

#checking to see whether I understand how the scores were calculated 
#most recent version
library(matrixStats)
#redo_ogsum_onescheme <- matrix(data = NA, nrow= dim(weights)[1], ncol = dim(reefs.df)[1])
redo_ogsum_onescheme <- rep(NA,dim(weights)[1])
redo_ogsum_onemet <- rep(NA, 47424)
#redo_ogsum_onemet <- matrix(data = NA, nrow= 47424, ncol = dim(reefs.df)[1])
redo_ogunscaledscore <- rep(NA, dim(reefs.df.stdized)[1])

for(i in 1:dim(reefs.df)[1]){ #for every reef
  print(paste("i = ",i))
  for(j in 1:dim(weights)[1]){ #across every weighting scheme #whole loop: 0.229 seconds 
redo_ogsum_onemet <- as.numeric(weights[j,1]*tsh_met[i,recalc_ogscores_perms[,2]] + weights[j,2]*tsp_met[i,recalc_ogscores_perms[,3]] + weights[j,3]*rtc_met[i,recalc_ogscores_perms[,4]] + weights[j,4]*conn_met[i,recalc_ogscores_perms[,5]]  + weights[j,5]*cyc_met[i,recalc_ogscores_perms[,6]])
redo_ogsum_onescheme[j] <- sum(redo_ogsum_onemet)
  }
  redo_ogunscaledscore[i] <- sum(redo_ogsum_onescheme)
}

save(redo_ogsum_onescheme,redo_ogunscaledscore,file="recalcscores_inprep_12.11.2019.RData")
redo_ogscaledscore <- redo_ogunscaledscore/max(abs(redo_ogunscaledscore)) #same as before, not honestly sure what im doing wrong?
save(redo_ogscaledscore,file="scaled_ogscores_12.11.2019.RData")


#long version
for(i in 1:dim(reefs.df)[1]){ #for every reef
  time.p <- proc.time()
  for(j in 1:dim(weights)[1]){ #across every weighting scheme
  print(paste("j = ",j))
 for(k in 1:dim(recalc_ogscores_perms)[1]){ #for each metric combo
   redo_ogsum_onemet[k] <- as.numeric(weights[j,1]*tsh_met[i,recalc_ogscores_perms[k,2]] + weights[j,2]*tsp_met[i,recalc_ogscores_perms[k,3]] + weights[j,3]*rtc_met[i,recalc_ogscores_perms[k,4]] + weights[j,4]*conn_met[i,recalc_ogscores_perms[k,5]]  + weights[j,5]*cyc_met[i,recalc_ogscores_perms[k,6]])
 }
  redo_ogsum_onescheme[j] <- sum(redo_ogsum_onemet)
  }
time.p <- proc.time()-time.p
redo_ogunscaledscore[i] <- sum(redo_ogsum_onescheme)
}


#this is wrong (okay well this got the same answer as the one above)
ogsum_onescheme <- matrix(data = NA, nrow= dim(weights)[1], ncol = dim(reefs.df)[1])
ogunscaledscore <- rep(NA, dim(orig_data)[1])

ogrstsh_met <- rowSums2(tsh_met[,recalc_ogscores_perms[,2]]) #Error: vector memory exhausted (limit reached) <- for all of them :(
#used the RStudio trick here and that seemed to work? https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
ogrstsp_met <- rowSums2(tsp_met[,recalc_ogscores_perms[,3]])
ogrsrtc_met <- rowSums2(rtc_met[,recalc_ogscores_perms[,4]])
ogrsconn_met <- rowSums2(conn_met[,recalc_ogscores_perms[,5]])
ogrscyc_met <- rowSums2(cyc_met[,recalc_ogscores_perms[,6]])

for(j in 1:dim(weights)[1]){
	print(paste("j = ",j))	
		#time.p <- proc.time()
		ogsum_onescheme[j,] <- weights[j,1]*ogrstsh_met + weights[j,2]*ogrstsp_met + weights[j,3]*ogrsrtc_met + weights[j,4]*ogrsconn_met  + weights[j,5]*ogrscyc_met 
		#time.p <- proc.time()-time.p				
}
ogunscaledscore <- colSums2(ogsum_onescheme) 
save(ogsum_onescheme,ogunscaledscore,file="recalcscores_inprep.RData")
ogscaledscore <- ogunscaledscore/max(abs(ogunscaledscore))
save(ogscaledscore,file="scaled_ogscores.RData")
load('~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/recalculating og scores_round1/scaled_ogscores.RData')
#hist(ogscaledscore) #looks like the re-scaled score hist and not like the original scores
#hist(reefs.df$score) #looks like before 
save(ogrstsh_met,ogrstsp_met,ogrsconn_met,ogrscyc_met,ogrstsp_met, file = "forcalc_ogscores.RData")
#slower versioon
sum_onescheme <- rep(NA, dim(weights)[1])
sum_onemet <- rep(NA,(13*8*19*4*2*3))
unscaledscore <- rep(NA, dim(orig_data)[1])

for(i in 1:dim(orig_data)[1]){
	for(j in 1:dim(weights)[1]){
	print(paste("j = ",j))	
		for(k in 1:(13*8*19*4*2*3)){
		sum_onemet[k] <- weights[j,1]*tsh_met[,recalc_ogscores_perms[k,2]] + weights[j,2]*tsp_met[,recalc_ogscores_perms[k,3]] + weights[j,3]*tsp_met[,recalc_ogscores_perms[k,4]] + weights[j,3]*tsp_met[,recalc_ogscores_perms[k,4]] + weights[j,4]*tsp_met[,recalc_ogscores_perms[k,5]] + weights[j,5]*tsp_met[,recalc_ogscores_perms[k,6]]
		}			
	sum_onescheme[j] <- sum(sum_onemet)	
	}
	unscaledscore[i] <- sum(sum_onescheme)
}
scaledscore <- unscaledscore/max(unscaledscore)
save(scaledscore,file="scaled_ogscores.RData")

#re-calculating the scores
#re-standardize the weights so they still add up to 1
reweighted$sum <- rowSums(reweighted[,1:4])
head(reweighted)
reweighted[,1] <- reweighted[,1]/reweighted$sum
reweighted[,2] <- reweighted[,2]/reweighted$sum
reweighted[,3] <- reweighted[,3]/reweighted$sum
reweighted[,4] <- reweighted[,4]/reweighted$sum
head(reweighted)

#shortest version - NEED TO RECALC
library(matrixStats)
sum_onescheme <- matrix(data = NA, nrow= dim(weights)[1], ncol = dim(reefs.df)[1])
#sum_onemet <- matrix(data = NA, nrow= (13*8*19*4*3), ncol = dim(reefs.df)[1])
unscaledscore <- rep(NA, dim(orig_data)[1])
rstsh_met <- rowSums2(tsh_met[,scores_perms[,2]]) #basically the sum of a big addition of things = a big addition of sums
rstsp_met <- rowSums2(tsp_met[,scores_perms[,3]])
rsrtc_met <- rowSums2(rtc_met[,scores_perms[,4]])
rscyc_met <- rowSums2(cyc_met[,scores_perms[,5]])

for(j in 1:dim(weights)[1]){
	print(paste("j = ",j))	
		time.p <- proc.time()
		sum_onescheme[j,] <- reweighted[j,1]*rstsh_met + reweighted[j,2]*rstsp_met + reweighted[j,3]*rsrtc_met + reweighted[j,4]*rscyc_met 
		time.p <- proc.time()-time.p				
}
unscaledscore <- colSums2(sum_onescheme) 
save(sum_onescheme,unscaledscore,file="newscores_inprep.RData")

scaledscore <- unscaledscore/max(abs(unscaledscore))
save(scaledscore,file="scaled_newscores.RData")
load('~/Dropbox/University of Toronto/Research Related/Hawthorne Beyer Data/rescaled_scores_v1_11.2019/scaled_newscores.RData')
hist(scaledscore)
reefs_newscores <- data.frame(ID = reefs.df$ID, longitude = reefs.df$longitude, latitude = reefs.df$latitude, score = scaledscore)
write.csv(scaledscore, file = "scaled_newscores.csv")
write.csv(reefs_newscores, file = "newscores_reefs.csv")

#longer version
sum_onescheme <- matrix(data = NA, nrow= dim(weights)[1], ncol = dim(reefs.df)[1])
sum_onemet <- matrix(data = NA, nrow= (13*8*19*4*3), ncol = dim(reefs.df)[1])
unscaledscore <- rep(NA, dim(orig_data)[1])

for(j in 1:dim(weights)[1]){
	print(paste("j = ",j))	
	for(k in 1:(13*8*19*4*3)){
		#time.p <- proc.time()
		sum_onemet[k,] <- reweighted[j,1]*tsh_met[,scores_perms[k,2]] + reweighted[j,2]*tsp_met[,scores_perms[k,3]] + reweighted[j,3]*rtc_met[,scores_perms[k,4]] + reweighted[j,4]*cyc_met[,scores_perms[k,5]]  
		#time.p <- proc.time()-time.p
	}			
	sum_onescheme[j,] <- colSums(sum_onemet)	 
}
unscaledscore <- colSums(sum_onescheme) 
save(sum_onescheme,unscaledscore,file="newscores_inprep.RData")

scaledscore <- unscaledscore/max(unscaledscore)
save(scaledscore,file="scaled_newscores.RData")

#old
sum_onescheme <- rep(NA, dim(weights)[1])
sum_onemet <- rep(NA,(13*8*19*4*3))
unscaledscore <- rep(NA, dim(reefs.df)[1])
#i = 347,  j = 27, k = 6385 <- that corresponds with what's in 'newscores_inprep.RData
#should restart at i = 347 bc that one didn't finish
#unscaledscore[101:500]...they're all the same, that is definitely wrong...maybe just stop this completely then and re-visit later? after 6 days it got through only 347/15k so...it needs to be run on compute canada anyways

for(i in 1:dim(reefs.df)[1]){
	for(j in 1:dim(weights)[1]){
	print(paste("j = ",j))	
		for(k in 1:(13*8*19*4*3)){
		sum_onemet[k] <- reweighted[j,1]*tsh_met[,scores_perms[k,2]] + reweighted[j,2]*tsp_met[,scores_perms[k,3]] + reweighted[j,3]*rtc_met[,scores_perms[k,4]] + reweighted[j,4]*cyc_met[,scores_perms[k,5]]  
		}			
	sum_onescheme[j] <- sum(sum_onemet)	
	}
	unscaledscore[i] <- sum(sum_onescheme)
	save(sum_onescheme,unscaledscore,file="newscores_inprep.RData")
}
scaledscore <- unscaledscore/max(unscaledscore)
save(scaledscore,file="scaled_newscores.RData")
