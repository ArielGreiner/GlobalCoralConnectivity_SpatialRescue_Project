library(rgdal)

#REPLACE LINE BELOW WITH THE DIRECTORY THAT YOU PUT THE DATA FROM THE BEYER ET AL 2018 PAPER
reefs.df <- read.csv(file = "~/.../50Reefs_data/tabular/50Reefs_dataset.csv")
head(reefs.df)
dim(reefs.df)

#standardize the dataset
reefs.df.stdized <- unname(reefs.df) #removes the column names
reefs.df.stdized <- as.matrix(reefs.df.stdized) 
for(i in 5:178){
	reefs.df.stdized[,i] <- (reefs.df.stdized[,i] - mean(reefs.df.stdized[,i]))/sd(reefs.df.stdized[,i]) #'mean 0 and unit variance standardization'
	if(sd(reefs.df.stdized[,i]) == 0){
		break;
	}
	reefs.df.stdized[,i] <- -reefs.df.stdized[,i] #'sign changed such that higher values are more desirable' 
}

#separate the columns representing the different metrics into their own matrices
#order of scores: thermal stress history metrics (13, col 5:17), thermal stress projection metrics (8*19), recent thermal conditions (4), connectivity (2), cyclones (3)
#names(reefs.df)[5:17] #thermal stress history metrics
tsh_met <- as.matrix(reefs.df.stdized[,5:17])
#names(reefs.df)[18:169] #thermal stress projection metrics
tsp_met <- as.matrix(reefs.df.stdized[,18:169])
#names(reefs.df)[170:173] #recent thermal conditions
rtc_met <- as.matrix(reefs.df.stdized[,170:173])
#names(reefs.df)[174:175] #connectivity
conn_met <- as.matrix(-reefs.df.stdized[,174:175]) 
#names(reefs.df)[176:178] #cyclones
cyc_met <- as.matrix(-reefs.df.stdized[,176:178]) 

#calculating scores without connectivity data
tshnum <- 13
tspnum <- 8*19
rtcnum <- 4
#connnum <- 2
cycnum <- 3

#take the sum of each row, across all of the columns of each metric (speeds up some later steps)
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
cyc_num <- c(NA,dim(reefs.df)[1])
for(i in 1:dim(reefs.df)[1]){
	cyc_num[i] <- sum(cyc_met[i,])
}

#calculate the re-weighted score, after the connectivity elements have been removed!

#REPLACE LINE BELOW WITH THE DIRECTORY THAT YOU PUT THE DATA FROM THE BEYER ET AL 2018 PAPER
weights<-read.csv("~/.../50Reefs_data/tabular/metric_weights.csv")
reweighted <- weights[,c(1,2,3,5)]
rescore_onescheme <- rep(NA,dim(reweighted)[1])
rescore_onemet <- rep(NA, 23712)
rescore_unscaledscore <- rep(NA, dim(reefs.df.stdized)[1])
for(i in 1:dim(reefs.df)[1]){
print(paste("i = ",i))
	for(j in 1:dim(reweighted)[1]){
rescore_onemet <- as.numeric(reweighted[j,1]*(tspnum*rtcnum*cycnum)*tsh_num[i] + reweighted[j,2]*(tshnum*rtcnum*cycnum)*tsp_num[i] + reweighted[j,3]*(tshnum*tspnum*cycnum)*rtc_num[i] + reweighted[j,4]*(tshnum*tspnum*rtcnum)*cyc_num[i])
rescore_onescheme[j] <- sum(rescore_onemet)
}
rescore_unscaledscore[i] <- sum(rescore_onescheme)
}
rescore_scaledscore <- rescore_unscaledscore/max(abs(rescore_unscaledscore))
head(rescore_scaledscore)
save(rescore_scaledscore,file="~/Github/GlobalCoralConnectivity_SpatialRescue_Project/OriginalDataConversion/GeneratingScoresfrom50ReefsScores/scaled_scores.RData")
scaledscoreswlatlong <- data.frame(longitude = reefs.df[,2],laitude = reefs.df[,3],rescore = rescore_scaledscore_dadm)
write.table(scaledscoreswlatlong,file="~/Github/GlobalCoralConnectivity_SpatialRescue_Project/OriginalDataConversion/GeneratingScoresfrom50ReefsScores/scaled_scores_latlong.csv",sep=",",row.names=FALSE)

