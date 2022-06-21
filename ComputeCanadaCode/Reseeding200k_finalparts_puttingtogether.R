setwd("/Users/arielgreiner/Documents/computecanada_randomremovalmorereps_june2020/Plotting200kreps_finalreseeding_outputs")
#numclust,meanclust,maxclust,avgsost,avgdeg

#load in i=1->50000 versions, rename as *_small
load("/Users/arielgreiner/Documents/computecanada_randomremovalmorereps_june2020/Plotting200kreps_finalreseeding_outputs/reseeding/meanfinalconnectivitymetrics_all_50000.RData")
numclust_small <- numclust
meanclust_small <- meanclust
maxclust_small <- maxclust
avgsost_small <- avgsost
avgdeg_small <- avgdeg 

#load in i=50001->end versions
load("/Users/arielgreiner/Documents/computecanada_randomremovalmorereps_june2020/Plotting200kreps_finalreseeding_outputs/reseeding/meanfinalconnectivitymetrics_all_2e+05.RData")

for(i in 1:50000){
  print(paste("i = ",i))
  numclust[i] <- numclust_small[i]
  meanclust[i] <- meanclust_small[i]
  maxclust[i] <- maxclust_small[i] 
  avgsost[i] <- avgsost_small[i]
  avgdeg[i] <- avgdeg_small[i]
}

mean_numclust <- mean(numclust) 
mean_gmeanclust <- mean(meanclust) 
mean_maxclust <- mean(maxclust)
mean_avgsost <- mean(avgsost) 
mean_avgdeg <- mean(avgdeg, na.rm = T)

save(mean_numclust,mean_gmeanclust,mean_maxclust,mean_avgsost,mean_avgdeg, file = "fixingmanually/meanfinalconnectivitymetrics_all_done_FULL.RData")