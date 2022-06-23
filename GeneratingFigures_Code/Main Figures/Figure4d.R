library(rgdal)
library(Matrix)
library(igraph)
library(fields)
library(scales)
library(colorRamps)
library(maps)
library(here)

here() #here is now set to your root directory with .Rproj file


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###REALIZED CONNECTIVITY###
#load in the shape file with all of the data
reefdata_full <- read.csv(file = "~/Dropbox/University of Toronto/Research Related/RecalculatingResultswithNewScores_2.2020/scaled_scores_latlong_reproj_SJ.csv")
head(reefdata_full)
dim(reefdata_full)
reefdata <- reefdata_full[,c(3:6,34,13:28)]
head(reefdata)
names(reefdata) <- c("TARGET_FID","PolyNo","Longitude", "Latitude", "score","reefpct","BCU_ID","BCU_Name","ReefArea","RfPctArea","BCUarea","BCUpctarea","Protectedarea","pct_Protectedarea","MarineRealm","MRealm_Name", "NumEEZ", "EEZName", "EEZ_Sovereign1","EEZ_Sovereign2","EEZ_Sovereign3")

load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/cells12292.RData') #this includes 'id.noReef' and 'b' (shapefile)
#reefdata_old <- b@data
#head(reefdata_old)

#load connectivity matrix
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/connmat_reduced.RData')
dim(connmat_reduced) #12292x12292 with >0 reef SA

connmat_thresholder <- connmat_reduced
#test <- matrix(c(1,3,2,4), nrow = 2, ncol = 2)
#test[test <= 2] <- 0
#     [,1] [,2]
#[1,]    0    0
#[2,]    3    4
connmat_thresholder[connmat_thresholder <= 0.033] <- 0
connmat_reduced <- connmat_thresholder

#what is the distribution of the scores?
#hist(reefdata$score) #centred on 0, -1 to 0.2

#correct the longitude
reefdata$Longitude[reefdata$Longitude > 180] <- reefdata$Longitude[reefdata$Longitude > 180] - 360

#for plotting pacific centred maps
reefdata$Longitude_corrected <- reefdata$Longitude
newcentre <- 180
#range(reefdata$Longitude_corrected)
reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] <- reefdata$Longitude_corrected[reefdata$Longitude > 0 & reefdata$Longitude  <  180] - newcentre

reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] <- reefdata$Longitude_corrected[reefdata$Longitude < 0 & reefdata$Longitude  >  -180] + newcentre
#reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] <- reefdata$Longitude_corrected[reefdata$Longitude_corrected > 180] - 360


world <- readOGR(dsn = "~/Dropbox/University of Toronto/Research Related/Code from Marco 6.2018 mediterranean larval connectivity/worldcountryshapeetc", layer = "ne_110m_admin_0_countries")

#note: plateaupoint_cr = 26, plateaupoint_50r = 23, plateaupoint_pr = 23 so final plateaupoint = 23
ultimateplateau <- 26
reefdata2 <- reefdata[,c(2,3,4)]
g_orig <- graph.adjacency(as.matrix(connmat_reduced), weighted = TRUE)
reefdata2$cl.m <- clusters(g_orig, mode = "weak")$membership
reseededclusters <- data.frame(clusterID <- seq(1,clusters(g_orig, mode = "weak")$no,by=1), clustersize <- clusters(g_orig, mode = "weak")$csize, num_reseeded_50r <- NA, per_reseeded_50r <- NA,num_reseeded_pr <- NA, per_reseeded_pr <- NA,num_reseeded_cr <- NA, per_reseeded_cr <- NA)

#REMOVED A BUNCH OF CODE THAT WAS UNNECESSARY FOR GOAL PLOT (12.15.2020)

#RESEEDING FROM 50 REEFS
#Start with just the 50 reefs
numsteps <-50
id.50reefs <- which(reefdata$BCU_ID > 0) #this IDs the rows with 50 reefs
length(id.50reefs) #3334
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = numsteps+1) #10 time steps
abundance_mat[,1][id.50reefs] <- 100 #time step 0: assign all the 50 reefs to an abundance of 100

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

stop <- "no"
#Plot the number that were re-seeded over time
numreseeded_50r <- rep(NA,(numsteps+1))
numreseeded_50r[1] <- 0
for(i in 2:(numsteps + 1)){
  #line below could be simplified to just be length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[id.50reefs,i]) but whatever, same for all of the other scenarios
  numreseeded_50r[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[id.50reefs,i][abundance_mat[id.50reefs,i]>0])	
  if((numreseeded_50r[i] <= numreseeded_50r[i-1]) & stop == "no"){ #added the 'or' condition bc of some of the other cases needing it
    plateaupoint <- i-1
    stop <- "yes"
  }
}
plateaupoint_50r <- plateaupoint


id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-id.50reefs]
}
#heatcols <- diverge_hcl(10)
#heatcols <- heat_hcl(12,h=c(0,-100), l = c(75,40), c = c(40,80), power = 1)
#heatcols <- heat_hcl(12,h=c(0,-100),l=c(75,40),c=c(40,80), power=1)
heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)


abundance_mat <- NA
plateaupoint <- 500


#RESEEDING FROM PREDICTED SCENARIO
ord <- order(reefdata$score) #lowest to highest
numsteps <- 50
num_remove <- 434
num_iter <- 20
final_reef <- (num_iter*num_remove)
reefdata$num <- seq(1,12292,by=1)
starterreefs <- reefdata$num[-ord[1:final_reef]]
reefdata$num <- NULL
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = (numsteps+1)) #50 time steps
abundance_mat[,1][starterreefs] <- 100 #time step 0: assign all the starterreefs to an abundance of 100

#calculate abundance after 1 time step
abundance_mat[,2] <- as.vector(t(abundance_mat[,1]) %*% connmat_reduced) #need to do this because it thinks of the right side as a matrix
#re-set abundance of all non-vacant reefs to 100
dummy <- abundance_mat[,2]
dummy[dummy > 0] <- 100 
#calculate abundance for the next time steps - note: the total abundance/sum of the subsequent columns should decrease every time step
for(i in 3:(numsteps+1)){
  abundance_mat[,i] <- as.vector(t(dummy) %*% connmat_reduced)
  dummy <- abundance_mat[,i]
  dummy[dummy > 0] <- 100
}


#Plot the number that were re-seeded over time
stop <- "no"
numreseeded_pr <- rep(NA,(numsteps+1))
numreseeded_pr[1] <- 0
for(i in 2:(numsteps+1)){
  numreseeded_pr[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[starterreefs,i][abundance_mat[starterreefs,i]>0])	
  if((numreseeded_pr[i] <= numreseeded_pr[i-1]) & stop == "no"){ #added 'or' condition to be consistent with the other scenarios
    plateaupoint <- i-1
    stop <- "yes"
  }
}
plateaupoint_pr <- plateaupoint

png(paste0("Realized_CU_Predicted_Reseed_",numsteps,"steps.png"))
plot(x=seq(1,(numsteps+1),by=1),y=numreseeded_pr,main="Realized Predicted Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "# of Reseeded Reefs")
dev.off()


id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)
#heatcols <- diverge_hcl(10)
#heatcols <- heat_hcl(12,h=c(0,-100), l = c(75,40), c = c(40,80), power = 1)
#heatcols <- heat_hcl(12,h=c(0,-100),l=c(75,40),c=c(40,80), power=1)
heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue

abundance_mat <- NA
plateaupoint <- 500


#RESEEDING FROM CATASTROPHE SCENARIO
ord <- order(reefdata$score, decreasing = TRUE) #highest to lowest
numsteps <- 50
num_remove <- 434
num_iter <- 20
final_reef <- (num_iter*num_remove)
reefdata$num <- seq(1,12292,by=1)
starterreefs <- reefdata$num[-ord[1:final_reef]]
reefdata$num <- NULL
abundance_mat <- matrix(data = 0, nrow = dim(connmat_reduced)[1], ncol = (numsteps+1)) #10 time steps
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


#Plot the number that were re-seeded over time
stop <- "no"
numreseeded_cr <- rep(NA,(numsteps+1))
numreseeded_cr[1] <- 0
for(i in 2:(numsteps+1)){
  numreseeded_cr[i] <- length(abundance_mat[,i][abundance_mat[,i]>0]) - length(abundance_mat[starterreefs,i][abundance_mat[starterreefs,i]>0])	
  if((numreseeded_cr[i] <= numreseeded_cr[i-1]) & stop == "no"){ #added in the 'or previous = bigger' condition bc of '0' issues
    plateaupoint <- i-1
    stop <- "yes"
  }
}
plateaupoint_cr <- plateaupoint


id.abund <- list()
for(i in 2:(numsteps+1)){
  id.abund[[i]] <- which((abundance_mat[,i] > 0) & (abundance_mat[,(i-1)] == 0))
  id.abund[[i]] <- id.abund[[i]][-starterreefs]
}
id.totalnotseeded <- which(abundance_mat[,(numsteps+1)] == 0)

heatcols <- matlab.like2((ultimateplateau+2))
#heatcols <- c('firebrick4', 'firebrick3', 'firebrick1','orangered1','darkorange1','gold1','yellow','aquamarine1','cyan1','cadetblue1', 'deepskyblue', 'dodgerblue1', 'royalblue1')#red -> yellow -> blue

abundance_mat <- NA


####

#load in numreseededrrlist from "reseeding/numreseeded.RData"
load(here("workingwith_computecanada_output","numreseeded.RData"))
#load in numreseededrr_mean from "reseeding/numreseededrr_all.RData"
load(here("workingwith_computecanada_output","numreseededrr_all.RData"))
random_iter <- 1000

#GOAL GRAPH
#NEW 7.27.2020 - adding thin lines around the random scenario one
#plotting all of the numreseeeded from the 4 different 'scenarios' on one graph - total seeded as a PERCENTTAGE OF TOTAL REEFS
#12292 - final_reef = number started with in the predicted and catastrophe situations
initialreefs = 12292 - final_reef
png(here("workingwith_computecanada_output",paste0("Realized_CU_ReseedFull_",numsteps,"steps_allscenarios_totalpercentages_allrandom_200kreps_no50r.png")))
#plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_50r+length(id.50reefs))/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "green", type = 'l', ylim = c(20,100), lwd = 3)
plot(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100),main="Realized Scenario Final Reefs Reseeded over Time",xlab = "Time Step", ylab = "%  of Reefs Reseeded", col = "#1d457f", type = 'l', ylim = c(20,100), lwd = 3) #blue
#lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_pr + initialreefs)/12292)*100), col = "blue",lwd = 3)
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseeded_cr + initialreefs)/12292)*100), col  = "#f9ad2a", lwd = 3) #darkgoldenrod1
for(i in 1:(random_iter*200)){
  lines(x = seq(1,(numsteps+1),by=1),y = (((numreseededrrlist[[i]] + initialreefs)/12292)*100),col = alpha("#e6b6c0", 0.2),lwd=1) ##c7c0d1
}
lines(x=seq(1,(numsteps+1),by=1),y=(((numreseededrr_mean + initialreefs)/12292)*100), col  = "#cc5c76",lty =3, lwd = 3) ##7f43d1
legend("bottomright", c("PCS Scenario", "Refuge-Loss Scenario","Mean(Random Scenario)", "Random Scenario"), col = c("#1d457f","#f9ad2a","#cc5c76","#e6b6c0"), lty = c(1,1,3,1), lwd = c(3,3,3,1))
#legend("bottomright", c("50 Reefs", "PCS Hyp", "Worst Hyp","Mean(Random Hyp)", "Random Hyp"), col = c("green","blue","orange","red","pink"), lty = c(1,1,1,3,1), lwd = c(3,3,3,3,1))
dev.off()

#need to figure out median # and IQR # of generations until plateau
load(here("workingwith_computecanada_output","plateaupoint_all.RData"))
summary(plateaupoint_rr)