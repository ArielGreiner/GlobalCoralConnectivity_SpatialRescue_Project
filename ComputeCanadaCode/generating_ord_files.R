#compute canada version of the Random Scenario Removal more reps - do 100 chunks of 2000 reps at a time?
#generate 200,000 'ords' first - split into 100 2000-column sets
#each set should also have a column with a number (1-100) so that the outputs don't save on top of one another
#then after, join all 100 lists together for each output
setwd('~/Documents/computecanada_randomremovalmorereps_june2020/ord_mats')
set.seed(2)

load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/cells12292.RData') #this includes 'id.noReef' and 'b' (shapefile)
reefdata <- b@data

for(i in 1:100){
	ord_mat <- matrix(NA, nrow = 12292, ncol = 2001)
	ord_mat[,1] <- i
	for(j in 1:2000){
	ord_mat[,(j+1)] <- sample(seq(1,length(reefdata$score),by=1))	
	}
print(paste("i =",i))
write.table(ord_mat, file=paste0("ord_mat_",i,".csv"),sep=",",row.names=FALSE,col.names=FALSE)
}

#generated 1-75, stopped

#just need to generate a bunch of setseed values, not 200k orderings
#this worked but it kept running out of memory
setwd('~/Documents/computecanada_randomremovalmorereps_june2020')
#i think i just need to generate 1-200
#1000 takes 8hrs on my computer, so just run 200 sets in parallel
iterator <- seq(1,200,1)
write.table(iterator[1:40], file="iterator1_randomrep.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(iterator[41:80], file="iterator2_randomrep.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(iterator[81:120], file="iterator3_randomrep.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(iterator[121:160], file="iterator4_randomrep.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(iterator[161:200], file="iterator5_randomrep.csv",sep=",",row.names=FALSE,col.names=FALSE)

#sets of 3 worked but some still ran out of memory
iterator <- seq(1,200,1)
for(i in 1:66){
write.table(iterator[(((i-1)*3)+1):(i*3)], file=paste0("iterator_randomrep",i,".csv"),sep=",",row.names=FALSE,col.names=FALSE)
}
write.table(iterator[199:200], file=paste0("iterator_randomrep67.csv"),sep=",",row.names=FALSE,col.names=FALSE)

#need to add two more iterators: 72,76,169,174
write.table(c(72,76), file=paste0("iterator_randomrep68.csv"),sep=",",row.names=FALSE,col.names=FALSE)
write.table(c(169,174), file=paste0("iterator_randomrep69.csv"),sep=",",row.names=FALSE,col.names=FALSE)


#new iterators with only one thing in them for the ones that failed early
setwd('~/Documents/computecanada_randomremovalmorereps_june2020/iterators')
iterator_redux <- c(16,17,18,37,38,39,46,47,48,67,68,69,145,146,147,151,152,153,160,161,162,187,188,189,196,197,198)
for(i in 1:27){
write.table(iterator_redux[i], file=paste0("iterator_randomrep",i+69,".csv"),sep=",",row.names=FALSE,col.names=FALSE)
}

