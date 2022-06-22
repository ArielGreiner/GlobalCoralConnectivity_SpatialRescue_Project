require(R.matlab)

#create empty matrices
sum_mat <- matrix(0,nrow=12397,ncol=12397)
sum_mat_avg <- matrix(0,nrow = 12397, ncol = 12397)

#2003 (because only 2 months)
for(month in 11:12){
	name.file <- paste0("~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Wood-et-al_2014_GEB_data/monthly-connectivity-matrices/SetMat_2003_",month,".mat") #paste0 removes spaces
a <- readMat(name.file)
sum_mat <- sum_mat + a$settle #$settle isolates the matrix part
sum_mat_avg <- sum_mat_avg + (a$settle)/100
print(paste("year = 2003, month = ", month))
flush.console() #prints as executing, takes less time
}
save(sum_mat, file = "sum_mat.RData")
save(sum_mat_avg, file = "sum_mat_avg.RData")

#2004 - 2010 (all 12 months)
for(year in 2004:2010){
	for(month in 1:12){ #used to be 11:12
	name.file <- paste0("~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Wood-et-al_2014_GEB_data/monthly-connectivity-matrices/SetMat_",year,"_",month,".mat")
a <- readMat(name.file)
sum_mat <- sum_mat + a$settle
sum_mat_avg <- sum_mat_avg + (a$settle)/100
print(paste("year = ",year,"month = ", month))
flush.console() #prints as executing, takes less time
}}
save(sum_mat, file = "sum_mat.RData")
save(sum_mat_avg, file = "sum_mat_avg.RData")

for(month in 1:7){
	name.file <- paste0("~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Wood-et-al_2014_GEB_data/monthly-connectivity-matrices/SetMat_2011_",month,".mat")
a <- readMat(name.file)
sum_mat <- sum_mat + a$settle
sum_mat_avg <- sum_mat_avg + (a$settle)/100
print(paste("year = 2011, month = ", month))
flush.console() #prints as executing, takes less time
}
save(sum_mat, file = "sum_mat.RData")
save(sum_mat_avg, file = "sum_mat_avg.RData")
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/sum_mat_avg.RData')
load('~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/sum_mat.RData')
#93 = the number of matrices being summed both times
#Conn_Mat_Sum = add all of the matrices together, then divide every element of the summed matrix by 93*100
#Conn_Mat_Avg = divide every element of every matrix by 100 to turn it into a connectivity matrix (bc otherwise it is just a number of particles that made it b/w two grid cells and not a proportion) and then take the average of all 93 matrices
#might turn out to be the same thing, if not -> use the Sum one because that one should inflate the rare occurrences (as opposed to deflate) which is more accurate since only 100 particles sent out per patch
Conn_Mat_Sum <- sum_mat/(93*100)
Conn_Mat_Avg <- sum_mat_avg/93
save(Conn_Mat_Sum, file = "Conn_Mat_Sum.RData")
save(Conn_Mat_Avg, file = "Conn_Mat_Avg.RData")

#create a sparse matrix - much smaller size, only works if matrix has lots of 0s because only changes the 0s...can still deal with sparse matrices as normal (add, subtract, etc) <- according to Marco
#save some sparse ones as well
library(Matrix)
Conn_Mat_Sum_sparse <- Matrix(Conn_Mat_Sum, sparse = T)
Conn_Mat_Avg_sparse <- Matrix(Conn_Mat_Avg, sparse = T)
save(Conn_Mat_Sum_sparse, file = "Conn_Mat_Sum_sparse.RData")
save(Conn_Mat_Avg_sparse, file = "Conn_Mat_Avg_sparse.RData")
