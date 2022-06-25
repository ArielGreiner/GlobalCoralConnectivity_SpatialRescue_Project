##NOTE: this code will not run because was unable to upload the Wood et al. 2014 monthly-connectivity-matrices to github as the folder was too large.

require(R.matlab)

#create empty matrices
sum_mat <- matrix(0,nrow=12397,ncol=12397)
sum_mat_avg <- matrix(0,nrow = 12397, ncol = 12397)

#2003 (because only 2 months)
for(month in 11:12){
  #REPLACE LINE BELOW WITH THE DIRECTORY THAT YOU PUT THE MONTHLY-CONNECTIVITY-MATRICES FOLDER IN
	name.file <- paste0("~/.../monthly-connectivity-matrices/SetMat_2003_",month,".mat") 
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
	for(month in 1:12){ 
	  #REPLACE LINE BELOW WITH THE DIRECTORY THAT YOU PUT THE MONTHLY-CONNECTIVITY-MATRICES FOLDER IN
	name.file <- paste0("~/.../monthly-connectivity-matrices/SetMat_",year,"_",month,".mat")
a <- readMat(name.file)
sum_mat <- sum_mat + a$settle
sum_mat_avg <- sum_mat_avg + (a$settle)/100
print(paste("year = ",year,"month = ", month))
flush.console() #prints as executing, takes less time
}}
save(sum_mat, file = "sum_mat.RData")
save(sum_mat_avg, file = "sum_mat_avg.RData")

for(month in 1:7){
  #REPLACE LINE BELOW WITH THE DIRECTORY THAT YOU PUT THE MONTHLY-CONNECTIVITY-MATRICES FOLDER IN
	name.file <- paste0("~/.../monthly-connectivity-matrices/SetMat_2011_",month,".mat")
a <- readMat(name.file)
sum_mat <- sum_mat + a$settle
sum_mat_avg <- sum_mat_avg + (a$settle)/100
print(paste("year = 2011, month = ", month))
flush.console() #prints as executing, takes less time
}
#save(sum_mat, file = "sum_mat.RData")
#save(sum_mat_avg, file = "sum_mat_avg.RData")

Conn_Mat_Sum <- sum_mat/(93*100)
Conn_Mat_Avg <- sum_mat_avg/93
#save(Conn_Mat_Sum, file = "Conn_Mat_Sum.RData")
#save(Conn_Mat_Avg, file = "Conn_Mat_Avg.RData")

#create a sparse matrix - much smaller size, only works if matrix has lots of 0s because only changes the 0s...can still deal with sparse matrices as normal (add, subtract, etc) <- according to Marco
#save some sparse ones as well
library(Matrix)

Conn_Mat_Sum_sparse <- Matrix(Conn_Mat_Sum, sparse = T)
Conn_Mat_Avg_sparse <- Matrix(Conn_Mat_Avg, sparse = T)
save(Conn_Mat_Sum_sparse, file = "~/Github/GlobalCoralConnectivity_SpatialRescue_Project/OriginalDataConversion/ConvertingWoodetal2014Data/Conn_Mat_Sum_sparse.RData")
#save(Conn_Mat_Avg_sparse, file = "Conn_Mat_Avg_sparse.RData")

#NOTE: the sum and average connectivity matrices are identical.
