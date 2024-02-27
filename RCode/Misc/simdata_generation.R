# Simulate data from SimFMRI123() function
set.seed(1997)
simData = SimFMRI123(noisyICA=F, nTR=50,var.inactive = 0, snr=0.2) 
xmat = simData$X
smat = simData$S
write.table(xmat,file = "Xmat.csv",row.names = F,col.names = F,sep = ",",quote = F)
write.table(smat,file = "Smat.csv",row.names = F,col.names = F,sep = ",",quote = F)