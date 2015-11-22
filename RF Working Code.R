# Load necessary packages:
library(randomForest)

# Read in the data
# Code to Read Data and Run RF OoB imp. measures and ranks
df.feat <- read.csv(file="/Users/lucasmentch/Desktop/Mahya DNA Methylation/Newest/TCGA_FEAT.csv",header=T)
nfeat <- dim(df.feat)[2]
df.resp <- read.csv(file="/Users/lucasmentch/Desktop/Mahya DNA Methylation/Newest/TCGA_RESP.csv",header=F)

# Separate Responses:
# Respose 1:  GSTP1; body
df.resp.1 <- df.resp[1,]
df.full.1 <- data.frame("y"=as.numeric(df.resp.1),df.feat)

# Respose 2:  RASSF1; promoter
df.resp.2 <- df.resp[2,]
df.full.2 <- data.frame("y"=as.numeric(df.resp.2),df.feat)

# Respose 3:  PITX2; body
df.resp.3 <- df.resp[3,]
df.full.3 <- data.frame("y"=as.numeric(df.resp.3),df.feat)

# Fix/Impute NAs
df.full.rough.1 <- na.roughfix(df.full.1)
df.full.rough.2 <- na.roughfix(df.full.2)
df.full.rough.3 <- na.roughfix(df.full.3)

# Create train and test sets for each response:
ntest <- 24
ind <- sample(1:424,ntest,replace=FALSE)
df1.test <- df.full.rough.1[ind,]
df1.train <- df.full.rough.1[-ind,]
df2.test <- df.full.rough.2[ind,]
df2.train <- df.full.rough.2[-ind,]
df3.test <- df.full.rough.3[ind,]
df3.train <- df.full.rough.3[-ind,]

# Create "bad" data
randfeat1 <- data.frame(matrix(rnorm(dim(df.feat)[1]*dim(df.feat)[2]),nrow=dim(df.feat)[1],ncol=dim(df.feat)[2]))
randfeat2 <- data.frame(matrix(rnorm(dim(df.feat)[1]*dim(df.feat)[2]),nrow=dim(df.feat)[1],ncol=dim(df.feat)[2]))
randfeat3 <- data.frame(matrix(rnorm(dim(df.feat)[1]*dim(df.feat)[2]),nrow=dim(df.feat)[1],ncol=dim(df.feat)[2]))
rand1.test <- data.frame("y"=as.numeric(df.resp.1)[ind],randfeat1[ind,])
rand1.train <- data.frame("y"=as.numeric(df.resp.1)[-ind],randfeat1[-ind,])
rand2.test <- data.frame("y"=as.numeric(df.resp.2)[ind],randfeat2[ind,])
rand2.train <- data.frame("y"=as.numeric(df.resp.2)[-ind],randfeat2[-ind,])
rand3.test <- data.frame("y"=as.numeric(df.resp.3)[ind],randfeat3[ind,])
rand3.train <- data.frame("y"=as.numeric(df.resp.3)[-ind],randfeat3[-ind,])


# Run the random forests
nsim <- 100
# Storing the predictions from each RF
pred1 <- matrix(NA,nrow=nsim,ncol=ntest)
pred2 <- matrix(NA,nrow=nsim,ncol=ntest)
pred3 <- matrix(NA,nrow=nsim,ncol=ntest)

# Storing the predictions from RF with random data
rpred1 <- matrix(NA,nrow=nsim,ncol=ntest)
rpred2 <- matrix(NA,nrow=nsim,ncol=ntest)
rpred3 <- matrix(NA,nrow=nsim,ncol=ntest) 

for (i in 1:nsim) {
	
	# Run RF on good data
	rf1 <- randomForest(y~.,data=df1.train,importance=FALSE,xtest=df1.test[,-1],ytest=df1.test[,1])
	rf2 <- randomForest(y~.,data=df2.train,importance=FALSE,xtest=df2.test[,-1],ytest=df2.test[,1])
	rf3 <- randomForest(y~.,data=df3.train,importance=FALSE,xtest=df3.test[,-1],ytest=df3.test[,1])
	
	# Run RF on random data
	rf1.rand <- randomForest(y~.,data=rand1.train,importance=FALSE,xtest=rand1.test[,-1],ytest=rand1.test[,1])
	rf2.rand <- randomForest(y~.,data=rand2.train,importance=FALSE,xtest=rand2.test[,-1],ytest=rand2.test[,1])
	rf3.rand <- randomForest(y~.,data=rand3.train,importance=FALSE,xtest=rand3.test[,-1],ytest=rand3.test[,1])
	
	# Save predictions from good data
	pred1[i,] <- as.numeric(rf1$test$predicted)
	pred2[i,] <- as.numeric(rf2$test$predicted)
	pred3[i,] <- as.numeric(rf3$test$predicted)
	
	# Save predictions from random data
	rpred1[i,] <- as.numeric(rf1.rand$test$predicted)
	rpred2[i,] <- as.numeric(rf2.rand$test$predicted)
	rpred3[i,] <- as.numeric(rf3.rand$test$predicted)
	
	# Verbose:
	cat("Iteration ",i," of ",nsim,"\n")
	
	
}

# Calculate average MSE on real data:
mse.1 <- apply((pred1 - matrix(df1.test[,1],byrow=TRUE,nrow=nsim,ncol=ntest))^2,1,mean)
mse.2 <- apply((pred2 - matrix(df2.test[,1],byrow=TRUE,nrow=nsim,ncol=ntest))^2,1,mean)
mse.3 <- apply((pred3 - matrix(df3.test[,1],byrow=TRUE,nrow=nsim,ncol=ntest))^2,1,mean)

# Calculate average MSE on random data:
mse.rand.1 <- apply((rpred1 - matrix(rand1.test[,1],byrow=TRUE,nrow=nsim,ncol=ntest))^2,1,mean)
mse.rand.2 <- apply((rpred2 - matrix(rand2.test[,1],byrow=TRUE,nrow=nsim,ncol=ntest))^2,1,mean)
mse.rand.3 <- apply((rpred3 - matrix(rand3.test[,1],byrow=TRUE,nrow=nsim,ncol=ntest))^2,1,mean)

# Plot the MSEs:
plot(mse.1,type='l',col='blue',ylim=c(0.0045,0.012),main="Random Forest Test Set MSE \n GSTP1",xlab="Random Forest Index",ylab="MSE",lwd=2)
lines(mse.rand.1,lwd=2)
legend("topright",c("Simulated Features","Real Features"),lty=c(1,1),col=c("black","blue"),lwd=c(2,2))

plot(mse.2,type='l',col='blue',ylim=c(0.00115,0.00235),main="Random Forest Test Set MSE \n RASSF1",xlab="Random Forest Index",ylab="MSE",lwd=2)
lines(mse.rand.2,lwd=2)
legend("topright",c("Simulated Features","Real Features"),lty=c(1,1),col=c("black","blue"),lwd=c(2,2))

plot(mse.3,type='l',col='blue',ylim=c(0.0035,0.0069),main="Random Forest Test Set MSE \n PITX2",xlab="Random Forest Index",ylab="MSE",lwd=2)
lines(mse.rand.3,lwd=2)
legend("topright",c("Simulated Features","Real Features"),lty=c(1,1),col=c("black","blue"),lwd=c(2,2))

# Plot MSE-good vs MSE-rand for each of the three pairs:
mse.mat <- matrix(c(mean(mse.1),mean(mse.rand.1),mean(mse.2),mean(mse.rand.2),mean(mse.3),mean(mse.rand.3)),byrow=FALSE,ncol=3)
colnames(mse.mat) <- c("GSTP1","RASSF1","PITX2")
barplot(mse.mat,beside=TRUE,col=c("blue","black"),main="Average MSE Across Forests")
legend("topright",c("Original Data","Simulated Data"),lwd=c(5,5),col=c("blue","black"))
box()

# Now generate data from a linear model:
beta <- matrix(runif(dim(df.feat)[2]),ncol=1)
X <- matrix(rnorm(dim(df.feat)[1]*dim(df.feat)[2]),ncol=dim(df.feat)[2])
eps <- rnorm(dim(df.feat)[1],sd=0.05)
y.lm <- X %*% beta + eps
df.lm <- data.frame("y"=y.lm,X)
X.rand <- matrix(rnorm(dim(df.feat)[1]*dim(df.feat)[2]),ncol=dim(df.feat)[2])
df.lm.rand <- data.frame("y"=y.lm,X.rand)

# Create train and test sets for each response:
ntest.lm <- 24
ind.lm <- sample(1:424,ntest.lm,replace=FALSE)
df.lm.test <- df.lm[ind.lm,]
df.lm.train <- df.lm[-ind.lm,]
df.lm.rand.test <- df.lm.rand[ind.lm,]
df.lm.rand.train <- df.lm.rand[-ind.lm,]

# Run the LM random forests
nsim.lm <- 100
# Storing the predictions from each RF
pred.lm <- matrix(NA,nrow=nsim.lm,ncol=ntest.lm)
rpred.lm <- matrix(NA,nrow=nsim.lm,ncol=ntest.lm)

for (i in 1:nsim.lm) {
	
	# Run RF on good data
	rf.lm <- randomForest(y~.,data=df.lm.train,importance=FALSE,xtest=df.lm.test[,-1],ytest=df.lm.test[,1])
	rf.lm.rand <- randomForest(y~.,data=df.lm.rand.train,importance=FALSE,xtest=df.lm.rand.test[,-1],ytest=df.lm.rand.test[,1])
	
	# Save predictions from good data
	pred.lm[i,] <- as.numeric(rf.lm$test$predicted)
	rpred.lm[i,] <- as.numeric(rf.lm.rand$test$predicted)
	
	# Verbose:
	cat("Iteration ",i," of ",nsim,"\n")
	
}

# Calculate the MSEs:
mse.lm <- apply((pred.lm - matrix(df.lm.test[,1],byrow=TRUE,nrow=nsim.lm,ncol=ntest.lm))^2,1,mean)
mse.lm.rand <- apply((rpred.lm - matrix(df.lm.rand.test[,1],byrow=TRUE,nrow=nsim.lm,ncol=ntest.lm))^2,1,mean)

# Plot the MSEs:
plot(mse.lm,type='l',col='blue',ylim=c(60,100),main="Random Forest Test Set MSE \n Underlying Linear Model",xlab="Random Forest Index",ylab="MSE",lwd=2)
lines(mse.lm.rand,lwd=2)
legend("topright",c("Simulated Features","Real Features"),lty=c(1,1),col=c("black","blue"),lwd=c(2,2))

# How well would a lm do?
# lm1.true <- lm(y~.,data=df.lm.train)
# mse.lm.true <- mean( (predict(lm1.true,df.lm.test[,-1]) - df.lm.test[,1])^2 )
# lm1.rand <- lm(y~.,data=df.lm.rand.train)
# mse.lm.true.rand <- mean( (predict(lm1.rand,df.lm.rand.test[,-1]) - df.lm.rand.test[,1])^2 )

# Plot the ratio of the MSEs:
mse.ratio.1 <- mse.rand.1/mse.1
mse.ratio.2 <- mse.rand.2/mse.2
mse.ratio.3 <- mse.rand.3/mse.3
mse.ratio.lm <- mse.lm.rand/mse.lm
plot(mse.ratio.lm,type='l',col='black',ylim=c(0.95,2.75),main="MSE Ratios across Underlying Models",xlab="Random Forest Index",ylab="MSE Ratio",lwd=2)
lines(mse.ratio.1,col='blue',lwd=2)
lines(mse.ratio.2,col='red',lwd=2)
lines(mse.ratio.3,col='green',lwd=2)
abline(h=1,lty=2)
legend("topright",c("Linear Model","GSTP1","RASSF1","PITX2"),lty=c(1,1,1,1),col=c("black","blue","red","green"),lwd=c(2,2,2,2))