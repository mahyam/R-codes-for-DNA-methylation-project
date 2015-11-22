# Load libraries:
library(MASS)
library(rpart)
library(randomForest)

# Load the data:
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

Bagged.Inf <- function(train=df.demo,test=test.demo,testvars=c(5,6),verbose=TRUE,k=75,nx1=50,nmc=250,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0) {

# Defining rpart Control Parameters
control.sim <- rpart.control(minsplit=minsplit,maxcompete=maxcompete,maxsurrogate=maxsurrogate,usesurrogate=usesurrogate)

# Defining the size of the training set and ensemble
n <- dim(train)[1]
m <- nx1*nmc

# Defining the reduced data:
train.red <- train[,-testvars]
test.red <- test[,-testvars]

# Build the trees and estimate the parameters:
pred.all <- matrix(0,nrow=1,ncol=dim(test)[1])
diff.all <- matrix(0,nrow=1,ncol=dim(test)[1])
cond.exp.full <- matrix(0,nrow=nx1,ncol=dim(test)[1])
cond.exp.diff <- matrix(0,nrow=nx1,ncol=dim(test)[1])
for (i in 1:nx1) {
	ind.x1 <- sample(1:dim(train)[1],size=1,replace=FALSE)
	pred.full <- matrix(0,nrow=nmc,ncol=dim(test)[1])
	pred.red <- matrix(0,nrow=nmc,ncol=dim(test)[1])
	pred.diff <- matrix(0,nrow=nmc,ncol=dim(test)[1])
	for (j in 1:nmc) {
		ind <- c(ind.x1,sample((1:dim(train)[1])[-ind.x1],k-1,replace=FALSE))
		ss.full <- train[ind,]	
		ss.red <- train.red[ind,]
		tree.full <- rpart(y~.,data=ss.full,control=control.sim)
		tree.red <- rpart(y~.,data=ss.red,control=control.sim)
		pred.full[j,] <- predict(tree.full,test)
		pred.red[j,] <- predict(tree.red,test.red)
		pred.diff[j,] <- pred.full[j,] - pred.red[j,]
		if (verbose) cat("nx1:  ",i,"          nmc:  ",j,"\n")
	}
	pred.all <- rbind(pred.all,pred.full)
	diff.all <- rbind(diff.all,pred.diff)
	cond.exp.full[i,] <- apply(pred.full,2,mean)
	cond.exp.diff[i,] <- apply(pred.diff,2,mean)
}
pred.all <- pred.all[-1,]
diff.all <- diff.all[-1,]

mean.full <- apply(pred.all,2,mean)
mean.diff <- apply(diff.all,2,mean)

zeta1.full <- apply(cond.exp.full,2,var)
zeta1.diff <- cov(cond.exp.diff)

zetak.full <- apply(pred.all,2,var)
zetak.diff <- cov(diff.all)

sd.full <- sqrt((m/n)*((k^2)/m)*zeta1.full + (1/m)*zetak.full)
lbounds <- qnorm(0.025,mean=mean.full,sd=sd.full)
ubounds <- qnorm(0.975,mean=mean.full,sd=sd.full)

tstat <- t(mean.diff) %*% ginv((m/n)*((k^2)/m)*zeta1.diff + (1/m)*zetak.diff) %*% mean.diff
#1-pchisq(t.stat,df=length(mean.diff))
pval <- 1-pchisq(tstat,df=dim(test)[1])
return(list("lbounds"=lbounds,"ubounds"=ubounds,"pred"=mean.full,"tstat"=tstat,"pval"=pval))
}

# group.names[10:13] <- "Methionine Cycle"
# group.names[152:201] <- "SGOC Network" # 152 - 201 
# ind <- sample(1:424,ntest,replace=FALSE)
# df1.test <- df.full.rough.1[ind,]
# df1.train <- df.full.rough.1[-ind,]
# df2.test <- df.full.rough.2[ind,]
# df2.train <- df.full.rough.2[-ind,]
# df3.test <- df.full.rough.3[ind,]
# df3.train <- df.full.rough.3[-ind,]

# pval MC1:  0.5941
# pval MC+SGOC1:  0.0274
#
#
#
MC1 <- Bagged.Inf(train=df1.train,test=df1.test,testvars=11:14,verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)
save(MC1,file="MC1.rda")
SGOCMC1 <- Bagged.Inf(train=df1.train,test=df1.test,testvars=c(11:14,153:202),verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)
save(SGOCMC1,file="SGOCMC1.rda")

MC2 <- Bagged.Inf(train=df2.train,test=df2.test,testvars=11:14,verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)
save(MC2,file="MC2.rda")
SGOCMC2 <- Bagged.Inf(train=df2.train,test=df2.test,testvars=c(11:14,153:202),verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)
save(SGOCMC2,file="SGOCMC2.rda")

MC3 <- Bagged.Inf(train=df3.train,test=df3.test,testvars=11:14,verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)
save(MC3,file="MC3.rda")
SGOCMC3 <- Bagged.Inf(train=df3.train,test=df3.test,testvars=c(11:14,153:202),verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)
save(SGOCMC3,file="SGOCMC3.rda")

LM1large <- Bagged.Inf(train=df.lm.train,test=df.lm.test,testvars=c(11:14,153:202),verbose=TRUE,k=50,nx1=50,nmc=500,minsplit=3,maxcompete=0,maxsurrogate=0,usesurrogate=0)