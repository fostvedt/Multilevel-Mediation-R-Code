library(ggplot2)
#install.packages("lme4")
library(lme4)
library(doBy)
library(mice)

############### ############### 
# MI adjustments for final estimates
############### ############### 
finalMI <- function(var,beta,m){
mbeta <- apply(beta,1,mean)
B1 <- t(beta)- matrix(rep(mbeta,m),nrow=m,byrow=T)
Bm <- t(B1)%*%B1
Um <- var/m
Tm <- Um + ((m+1)/m)*Bm
y <- cbind(Estimate = mbeta,SE = sqrt(diag(Tm)), t.stat = mbeta/sqrt(diag(Tm)), p.value= round(2 * pnorm(abs(mbeta/sqrt(diag(Tm))), lower.tail = FALSE),4 ))
return(y)
}

data <- read.csv("/Users/lukefostvedt/Documents/swh/Teachers/data/VideoScores_Yr1.csv")
itbs <- read.csv("/Users/lukefostvedt/Documents/swh/data/itemized itbs csv files/itpforme11.csv")
itbs <- itbs[,1:45]
cct <- read.csv("/Users/lukefostvedt/Documents/swh/Teachers/data/ies_2011_2012_ctt.csv")




head(cct)
ind <- match(cct$Pre.Teacher.Code, data$New.Code)
cct$Rating <- data$Total.Score[ind]
cct$Improve <- cct$Post.TotalScore - cct$Pre.TotalScore
cct$TRT <- as.numeric(-1*(cct$TRT.CTL-2))
#cct$TRT <- factor(cct$TRT , levels=c(0,1), labels=c("Control","SWH"))
cct$Teacher <- factor(cct$Pre.Teacher.Code)
cct$School <- factor(cct$pre.SchoolName)


#selecting useful columns
indc <- match(cct$Student.ID,itbs$STATEID)
var <- c("SED","GAT","ELL","FRL","SEX")
dat <- cbind(cct,itbs[indc,var])

var <- names(dat)[c(1,26,27,25,23,6,10,18,28:32)]
df <- dat[,var]
df$Improve <- df$Post.TotalScore - df$Pre.TotalScore
apply(df,2,function(x) sum(is.na(x)))
df <- df[-which(df$Pre.TotalScore < 5 | df$Post.TotalScore < 5),]
indt <- match(unique(df$Teacher), df$Teacher)

level2imp <- summaryBy(data=df, Pre.TotalScore + Post.TotalScore + Rating + TRT + SEX+SED+FRL+GAT+ELL ~ Teacher, FUN=mean, na.rm=T )
names(level2imp) <- c("Teacher", "Pre", "Post", "Rating","TRT","SEX","SED","FRL","GAT","ELL")

level2imp$Rating[which(level2imp$Rating == "NaN")] <- NA
level2imp$School <- dat[ match( level2imp$Teacher,dat$Teacher)  ,"School"] 

md.pattern(level2imp)

udata <- list()
tdata <- list()

m <- 5
for (i in 1:(m+1)){

# imputing teacher ratings
ini <- mice(as.matrix(level2imp[,2:10]),maxit=5,m=1)
d1 <- complete(ini,1)
#d1$TRT <- factor(d1$TRT , levels=c(1,2), labels=c("SWH","Control"))
d1$Teacher <- level2imp$Teacher
d1$School <- level2imp$School
d1$Rating <- as.numeric(as.character(d1$Rating))
tdata[[i]] <- d1

#d1$Pre <- as.numeric(d1$Pre)
#d1$Post <- as.numeric(d1$Post)



# adding imputed ratings to file with missing demographics
df1 <- df
df1$Rating <- d1$Rating[match(df1$Teacher,d1$Teacher)]
ins <- mice(df1[,c(4,5,7:13)],maxit=5,m=1)
dat2 <- complete(ins,1)
dat2$Improve <- dat2$Post.TotalScore - dat2$Pre.TotalScore
udata[[i]] <- cbind(df[,1:3],dat2)

# reinserting the missing ratings with imputed demographics

dat <- cbind(df[,c(1:5)],dat2[,-c(1,2)])
level2imp <- summaryBy(data=dat, Pre.TotalScore + Post.TotalScore + Rating + TRT + SEX+SED+FRL+GAT+ELL ~ Teacher, FUN=mean, na.rm=T )
names(level2imp) <- c("Teacher", "Pre", "Post", "Rating","TRT","SEX","SED","FRL","GAT","ELL")
level2imp$Rating[which(level2imp$Rating == "NaN")] <- NA
level2imp$School <- d1$School

}




m1 <- lmer(Post.TotalScore ~  Rating +TRT + SED+GAT+ELL+FRL+SEX+ (1|Teacher)+(Rating|School),data=udata[[2]])
m2 <- lmer(Post.TotalScore ~  TRT +SED+GAT+ELL+FRL+SEX+  (1|Teacher)+(1|School), data=udata[[2]])
m3 <- lmer(Rating ~ TRT+ (1|School), data=tdata[[2]])

n1 <- length(fixef(m1))
n2 <- length(fixef(m2))
n3 <- length(fixef(m3))

var1 <- matrix(0,n1,n1)
beta1 <- matrix(0,n1,m)
var2 <- matrix(0,n2,n2)	
beta2 <- matrix(0,n2,m)
var3 <- matrix(0,n3,n3)	
beta3 <-matrix(0,n3,m)

varcom1 <- NULL
varcom2 <- NULL
varcom3 <- NULL


for(i in 2:6){
ds <- udata[[i]]
ds$Pre.TotalScore <- as.numeric(ds$Pre.TotalScore ) - mean(ds$Pre.TotalScore)
#ds$Rating <- as.numeric(ds$Rating)
dt <- tdata[[i]]
#dt$Rating <- as.numeric(dt$Rating)

m1 <- lmer(Improve~ TRT + Rating +SED+GAT+ELL+FRL+SEX + (1|Teacher)+(Rating|School),data=ds)
m2 <- lmer(Improve~ TRT +SED+GAT+ELL+FRL+SEX +  (1|Teacher)+(1|School), data=ds)
m3 <- lmer(Rating ~ TRT+ (1|School), data=dt)

var1 <- var1 + vcov(m1)
beta1[,i-1] <- fixef(m1)

var2 <- var2 + vcov(m2)
beta2[,i-1] <- fixef(m2)

var3 <- var3 + vcov(m3)
beta3[,i-1] <- fixef(m3)

varcom1 <- rbind(varcom1,c( attr(VarCorr(m1)$Teacher, "stddev")^2, attr(VarCorr(m1)$School, "stddev")^2, attr(VarCorr(m1), "sc")^2))
varcom2 <- rbind(varcom2,c(as.numeric(VarCorr(m2)),attr(VarCorr(m2), "sc")^2))
varcom3 <- rbind(varcom3,c(as.numeric(VarCorr(m3)),attr(VarCorr(m3), "sc")^2))

}


vc1 <- apply(varcom1,2,mean) 
vc2 <- apply(varcom2,2,mean) 
vc3 <- apply(varcom3,2,mean) 

l1 <-finalMI(var1,beta1,m)
l2 <-finalMI(var2,beta2,m)
l3 <-finalMI(var3,beta3,m)

row.names(l1) <- names(fixef(m1))
row.names(l2) <- names(fixef(m2))
row.names(l3) <- names(fixef(m3))

res <- list(l1,vc1,l2,vc2,l3,vc3)

write.matrix(round(res[[1]],3),sep=" &")
write.matrix(round(res[[2]],3),sep=" &")
write.matrix(round(res[[3]],3),sep=" &")
write.matrix(round(res[[4]],3),sep=" &")
write.matrix(round(res[[5]],3),sep=" &")
write.matrix(round(res[[6]],3),sep=" &")







############### ############### 
# MI adjustments for final estimates
############### ############### 
finalMI <- function(var,beta,m){
mbeta <- apply(beta,1,mean)
B1 <- t(beta)- matrix(rep(mbeta,m),nrow=m,byrow=T)
Bm <- t(B1)%*%B1
Um <- var/m
Tm <- Um + ((m+1)/m)*Bm
y <- cbind(Estimate = mbeta,SE = sqrt(diag(Tm)), t.stat = mbeta/sqrt(diag(Tm)), p.value= round(2 * pnorm(abs(mbeta/sqrt(diag(Tm))), lower.tail = FALSE),4 ))
return(y)
}


finalMI(var1,beta1,5)
finalMI(var2,beta2,5)
finalMI(var3,beta3,5)






cct <- read.csv("/Users/lukefostvedt/Documents/swh/Teachers/data/ies_2011_2012_ctt.csv")
# model without imputing demographics
head(cct)
ind <- match(cct$Pre.Teacher.Code, data$New.Code)
cct$Rating <- data$Total.Score[ind]
cct$Improve <- cct$Post.TotalScore - cct$Pre.TotalScore
cct$TRT <- as.numeric(-1*(cct$TRT.CTL-2))
#cct$TRT <- factor(cct$TRT , levels=c(0,1), labels=c("Control","SWH"))
cct$Teacher <- factor(cct$Pre.Teacher.Code)
cct$School <- factor(cct$pre.SchoolName)


#selecting useful columns
indc <- match(cct$Student.ID,itbs$STATEID)
var <- c("SED","GAT","ELL","FRL","BLK","HSP","SEX")
dat <- cbind(cct,itbs[indc,var])

var <- names(dat)[c(1,26,27,25,23,6,10,18,24,28:34)]
df <- dat[,var]
apply(df,2,function(x) sum(is.na(x)))
df <- df[-which(df$Pre.TotalScore < 5 | df$Post.TotalScore < 5),]
indt <- match(unique(df$Teacher), df$Teacher)


level2imp <- summaryBy(data=df, Pre.TotalScore + Post.TotalScore + Rating + TRT ~ Teacher, FUN=mean, na.rm=T )
names(level2imp) <- c("Teacher", "Pre", "Post", "Rating","TRT")
level2imp$Rating[which(level2imp$Rating == "NaN")] <- NA
level2imp$School <- df[ match( level2imp$Teacher,df$Teacher)  ,"School"] 
#level2imp$TRT  <- factor(level2imp$TRT , levels=c(0,1), labels=c("Control","SWH"))





miiter <- function(seed,iter){
set.seed(seed)
m <- iter 
ins <- mice(level2imp[,-1],m=iter,maxit=iter)
df1 <- df
td <- list()
sd <- list()
for (i in 1:m ){
	tdat <- cbind(Teacher=level2imp$Teacher,complete(ins,i))
	td[[i]] <- tdat
	df1$Rating <- tdat$Rating[match(df1$Teacher,tdat$Teacher)]
	sd[[i]] <- df1
}


m1 <- lmer(Post.TotalScore~ Pre.TotalScore+ Rating +TRT + SED+GAT+ELL+FRL+ (1|Teacher)+(Rating|School),data=sd[[1]])
m2 <- lmer(Post.TotalScore~Pre.TotalScore+ TRT +SED+GAT+ELL+FRL+  (1|Teacher)+(1|School), data=sd[[1]])
m3 <- lmer(Rating ~ TRT+ (1|School), data=td[[1]])

n1 <- length(fixef(m1))
n2 <- length(fixef(m2))
n3 <- length(fixef(m3))

var1 <- matrix(0,n1,n1)
beta1 <- matrix(0,n1,m)
var2 <- matrix(0,n2,n2)	
beta2 <- matrix(0,n2,m)
var3 <- matrix(0,n3,n3)	
beta3 <-matrix(0,n3,m)

varcom1 <- NULL
varcom2 <- NULL
varcom3 <- NULL

for(i in 1:m){
ds <- sd[[i]]
ds$Pre.TotalScore <- as.numeric(ds$Pre.TotalScore ) - mean(as.numeric(ds$Pre.TotalScore ))
ds$Post.TotalScore <- as.numeric(ds$Post.TotalScore ) - mean(as.numeric(ds$Post.TotalScore ))
ds$Rating <- as.numeric(ds$Rating)
dt <- td[[i]]
dt$Rating <- as.numeric(dt$Rating)


#m1 <- lmer(Improve~ Rating +TRT+SED+GAT+ELL+FRL+ (1|Teacher)+(Rating |School),data=ds)
#m2 <- lmer(Improve~   TRT + SED+GAT+ELL+FRL+ (1|Teacher)+(1|School), data=ds)
#m3 <- lmer(Rating ~ TRT+ (1|School), data=dt)

m1 <- lmer(Improve~Pre.TotalScore +  TRT + Rating +SED+GAT+ELL+FRL+ (1|Teacher)+(Rating|School),data=ds)
m2 <- lmer(Improve ~Pre.TotalScore +  TRT + SED+GAT+ELL+FRL+ (1|Teacher)+(1|School), data=ds)
m3 <- lmer(Rating ~ TRT+ (1|School), data=dt)

#m1 <- lmer(Improve~Pre.TotalScore +  TRT + Rating  + (1|Teacher)+(1|School),data=ds)
#m2 <- lmer(Improve~Pre.TotalScore +  TRT +  (1|Teacher)+(1|School), data=ds)
#m3 <- lmer(Rating ~ TRT+ (1|School), data=dt)

var1 <- var1 + vcov(m1)
beta1[,i] <- fixef(m1)
var2 <- var2 + vcov(m2)
beta2[,i] <- fixef(m2)
var3 <- var3 + vcov(m3)
beta3[,i] <- fixef(m3)

varcom1 <- rbind(varcom1,c( attr(VarCorr(m1)$Teacher, "stddev")^2, attr(VarCorr(m1)$School, "stddev")^2, attr(VarCorr(m1), "sc")^2))
varcom2 <- rbind(varcom2,c(as.numeric(VarCorr(m2)),attr(VarCorr(m2), "sc")^2))
varcom3 <- rbind(varcom3,c(as.numeric(VarCorr(m3)),attr(VarCorr(m3), "sc")^2))
}


vc1 <- apply(varcom1,2,mean) 
vc2 <- apply(varcom2,2,mean) 
vc3 <- apply(varcom3,2,mean) 

l1 <-finalMI(var1,beta1,m)
l2 <-finalMI(var2,beta2,m)
l3 <-finalMI(var3,beta3,m)

row.names(l1) <- names(fixef(m1))
row.names(l2) <- names(fixef(m2))
row.names(l3) <- names(fixef(m3))


return(list(l1,vc1,l2,vc2,l3,vc3))
}

miiter(235,5) -> res



write.matrix(round(res[[1]],3),sep=" &")
write.matrix(round(res[[2]],3),sep=" &")
write.matrix(round(res[[3]],3),sep=" &")
write.matrix(round(res[[4]],3),sep=" &")
write.matrix(round(res[[5]],3),sep=" &")
write.matrix(round(res[[6]],3),sep=" &")

