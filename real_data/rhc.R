library(dsmatch1)


bootamw<-function(A,X,Y,M,bn=100){

  sample_size=as.numeric(dim(X)[1])
  tsize=sum(A)

  #n=build model
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)

  zp = model.matrix(Y ~. ,data = Z)
  Trt=rep(0,sample_size)
  Trt[1:tsize]=rep(1,tsize)
  bs=rep(0,bn)

  i=1

  while(i <=bn){
    print(i)

    sampledata1=sample(loc.1,sum(A),replace = TRUE)
    sampledata0=sample(loc.0,(sample_size-sum(A)),replace = TRUE)
    ## make new
    ya=rep(0,sample_size)
    ya[1:tsize]=Y[sampledata1]
    ya[(tsize+1):sample_size]<-Y[sampledata0]
    ####z
    Z = data.frame(X,ya)

    Z1=Z[sampledata1,]
    pg1m = try(model.matrix(ya ~.,data = Z1)[,-1],silent=TRUE)
    Z0=Z[sampledata0,]
    pg0m = try(model.matrix(ya ~. ,data = Z0)[,-1],silent=TRUE)
    if ('try-error' %in% class(pg1m)){
      next}
    if ('try-error' %in% class(pg0m)){
      next}

    Xmat=X
    Xmat[1:tsize,]=X[sampledata1,]
    Xmat[(tsize+1):sample_size,]<-X[sampledata0,]
    Zmat=data.frame(Xmat,ya)
    zp=model.matrix(ya ~. ,data = Zmat)

    ###glm
    lm1.out <- try(glm(ya ~ pg1m,data=Z1,family = binomial(link = "logit")),silent=TRUE)
    lm0.out <- try(glm(ya ~ pg0m ,data=Z0,family = binomial(link = "logit")),silent=TRUE)
    if ('try-error' %in% class(lm1.out)){
      next}
    if ('try-error' %in% class(lm0.out)){
      next}

    glm_pg1 <-1/(1/(exp(zp[,which(!is.na(lm1.out$coefficients))] %*%
                          lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]))+1)
    glm_pg0 <-1/(1/(exp(zp[,which(!is.na(lm0.out$coefficients))] %*%
                          lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]))+1)
    glm_pg<-cbind(glm_pg0,glm_pg1)

    ####ps
    Z3 = data.frame(Xmat,Trt)
    zp1 = model.matrix(Trt ~.,data = Z3)

    glm.out <- try(glm(Trt ~ zp1,data=Z3, family = binomial(link = "logit")),silent=TRUE)
    if ('try-error' %in% class(glm.out)){
      next}
    glm_ps <- predict(glm.out,data.frame(x=zp1),type="response")
    glm_psr=(glm_ps-mean(glm_ps))/sd(glm_ps)



    k=dsmatchATE2(ya,Trt, ps = glm_psr, M=M,model = "psm")$Kiw.ds
    k=(k-1)/M+1

    ate=0
    for(ii in 1:sample_size){
      ate=ate+(Trt[ii]*ya[ii]*k[ii])+((1-(Trt[ii]*k[ii]))*glm_pg[ii,2])-
        (((1-Trt[ii])*k[ii])*ya[ii])-
        ((1-((1-Trt[ii])*k[ii]))*glm_pg[ii,1])
    }
    bs[i]=ate/sample_size
    i=i+1

  }
  return(var(bs))

}
cvtau<-function(A,X,Y,M,sample_size){

  Z = data.frame(X,Y)
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)
  Z1=Z[loc.1,]
  pg1m = model.matrix(Y ~.,data = Z1)[,-1]
  Z0=Z[loc.0,]
  pg0m = model.matrix(Y ~. ,data = Z0)[,-1]
  zp = model.matrix(Y ~. ,data = Z)

  lm1.out <- glm(Y ~ pg1m,data=Z1,family = binomial(link = "logit"))
  lm0.out <- glm(Y ~ pg0m ,data=Z0,family = binomial(link = "logit"))
  glm_pg1 <-1/(1/(exp(zp[,which(!is.na(lm1.out$coefficients))] %*%
                        lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]))+1)
  glm_pg0 <-1/(1/(exp(zp[,which(!is.na(lm0.out$coefficients))] %*%
                        lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]))+1)
  glm_pg<-cbind(glm_pg0,glm_pg1)

  Z3 = data.frame(X,A)
  zp1 = model.matrix(A ~.,data = Z3)
  glm.out <- glm(A ~ zp1,data=Z3, family = binomial(link = "logit"))
  glm_ps <- predict(glm.out,data.frame(x=zp1),type="response")
  glm_psr=(glm_ps-mean(glm_ps))/sd(glm_ps)


  k=dsmatchATE2(Y,A, ps = glm_psr, M=M, model="psm")$Kiw.ds
  k=(k-1)/M+1
  ate=0
  for(ii in 1:sample_size){
    ate=ate+(A[ii]*Y[ii]*k[ii])+((1-(A[ii]*k[ii]))*glm_pg[ii,2])-
      (((1-A[ii])*k[ii])*Y[ii])-
      ((1-((1-A[ii])*k[ii]))*glm_pg[ii,1])
  }
  return(ate)


}
cvnew<-function(A,X,Y,M,var,repeats=25){

  bias1=rep(0,repeats)
  bias0=rep(0,repeats)

  kk=1
  while(kk<=repeats){
    #  print(kk)
    d1=sample(seq(1:length(A)),2867,replace = FALSE)
    atek=try(cvtau(A=A[d1],X=X[d1,],Y=Y[d1],M=M,sample_size=2867),silent = TRUE)
    if ('try-error' %in% class(atek)){
      next}
    tau=atek/2867

    atek=try(cvtau(A=A[-d1],X=X[-d1,],Y=Y[-d1],M=1,sample_size=2868),silent = TRUE)
    if ('try-error' %in% class(atek)){
      next}
    tau0=atek/2868

    bias0[kk]=tau0-tau

    kk=kk+1
  }

  mse0=mean(bias0^2)+var
  return(mse0)


}


# https://hbiostat.org/data/
rhc <- read.csv("~/Desktop/rhc.csv")

### data clean process
#delete missing value
X=rhc[,c("age" , "sex" , "race" , "edu" , "income" , "ninsclas" ,
         "cat1" , "das2d3pc" , "dnr1" , "ca" , "surv2md1" , "aps1" , "scoma1" ,
         "wtkilo1" , "temp1" , "meanbp1" , "resp1" , "hrt1" , "pafi1" , "paco21" ,
         "ph1" , "wblc1" , "hema1" , "sod1" , "pot1" , "crea1" , "bili1" ,
         "alb1" , "resp" , "card" , "neuro" , "gastr" , "renal" , "meta" ,
         "hema" , "seps" , "trauma" , "ortho" , "cardiohx" , "chfhx" , "dementhx" , "psychhx" , "chrpulhx" ,
         "renalhx" , "liverhx" , "gibledhx" , "malighx" , "immunhx" , "transhx" , "amihx")]

sample_size=as.numeric(dim(X)[1])
A=rep(0,sample_size)
for(i in 1:sample_size){
  if(rhc$swang1[i]=="RHC" ){A[i]=1}}
Y=rep(1,sample_size)
for(i in 1:sample_size){
  if(rhc$death[i]=="Yes" ){Y[i]=0}}

Z = data.frame(X,Y)
#n=build model
loc.1 <- which(A == 1)
loc.0 <- which(A == 0)
Z1=Z[loc.1,]
pg1m = model.matrix(Y ~.,data = Z1)[,-1]
Z0=Z[loc.0,]
pg0m = model.matrix(Y ~. ,data = Z0)[,-1]
zp = model.matrix(Y ~. ,data = Z)

lm1.out <- glm(Y ~ pg1m,data=Z1,family = binomial(link = "logit"))
lm0.out <- glm(Y ~ pg0m ,data=Z0,family = binomial(link = "logit"))
glm_pg1 <-1/(1/(exp(zp[,which(!is.na(lm1.out$coefficients))] %*%
                      lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]))+1)
glm_pg0 <-1/(1/(exp(zp[,which(!is.na(lm0.out$coefficients))] %*%
                      lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]))+1)
glm_pg<-cbind(glm_pg0,glm_pg1)

Z3 = data.frame(X,A)
zp1 = model.matrix(A ~.,data = Z3)
glm.out <- glm(A ~ zp1,data=Z3, family = binomial(link = "logit"))
glm_ps <- predict(glm.out,data.frame(x=zp1),type="response")


glm_pg0=(glm_pg0-mean(glm_pg0))/sd(glm_pg0)
glm_pg1=(glm_pg1-mean(glm_pg1))/sd(glm_pg1)
glm_psr=(glm_ps-mean(glm_ps))/sd(glm_ps)
glm_pgr<-cbind(glm_pg0,glm_pg1)

set.seed(4444)
Mset=seq(1,30,2)
Muse=1
Mmin=100000000000
for(mm in 1:15){
  M=Mset[mm]
  var2=bootamw(A,X,Y,M,bn=100)
  r2=cvnew(A,X,Y,M,var=var2)
  #  print(r2)
  if(r2<Mmin){
    Muse=M
    Mmin=r2
  }
}

#
Muse=5
causaleffect=c()
k=dsmatchATE2(Y,A,ps = glm_psr,pg = glm_pgr,M=Muse,model="psm")$Kiw.ds
k=(k-1)/Muse+1

ate=0
for(ii in 1:sample_size){
  ate=ate+(A[ii]*Y[ii]*k[ii])+((1-(A[ii]*k[ii]))*glm_pg[ii,2])-
    (((1-A[ii])*k[ii])*Y[ii])-
    ((1-((1-A[ii])*k[ii]))*glm_pg[ii,1])
}

causaleffect[1]=ate/sample_size
bootr=bootamw(A,X,Y,M=Muse,bn=500)
causaleffect[2]=sqrt(bootr[1])
causaleffect[3]=Muse

# 
for (c in 2:ncol(zp)) {
  cat("\n***** (V", c, ") ", colnames(zp)[c], " *****\n", sep = "")
  res = matrix(0, 3, 2)
  colnames(res) <- c("Before Matching", "After Matching")
  rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
  res[1, 1] <- mean(zp[loc.1, c])
  res[2, 1] <- mean(zp[loc.0, c])
  res[3, 1] <- (res[1, 1] - res[2, 1]) /  sd(zp[, c])
  res[1, 2] <- sum(zp[loc.1, c]*k[loc.1])/sample_size
  res[2, 2] <- sum(zp[loc.0, c]*k[loc.0])/sample_size
  res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(zp[, c])
  print(round(res,3))
}

aa=matrix(rep(1,65*3),ncol  = 5)
for (c in 1:13) {
  for(j in 1:5){
    d=(c-1)*5+j+1
    aa[3*(c-1)+1,j]=colnames(zp)[d]
    aa[3*(c-1)+2,j]=round(mean(zp[loc.1, d])-mean(zp[loc.0, d])/sd(zp[, d]),3)
    aa[3*(c-1)+3,j]=round((sum(zp[loc.1, d]*k[loc.1])/sample_size-
      sum(zp[loc.0, d]*k[loc.0])/sample_size)/sd(zp[, d]),3)
  }
}
aa=as.data.frame(aa)
library(writexl)
write_xlsx(aa,"~/Desktop/rr.xlsx")



#### forest plot

#### make data framework

library(tidyverse)

Group=c(66,1)
for (c in 2:34) {
  d=c-1
  Group[2*d-1]="before matching"
  Group[2*d]="after matching"
}

Variables=c(66,1)
for (c in 2:34) {
  d=c-1
      Variables[2*d-1]=colnames(zp)[c]
      Variables[2*d]=colnames(zp)[c]
}

difference=c(66,1)
LowerLimit=c(66,1)
UpperLimit=c(66,1)
for (c in 2:34) {
  d=c-1
  difference[2*d-1]=round((mean(zp[loc.1, c])-mean(zp[loc.0, c]))/sd(zp[, c]),3)
  difference[2*d]=round((sum(zp[loc.1, c]*k[loc.1])/sample_size-
                           sum(zp[loc.0, c]*k[loc.0])/sample_size)/sd(zp[, c]),3)
  LowerLimit[2*d-1]=round((mean(zp[loc.1, c])-mean(zp[loc.0, c]))/sd(zp[, c]),3)-0.03
  LowerLimit[2*d]=round((sum(zp[loc.1, c]*k[loc.1])/sample_size-
                           sum(zp[loc.0, c]*k[loc.0])/sample_size)/sd(zp[, c]),3)-0.03
  UpperLimit[2*d-1]=round((mean(zp[loc.1, c])-mean(zp[loc.0, c]))/sd(zp[, c]),3)+0.03
  UpperLimit[2*d]=round((sum(zp[loc.1, c]*k[loc.1])/sample_size-
                           sum(zp[loc.0, c]*k[loc.0])/sample_size)/sd(zp[, c]),3)+0.03
}


aa=data.frame(Group,Variables,difference,LowerLimit,UpperLimit)

p = ggplot(data=aa,
           aes(x = Group,y =difference, ymin = LowerLimit, ymax = UpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab('Variable')+ ylab("Standard Difference")+
  geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.1,cex=0.1)+ 
  facet_wrap(~Variables,strip.position="left",nrow=11,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p

Group=c(64,1)
for (c in 35:66) {
  d=c-34
  Group[2*d-1]="before matching"
  Group[2*d]="after matching"
}

Variables=c(64,1)
for (c in 35:66) {
  d=c-34
  Variables[2*d-1]=colnames(zp)[c]
  Variables[2*d]=colnames(zp)[c]
}

####
difference=c(64,1)
LowerLimit=c(64,1)
UpperLimit=c(64,1)
for (c in 35:66) {
  d=c-34
  difference[2*d-1]=round((mean(zp[loc.1, c])-mean(zp[loc.0, c]))/sd(zp[, c]),3)
  difference[2*d]=round((sum(zp[loc.1, c]*k[loc.1])/sample_size-
                           sum(zp[loc.0, c]*k[loc.0])/sample_size)/sd(zp[, c]),3)
  LowerLimit[2*d-1]=round((mean(zp[loc.1, c])-mean(zp[loc.0, c]))/sd(zp[, c]),3)-0.003
  LowerLimit[2*d]=round((sum(zp[loc.1, c]*k[loc.1])/sample_size-
                           sum(zp[loc.0, c]*k[loc.0])/sample_size)/sd(zp[, c]),3)-0.003
  UpperLimit[2*d-1]=round((mean(zp[loc.1, c])-mean(zp[loc.0, c]))/sd(zp[, c]),3)+0.003
  UpperLimit[2*d]=round((sum(zp[loc.1, c]*k[loc.1])/sample_size-
                           sum(zp[loc.0, c]*k[loc.0])/sample_size)/sd(zp[, c]),3)+0.003
}


aa=data.frame(Group,Variables,difference,LowerLimit,UpperLimit)

p = ggplot(data=aa,
           aes(x = Group,y =difference, ymin = LowerLimit, ymax = UpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab('Variable')+ ylab("Standard Difference")+
  geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.1,cex=0.1)+ 
  facet_wrap(~Variables,strip.position="left",nrow=16,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=10,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p