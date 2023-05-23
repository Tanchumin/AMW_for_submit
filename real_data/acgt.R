library(dsmatch1)
library(gam)



bootamw<-function(Xmator,yaor,Trtor,zor,ps,pg,sample_size,M,bn=100){


  tsize=sum(Trtor)
  Xmat=Xmator
  ya=yaor
  z=zor
  bs=rep(0,bn)
  loc.1or <- which(Trtor == 1)
  loc.0or <- which(Trtor == 0)

  for(i in 1:bn){

    sampledata1=sample(loc.1or,tsize,replace = TRUE)
    sampledata0=sample(loc.0or,sample_size-tsize,replace = TRUE)
    ## make new
    Xmat[1:tsize,]=Xmator[sampledata1,]
    z[1:tsize,]=zor[sampledata1,]
    ya[1:tsize]=yaor[sampledata1]
    # print(Xmator[sampledata0,])
    #print(Xmat[(tsize+1):sample_size,])
    Xmat[(tsize+1):sample_size,]<-Xmator[sampledata0,]
    z[(tsize+1):sample_size,]<-zor[sampledata0,]
    ya[(tsize+1):sample_size]<-yaor[sampledata0]

    if(ps==TRUE){psm=z
    }else{
      psm=Xmat}

    if(pg==TRUE){pgm=z
    }else{pgm=Xmat}

    Trt=rep(0,sample_size)
    Trt[1:tsize]=rep(1,tsize)

    loc.1 <- which(Trt == 1)
    loc.0 <- which(Trt == 0)

    lm1.out <- glm(ya[loc.1] ~ pgm[loc.1, ])
    lm0.out <- glm(ya[loc.0] ~ pgm[loc.0, ])
    glm_pg1 <- cbind(1, pgm)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
    glm_pg0 <- cbind(1, pgm)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    glm_pg<-cbind(glm_pg0,glm_pg1)


    glm.out <- glm(Trt ~ psm, family = binomial(link = "logit"))
    #glm_ps<-cbind(1,Xmat)%*%glm.out$coefficients
    glm_ps <- predict(glm.out,data.frame(x=psm),type="response")

    glm_ps1<- cbind(1, psm)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
    glm_pg0=(glm_pg0-mean(glm_pg0))/sd(glm_pg0)
    glm_pg1=(glm_pg1-mean(glm_pg1))/sd(glm_pg1)
    glm_psr=(glm_ps1-mean(glm_ps1))/sd(glm_ps1)
    glm_pgr<-cbind(glm_pg0,glm_pg1)


    k=dsmatchATE2(ya,Trt, ps = glm_psr, pg = glm_pgr,M=M,model = "psm")$Kiw.ds
    k=(k-1)/M+1

    ate=0
    for(ii in 1:sample_size){
      ate=ate+(Trt[ii]*ya[ii]*k[ii])+((1-(Trt[ii]*k[ii]))*glm_pg[ii,2])-
        (((1-Trt[ii])*k[ii])*ya[ii])-
        ((1-((1-Trt[ii])*k[ii]))*glm_pg[ii,1])
    }
    bs[i]=ate/sample_size

  }
  return(var(bs))

}
bootamw1<-function(Xmator,yaor,Trtor,zor,ps,pg,sample_size,M,bn=100,model="amw"){


  tsize=sum(Trtor)
  Xmat=Xmator
  ya=yaor
  z=zor
  bs=rep(0,bn)
  loc.1or <- which(Trtor == 1)
  loc.0or <- which(Trtor == 0)

  for(i in 1:bn){

    sampledata1=sample(loc.1or,tsize,replace = TRUE)
    sampledata0=sample(loc.0or,sample_size-tsize,replace = TRUE)
    ## make new
    Xmat[1:tsize,]=Xmator[sampledata1,]
    z[1:tsize,]=zor[sampledata1,]
    ya[1:tsize]=yaor[sampledata1]
    # print(Xmator[sampledata0,])
    #print(Xmat[(tsize+1):sample_size,])
    Xmat[(tsize+1):sample_size,]<-Xmator[sampledata0,]
    z[(tsize+1):sample_size,]<-zor[sampledata0,]
    ya[(tsize+1):sample_size]<-yaor[sampledata0]

    if(ps==TRUE){psm=z
    }else{
      psm=Xmat}

    if(pg==TRUE){pgm=z
    }else{pgm=Xmat}

    Trt=rep(0,sample_size)
    Trt[1:tsize]=rep(1,tsize)

    loc.1 <- which(Trt == 1)
    loc.0 <- which(Trt == 0)

    lm1.out <- glm(ya[loc.1] ~ pgm[loc.1, ])
    lm0.out <- glm(ya[loc.0] ~ pgm[loc.0, ])
    glm_pg1 <- cbind(1, pgm)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
    glm_pg0 <- cbind(1, pgm)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    glm_pg<-cbind(glm_pg0,glm_pg1)


    glm.out <- glm(Trt ~ psm, family = binomial(link = "logit"))
    #glm_ps<-cbind(1,Xmat)%*%glm.out$coefficients
    glm_ps <- predict(glm.out,data.frame(x=psm),type="response")

    glm_ps1<- cbind(1, psm)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
    glm_pg0=(glm_pg0-mean(glm_pg0))/sd(glm_pg0)
    glm_pg1=(glm_pg1-mean(glm_pg1))/sd(glm_pg1)
    glm_psr=(glm_ps1-mean(glm_ps1))/sd(glm_ps1)
    glm_pgr<-cbind(glm_pg0,glm_pg1)

    if (model=="amw"){
      k=dsmatchATE2(ya,Trt, ps = glm_psr, pg = glm_pgr,M=M,model = "psm")$Kiw.ds
      k=(k-1)/M+1
      ate=0
      for(ii in 1:sample_size){
        ate=ate+(Trt[ii]*ya[ii]*k[ii])+((1-(Trt[ii]*k[ii]))*glm_pg[ii,2])-
          (((1-Trt[ii])*k[ii])*ya[ii])-
          ((1-((1-Trt[ii])*k[ii]))*glm_pg[ii,1])
      }
      bs[i]=ate/sample_size
    }
    else if (model=="ipw"){
      ate=0
      for(ii in 1:sample_size){
        ate=ate+(Trt[ii]*ya[ii]/glm_ps[ii])-((1-Trt[ii])*ya[ii]/(1-glm_ps[ii]))
      }
      bs[i]=ate/sample_size
    }

    else if (model=="aipw"){
      ate=0
      for(ii in 1:sample_size){
        ate=ate+(Trt[ii]*ya[ii]/glm_ps[ii])+((1-(Trt[ii]/glm_ps[ii]))*glm_pg[ii,2])-
          (((1-Trt[ii])/(1-glm_ps[ii]))*ya[ii])-
          ((1-((1-Trt[ii])/(1-glm_ps[ii])))*glm_pg[ii,1])
      }
      bs[i]=ate/sample_size
    }

    else if (model=="psm"){
      ps1=dsmatchATE(ya, Xmat,Trt, method = "ps",ps = glm_psr)
      bs[i]=ps1$est.ps
    }


  }
  lq=as.numeric(quantile(bs, probs = 0.025, na.rm = FALSE))
  hq=as.numeric(quantile(bs, probs = 0.975, na.rm = FALSE))
  return(c(var(bs),lq,hq))

}
cvnew<-function(ya,Xmat,Trt,z,ps,pg,M,var,repeats=25){

  size=dim(Xmat)[1]
  p=dim(Xmat)[2]

  if(ps==TRUE){psm=z
  }else{
    psm=Xmat}

  if(pg==TRUE){pgm=z
  }else{
    pgm=Xmat}


  bias1=rep(0,repeats)
  bias0=rep(0,repeats)

  for(kk in 1:repeats){
    cvData=cbind(Trt,ya,psm,pgm,Xmat)
    cvData<-cvData[sample(nrow(cvData)),]
    #Create 2 equally size folds
    yourdata <- cut(seq(1,nrow(cvData)),breaks=2,labels=FALSE)
    #Perform 5 fold cross validation

    iii=1
    #Segement your data by fold using the which() function
    testIndexes <- which(yourdata==iii,arr.ind=TRUE)
    testData <- cvData[testIndexes, ]
    trainData <- cvData[-testIndexes, ]

    train_ya=trainData[,2]
    train_Trt=trainData[,1]
    train_loc.1=which(train_Trt == 1)
    train_loc.0=which(train_Trt == 0)
    train_psm=trainData[,3:(3+p-1)]
    train_pgm=trainData[,(3+p):(3+2*p-1)]

    test_ya=testData[,2]
    test_Trt=testData[,1]
    test_loc.1=which(test_Trt == 1)
    test_loc.0=which(test_Trt == 0)
    test_psm=testData[,3:(3+p-1)]
    test_pgm=testData[,(3+p):(3+2*p-1)]
    test_X=testData[,(3+2*p):(3+3*p-1)]



    lm1.out <- glm(train_ya[train_loc.1] ~ train_pgm[train_loc.1, ])
    lm0.out <- glm(train_ya[train_loc.0] ~ train_pgm[train_loc.0, ])
    glm_pg1 <- cbind(1,  train_pgm)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
    glm_pg0 <- cbind(1,  train_pgm)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    glm_pg<-cbind(glm_pg0,glm_pg1)

    glm.out <- glm(train_Trt ~ train_psm, family = binomial(link = "logit"))
    pro=exp(cbind(1,train_psm)%*%glm.out$coefficients)
    glm_ps=pro/(pro+1)
    #  glm_ps <- predict(glm.out,data.frame(x=full_psm),type="response")
    # regularization for scores which can be used in matching
    glm_ps1<- cbind(1, train_psm)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
    glm_pg0=(glm_pg0-mean(glm_pg0))/sd(glm_pg0)
    glm_pg1=(glm_pg1-mean(glm_pg1))/sd(glm_pg1)
    glm_psr=(glm_ps1-mean(glm_ps1))/sd(glm_ps1)
    glm_pgr<-cbind(glm_pg0,glm_pg1)
    k=dsmatchATE2(train_ya,train_Trt, ps = glm_psr, pg = glm_pgr,M=1,model = "psm")$Kiw.ds

    atek0=0
    for(ii in 1:(size/2)){
      atek0=atek0+(train_Trt[ii]*train_ya[ii]*k[ii])+((1-(train_Trt[ii]*k[ii]))*glm_pg[ii,2])-
        (((1-train_Trt[ii])*k[ii])*train_ya[ii])-
        ((1-((1-train_Trt[ii])*k[ii]))*glm_pg[ii,1])
    }
    tau0=atek0/(size/2)

    lm1.out <- glm(test_ya[test_loc.1] ~ test_pgm[test_loc.1, ])
    lm0.out <- glm(test_ya[test_loc.0] ~ test_pgm[test_loc.0, ])
    glm_pg1 <- cbind(1,  test_pgm)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
    glm_pg0 <- cbind(1,  test_pgm)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    glm_pg<-cbind(glm_pg0,glm_pg1)

    glm.out <- glm(test_Trt ~ test_psm, family = binomial(link = "logit"))
    pro=exp(cbind(1,test_psm)%*%glm.out$coefficients)
    glm_ps=pro/(pro+1)
    #  glm_ps <- predict(glm.out,data.frame(x=full_psm),type="response")
    # regularization for scores which can be used in matching
    glm_ps1<- cbind(1, test_psm)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
    glm_pg0=(glm_pg0-mean(glm_pg0))/sd(glm_pg0)
    glm_pg1=(glm_pg1-mean(glm_pg1))/sd(glm_pg1)
    glm_psr=(glm_ps1-mean(glm_ps1))/sd(glm_ps1)
    glm_pgr<-cbind(glm_pg0,glm_pg1)
    k=dsmatchATE2(test_ya,test_Trt, ps = glm_psr, pg = glm_pgr,M=M,model="psm")$Kiw.ds
    k=(k-1)/M+1

    atek=0
    for(ii in 1:(size/2)){
      atek=atek+(test_Trt[ii]*test_ya[ii]*k[ii])+((1-(test_Trt[ii]*k[ii]))*glm_pg[ii,2])-
        (((1-test_Trt[ii])*k[ii])*test_ya[ii])-
        ((1-((1-test_Trt[ii])*k[ii]))*glm_pg[ii,1])
    }
    tau=atek/(size/2)

    bias0[kk]=tau0-tau
  }

  mse0=mean(bias0^2)+var
  return(mse0)


}
cover<-function(true,lb,hb){
  a=0
  if(lb<=true & true<=hb){a=1}
  return(a)
}

##############


library(speff2trial)
data("ACTG175")

### data clean process
#delete missing value
delete_na=which(ACTG175$r==0)
datanew=ACTG175[-delete_na,]
Y = datanew$cd496-datanew$cd40
X = datanew[,c("age","wtkg","hemo","homo","drugs",
               "karnof","oprior","z30","preanti",
               "race","gender","str2","strat","symptom",
               "cd40","cd80")]
X = as.matrix(X)
sample_size=as.numeric(dim(X)[1])
A=rep(0,sample_size)
for(i in 1:sample_size){if(datanew$arms[i]==1 || datanew$arms[i]==2){A[i]=1}}
Z = scale(X)

loc.1 <- which(A == 1)
loc.0 <- which(A == 0)

lm1.out <- glm(Y[loc.1] ~ Z[loc.1, ])
lm0.out <- glm(Y[loc.0] ~ Z[loc.0, ])
glm_pg1 <- cbind(1, Z)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
glm_pg0 <- cbind(1, Z)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
glm_pg<-cbind(glm_pg0,glm_pg1)

glm.out <- glm(A ~ Z, family = binomial(link = "logit"))
#glm_ps<-cbind(1,Xmat)%*%glm.out$coefficients
glm_ps <- predict(glm.out,data.frame(x=Z),type="response")

glm_ps1<- cbind(1, Z)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
glm_pg0=(glm_pg0-mean(glm_pg0))/sd(glm_pg0)
glm_pg1=(glm_pg1-mean(glm_pg1))/sd(glm_pg1)
glm_psr=(glm_ps1-mean(glm_ps1))/sd(glm_ps1)
glm_pgr<-cbind(glm_pg0,glm_pg1)

set.seed(444)
Mset=seq(1,30,1)
Muse=1
Mmin=100000000000
for(mm in 1:30){
  M=Mset[mm]
  # var1=bootamwwild(Xmat,ya,Trt,z,ps,pg,sample_size,M=M,bn=100)
  #  r1=cvnew(ya,Xmat,Trt,z,ps=ps,pg=pg,M=M,var=var1)
  var2=bootamw(Z,Y,A,Z,ps=TRUE,pg=TRUE,sample_size=sample_size,M=M,bn=100)
  r2=cvnew(Y,Z,A,Z,ps=TRUE,pg=TRUE,M=M,var=var2)
  #  print(r2)
  if(r2<Mmin){
    Muse=M
    Mmin=r2
  }
}

#M=12
Muse=12
causaleffect=c()
k=dsmatchATE2(Y,A, ps = glm_psr, pg = glm_pgr,M=Muse,model="psm")$Kiw.ds
k=(k-1)/Muse+1
ate=0
for(ii in 1:sample_size){
  ate=ate+(A[ii]*Y[ii]*k[ii])+((1-(A[ii]*k[ii]))*glm_pg[ii,2])-
    (((1-A[ii])*k[ii])*Y[ii])-
    ((1-((1-A[ii])*k[ii]))*glm_pg[ii,1])
}
causaleffect[1]=ate/sample_size
bootr=bootamw1(Z,Y,A,Z,ps=TRUE,pg=TRUE,sample_size,M=Muse,bn=100,model="amw")
causaleffect[2]=sqrt(bootr[1])
causaleffect[3]=Muse

for (c in 1:ncol(X)) {
  cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
  res = matrix(0, 3, 2)
  colnames(res) <- c("Before Matching", "After Matching")
  rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
  res[1, 1] <- mean(X[loc.1, c])
  res[2, 1] <- mean(X[loc.0, c])
  res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
  res[1, 2] <- sum(X[loc.1, c]*k[loc.1])/sample_size
  res[2, 2] <- sum(X[loc.0, c]*k[loc.0])/sample_size
  res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
  print(round(res,3))
}

Group=rep(1,32)
for (c in 1:16) {
  d=c
  Group[2*d-1]="before"
  Group[2*d]="after"
}

Variables1=rep(1,32)
for (c in 1:16) {
  d=c
  Variables1[2*d-1]=colnames(X)[c]
  Variables1[2*d]=colnames(X)[c]
}

####
difference1=rep(1,32)
LowerLimit1=rep(1,32)
UpperLimit1=rep(1,32)
for (c in 1:16) {
  d=c
  difference1[2*d-1]=abs(round((mean(X[loc.1, c])-mean(X[loc.0, c]))/sd(X[, c]),3))
  difference1[2*d]=abs(round((sum(X[loc.1, c]*k[loc.1])/sample_size-
                           sum(X[loc.0, c]*k[loc.0])/sample_size)/sd(X[, c]),3))
  LowerLimit1[2*d-1]=abs(difference1[2*d-1])-0.001
  LowerLimit1[2*d]=abs(difference1[2*d])-0.001
  UpperLimit1[2*d-1]=abs(difference1[2*d-1])+0.001
  UpperLimit1[2*d]=abs(difference1[2*d])+0.001
}


aa1=data.frame(Group1,Variables1,difference1,LowerLimit1,UpperLimit1)

#p1 = ggplot(data=aa1,
  #         aes(x = Group1,y =difference1, ymin = LowerLimit1, ymax = UpperLimit1 ))+
 # geom_pointrange(aes(col=Group1),size=0.3)+
 # geom_hline(aes(fill=Group1),yintercept =0, linetype=1)+
 # geom_hline(aes(fill=Group1),yintercept =0.2, linetype=2)+
 # geom_hline(aes(fill=Group1),yintercept =-0.2, linetype=2)+
 # xlab('Variable')+ ylab("Difference ACGT")+
#  geom_errorbar(aes(ymin=LowerLimit1, ymax=UpperLimit1,col=Group1),width=0.01,cex=0.01,size=1)+ 
 # facet_wrap(~Variables1,strip.position="left",nrow=8,scales = "free_y") +
 # theme(plot.title=element_text(size=16,face="bold"),
 #       axis.text.y=element_blank(),
 #       axis.ticks.y=element_blank(),
 #       axis.text.x=element_text(face="bold"),
        # xlab ylab
 #       axis.title=element_text(size=7,face="bold"),
  #      legend.key.size = unit(0.1, 'cm'),
  #      strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
  #      legend.position="none")+
  #        coord_flip()


ggplot(data=aa1,
       aes(x = Variables1,y =difference1, ymin = LowerLimit1, ymax = UpperLimit1 ,label=Group))+
geom_pointrange(aes(col=Group)) + 
  geom_errorbar(aes(ymin=LowerLimit1, ymax=UpperLimit1,col=Group1),width=0.01,cex=0.01,size=1)+
  geom_hline(aes(fill=Group),yintercept =0, linetype=1)+
  geom_hline(aes(fill=Group),yintercept =0.05, linetype=2)+
 # coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab('Variable')+ ylab("Absolute Standardized Differences")+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(size=6,face="bold"),
        # xlab ylab
        axis.title=element_text(size=7,face="bold"),
        legend.key.size = unit(0.1, 'cm'),
        strip.text.x = element_text(hjust=0,vjust = 1,angle=90,face="bold")
      )

p1

solve(matrix(c(210.90901164 ,-32.19136858,
-32.19136858,   6.15359157),nrow = 2))
