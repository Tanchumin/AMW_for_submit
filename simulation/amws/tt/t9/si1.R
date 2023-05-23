#! /usr/bin/Rscript

library(dsmatch1)

sample_size=1000
p=12
re=100
a0=matrix(data=c(1,1,1,1,-1,-1,-1,-1,1,-1,1,-1)*1,nrow=p)
b0=matrix(data=c(1,1,1,1,-1,-1,-1,-1,1,-1,1,-1)*0.2,nrow=p)
set.seed(4449)

poly3=function(a,n,p){
  z=matrix(data=runif(n*p,-sqrt(3),sqrt(3)),nrow=n)
  z[,1]=exp(a[,1])
  z[,2]=exp(a[,2])
  z[,3]=log((a[,3]+1)^2)
  z[,4]=log((a[,4]+1)^2)
  z[,5]=sin((a[,5]-a[,6]))
  z[,6]=cos((a[,5]+a[,6]) )
  z[,7]=sin(a[,7])
  z[,8]=cos(a[,7]-1)
  
  
  
  for(i in 1:n){
    
    #  if(a[i,7]>=0){
    #   z[i,9]=1}
    #  else{z[i,8]=0}
    
    if(a[i,8]>=0.4){
      z[i,9]=1}
    else{z[i,9]=0}
    
    if(a[i,8]>=-0.4){
      z[i,10]=1}
    else{z[i,10]=0}
    
    if(a[i,9]>=0.3){
      z[i,11]=1}
    else{z[i,11]=0}
    
    if(a[i,10]>=-0.3){
      z[i,12]=1}
    else{z[i,12]=0}
    
  }
  
  for(i in 1:p){
    z[,i]=(z[,i]-mean(z[,i]))/sd(z[,i])
  }
  
  return(z)
  
}

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

ps=FALSE
pg=FALSE
causaleffect=matrix(data=rep(1,21*re),ncol=21)
for(i in 1:re){
  print(i)
  Xmat=matrix(data=runif(sample_size*p,1-sqrt(3),1+sqrt(3)),nrow=sample_size)
  z=poly3(Xmat,sample_size,p)
  
  Trt=rep(1,sample_size)
  for(ii in 1:sample_size){
    odd=exp(-z[ii,]%*%a0)
    prob= 1 / (1 + odd)
    Trt[ii]=rbinom(1,size=1,prob=prob )
  }
  
  
  ya=c()
  for (j in 1:sample_size){
    if(Trt[j]==0){
      #ya[i]=rnorm(1,poly1(Xmat[i,],p),1.5)
      ya[j]=rnorm(1,z[j,]%*%b0,4)
    }
    else{
      #ya[i]=rnorm(1,poly1(Xmat[i,],p)+4,1.5)
      ya[j]=rnorm(1,(z[j,]%*%b0),1)
    }
  }
  
  if(ps==TRUE){psm=z
  }else{
    psm=Xmat}
  
  if(pg==TRUE){pgm=z
  }else{pgm=Xmat}
  
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
  
  Mset=seq(1,19,2)
  Muse=1
  Mmin=100000000000
  for(mm in 1:10){
    M=Mset[mm]
    # var1=bootamwwild(Xmat,ya,Trt,z,ps,pg,sample_size,M=M,bn=100)
    #  r1=cvnew(ya,Xmat,Trt,z,ps=ps,pg=pg,M=M,var=var1)
    var2=bootamw(Xmat,ya,Trt,z,ps,pg,sample_size,M=M,bn=100)
    r2=cvnew(ya,Xmat,Trt,z,ps=ps,pg=pg,M=M,var=var2)
    #  print(r2)
    if(r2<Mmin){
      Muse=M
      Mmin=r2
    }
  }
  # print(Muse)
  
  ##### compute tau var
  
  #amw
  k=dsmatchATE2(ya,Trt, ps = glm_psr, pg = glm_pgr,M=Muse,model="psm")$Kiw.ds
  k=(k-1)/Muse+1
  ate=0
  for(ii in 1:sample_size){
    ate=ate+(Trt[ii]*ya[ii]*k[ii])+((1-(Trt[ii]*k[ii]))*glm_pg[ii,2])-
      (((1-Trt[ii])*k[ii])*ya[ii])-
      ((1-((1-Trt[ii])*k[ii]))*glm_pg[ii,1])
  }
  causaleffect[i,1]=ate/sample_size
  bootr=bootamw1(Xmat,ya,Trt,z,ps,pg,sample_size,M=Muse,bn=100,model="amw")
  causaleffect[i,2]=sqrt(bootr[1])
  causaleffect[i,3]=Muse
  causaleffect[i,4]=cover(0,lb=bootr[2],hb=bootr[3])
  causaleffect[i,5]=cover(0,lb=(-1.96*sqrt(bootr[1])+ate/sample_size),
                          hb=(1.96*sqrt(bootr[1])+ate/sample_size))
  
  #ipw
  ate=0
  for(ii in 1:sample_size){
    ate=ate+(Trt[ii]*ya[ii]/glm_ps[ii])-((1-Trt[ii])*ya[ii]/(1-glm_ps[ii]))
  }
  causaleffect[i,6]=ate/sample_size
  bootr=bootamw1(Xmat,ya,Trt,z,ps,pg,sample_size,M=Muse,bn=100,model="ipw")
  causaleffect[i,7]=sqrt(bootr[1])
  causaleffect[i,8]=cover(0,lb=bootr[2],hb=bootr[3])
  causaleffect[i,9]=cover(0,lb=(-1.96*sqrt(bootr[1])+ate/sample_size),
                          hb=(1.96*sqrt(bootr[1])+ate/sample_size))
  
  #aipw
  ate=0
  for(ii in 1:sample_size){
    ate=ate+(Trt[ii]*ya[ii]/glm_ps[ii])+((1-(Trt[ii]/glm_ps[ii]))*glm_pg[ii,2])-
      (((1-Trt[ii])/(1-glm_ps[ii]))*ya[ii])-
      ((1-((1-Trt[ii])/(1-glm_ps[ii])))*glm_pg[ii,1])
  }
  causaleffect[i,10]=ate/sample_size
  bootr=bootamw1(Xmat,ya,Trt,z,ps,pg,sample_size,M=Muse,bn=100,model="aipw")
  causaleffect[i,11]=sqrt(bootr[1])
  causaleffect[i,12]=cover(0,lb=bootr[2],hb=bootr[3])
  causaleffect[i,13]=cover(0,lb=(-1.96*sqrt(bootr[1])+ate/sample_size),
                           hb=(1.96*sqrt(bootr[1])+ate/sample_size))
  
  #psm
  ps1=dsmatchATE(ya, Xmat,Trt, method = "ps",ps = glm_psr)
  causaleffect[i,14]=ps1$est.ps
  bootr=bootamw1(Xmat,ya,Trt,z,ps,pg,sample_size,M=Muse,bn=100,model="psm")
  causaleffect[i,15]=sqrt(bootr[1])
  causaleffect[i,16]=cover(0,lb=bootr[2],hb=bootr[3])
  causaleffect[i,17]=cover(0,lb=(-1.96*sqrt(bootr[1])+ate/sample_size),
                           hb=(1.96*sqrt(bootr[1])+ate/sample_size))
  
  #AMW2
  k=dsmatchATE2(ya,Trt, ps = glm_psr, pg = glm_pgr,M=1,model="psm")$Kiw.ds
  k=(k-1)/Muse+1
  ate=0
  for(ii in 1:sample_size){
    ate=ate+(Trt[ii]*ya[ii]*k[ii])+((1-(Trt[ii]*k[ii]))*glm_pg[ii,2])-
      (((1-Trt[ii])*k[ii])*ya[ii])-
      ((1-((1-Trt[ii])*k[ii]))*glm_pg[ii,1])
  }
  causaleffect[i,18]=ate/sample_size
  bootr=bootamw1(Xmat,ya,Trt,z,ps,pg,sample_size,M=1,bn=100,model="amw")
  causaleffect[i,19]=sqrt(bootr[1])
  causaleffect[i,20]=cover(0,lb=bootr[2],hb=bootr[3])
  causaleffect[i,21]=cover(0,lb=(-1.96*sqrt(bootr[1])+ate/sample_size),
                          hb=(1.96*sqrt(bootr[1])+ate/sample_size))
}
write.csv(causaleffect,file="ac")


