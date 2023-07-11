require(mvtnorm);require(JPEN);require(glasso);require(PDSCE);require(flare)
# methods included: JPEN, GLASSO, PDSCE, BLT, LW
# lamvec modified
lamvec=function(Str)
{
  l2=5; lambda.min.ratio = 0.4;lambda.max.tmp1 = min(max(Str-diag(diag(Str))), -min(Str-diag(diag(Str))))
  lambda.max.tmp2 = max(max(Str - diag(diag(Str))), -min(Str - diag(diag(Str))))
  if (lambda.max.tmp1 == 0) 
     lambda.max = lambda.max.tmp2
    else lambda.max = lambda.max.tmp1
    lambda.min = lambda.min.ratio * lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min),length = l2));
   return(lambda)
}

# trace function.
tr=function(A)
{return(sum(diag(A)))};

# Loss functions.
loss1=function(Sig,Est) # here Omg is inverse and Est is cov matrix.
{
  Arb=tr(solve(Sig)%*%Est)-log(max(det(solve(Sig)%*%Est),10^(-100)))-dim(Sig)[1];
  return(Arb)
}
# loss 2
loss2=function(Sig,Est)
{
  Arb=solve(Sig)%*%Est-diag(dim(Sig)[1]);return(tr(Arb^2)); 
}
# code for mul normal lok like
loglik=function(Y,sigma)
{
  n=dim(Y)[1];p=dim(Y)[2];
  Q=0;
  Q=Y%*%solve(sigma)%*%t(Y);
  Q=0.5*sum(diag(Q));e11=eigen(sigma);ldet=sum(log(e11$values));
  Q=-(n/2)*ldet-n*p/2*log(2*pi)-Q;
  Q=-Q;
  return(Q)
}

#### return the orginal indicies of sorted data #########

whichpart <- function(x, n) {  x=-x;  nx <- length(x);  p <- nx-n;  xp <- sort(x, partial=p)[p];  which(x >= xp);}

# data splitiing in k-fold
f.k.fold <- function(Nobs,K=5){
  rs <- runif(Nobs);  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

#1. JPEN 
jpen_h0=function(Ytr,lambda00,gamma00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  n=nrow(Ytr);m=ncol(Ytr); S1=cov(Ytr);S1=S1+0.0001*diag(m);d0=diag(chol(S1)^2);D11=diag(sqrt(diag(S1)));R1=cov2cor(S1);
  new_cov00=jpen(R1,gam=gamma00,lam=lambda00);new_cov0=D11%*%new_cov00%*%D11;   d=d0;d_old=5*d0;
  while(sum(abs(d-d_old))>.1)
    {
      d_old=d; S1=new_cov0;D11=diag(sqrt(diag(S1)));R1=cor(Ytr); new_cov00=jpen(R1,gam=gamma00,lam=lambda00);new_cov0=D11%*%new_cov00%*%D11; d=diag(chol(new_cov0)^2);
    }
    return(as.matrix(new_cov0))
}

# tuning parameter selection

jpen.tune0=function(Ytr, gama=NULL,lambda=NULL)
{
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];  Str=var(Ytr);R=cov2cor(Str);
  if(is.null(lambda)){lambda = lamvec(R);l2=5} else {l2=length(lambda)}
  if(is.null(gama)){L=5; errl=matrix(0,ncol=3,nrow=L*l2)}   else {L = length(gama); errl = matrix(0, ncol = 3, nrow = L * l2); }
  for (l in 1:l2) {
    st = (l - 1) * L + 1;en = st + L - 1;
    if(is.null(lambda)) { errl[st:en, 1] = 2*lambda[l];c=lambda[l]} else {errl[st:en,1]=lambda[l];c=lambda[l]/2}
    if(is.null(gama)) { c1=(seq(1:5)/5)*c/(1+pi*sqrt(log(p)/n));gamm1=(c/c1-1);     errl[st:en, 2] = c(gamm1);}
    else {errl[st:en,2]=gama}
    }
  for (r in 1:(L * l2)) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];       gamma00 = errl[r, 2];       Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       
      est1=jpen_h0(Ytr1,lambda00,gamma00);  ei1=eigen(est1)$values;
      if (min(ei1) < 1e-06) 
      {
        errl[r, 3] = 10^20
      } else 
      {
        errl[r, 3] =errl[r,3] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=est1,log=TRUE));
      }
    }
  }
  k = which.min(errl[, 3]);   opt = errl[k, -3];
  return(c(opt))
}


#2. GLASOO 
glasso_h0=function(chrom2,lambda00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  n=nrow(chrom2);m=ncol(chrom2); Ytr=chrom2; S1=cov(Ytr);D11=diag(sqrt(diag(S1)));R1=cor(Ytr);
  new_cov00=glasso(R1,rho=lambda00)$w;new_cov0=D11%*%new_cov00%*%D11; return(as.matrix(new_cov0))
}

# tuning parameter selection

glasso.tune0=function(Ytr, lambda=NULL)
{
  p = dim(Ytr)[2];Str=var(Ytr);
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];
  if(length(lambda)==0)   {   lambda=lamvec(cov2cor(Str))}
  L=length(lambda);  errl=matrix(0,ncol=2,nrow=L);
  errl[,1]=lambda;
  for (r in 1:L) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];      Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       est1=glasso_h0(Ytr1,lambda00);cov0=est1;
      errl[r, 2] =errl[r,2] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
    }
  k = which.min(errl[, 2]);   opt = errl[k,1];
  return(c(opt))
}

#3. PDSCE
pdsce_h0=function(chrom2,lambda00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  n=nrow(chrom2);m=ncol(chrom2); Ytr=chrom2;   S1=cov(Ytr);D11=diag(sqrt(diag(S1)));R1=cor(Ytr);
  new_cov00=pdsoft(R1,lambda00)$sigma;new_cov0=D11%*%new_cov00%*%D11;return(as.matrix(new_cov0))
}

# tuning parameter selection
pdsce.tune0=function(Ytr, lambda=NULL)
{
  p = dim(Ytr)[2];Str=var(Ytr);
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];
  if(length(lambda)==0)   {   lambda=lamvec(cov2cor(Str))}
  L=length(lambda);  errl=matrix(0,ncol=2,nrow=L);
  errl[,1]=lambda;
  for (r in 1:L) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];      Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       est1=pdsce_h0(Ytr1,lambda00);cov0=est1;
      errl[r, 2] =errl[r,2] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
    }
  k = which.min(errl[, 2]);   opt = errl[k,1];
  return(c(opt))
}

#4. BLT_thresh 
BLT_h0=function(chrom2,s00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  n=nrow(chrom2);m=ncol(chrom2); Ytr=chrom2;   S1=cov(Ytr);
  new_cov0= abs(sign(sign(S1)*pmax(abs(S1)-s00,0)))*S1;
  Arb=eigen(new_cov0);V=Arb$vectors;Dv=diag(Arb$values);
  if(min(Arb$values)<0.01)
    {
      diag(Dv)=pmax(diag(Dv),0.01);new_cov0=V%*%Dv%*%t(V);
    }
  return(as.matrix(new_cov0))
}

# tuning parameter selection
BLT.tune0=function(Ytr, lambda=NULL)
{
  p = dim(Ytr)[2];Str=var(Ytr);
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];
  if(length(lambda)==0)   {   lambda=lamvec(cov2cor(Str))}
  L=length(lambda);  errl=matrix(0,ncol=2,nrow=L);
  errl[,1]=lambda;
  for (r in 1:L) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];      Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       est1=BLT_h0(Ytr1,s00=lambda00);cov0=est1;
      errl[r, 2] =errl[r,2] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
  }
  k = which.min(errl[, 2]);   opt = errl[k,1];
  return(c(opt))
}

#5. Ledoit-Wolf
lw_h0=function(chrom2,s00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  n=nrow(chrom2);m=ncol(chrom2);Ytr=chrom2;  cov0=cov(Ytr); shrinkage=max(min(s00,1),0);
  meanvar=mean(diag(cov0));prior=meanvar*diag(p);new_cov0=shrinkage*prior+(1-shrinkage)*cov0;
  return(as.matrix(new_cov0))
}

# tuning parameter selection
lw.tune0=function(Ytr, lambda=NULL)
{
  p = dim(Ytr)[2];Str=var(Ytr);
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];
  if(length(lambda)==0)   {   lambda=lamvec(cov2cor(Str))}
  L=length(lambda);  errl=matrix(0,ncol=2,nrow=L);
  errl[,1]=lambda;
  for (r in 1:L) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];      Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       est1= lw_h0(Ytr1,lambda00);cov0=est1;
      errl[r, 2] =errl[r,2] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
    }
  k = which.min(errl[, 2]);   opt = errl[k,1];
  return(c(opt))
}

############################################################
######### for INVERSE COV matrix estimation  ###############
############################################################
# Warm Start estimator for JPEN estimation
jpen.inv.init=function(Ytr)   # returns the JPEN cov matrix
{
  p=dim(Ytr)[2];n=nrow(Ytr);R=cor(Ytr);  c=max(abs(R)-diag(diag(R)));c1=(seq(1:5)/5)*c/(1+pi*sqrt(log(p)/n));gamm1=(c/c1-1);lamb1=c/2;
  jp_opt1=jpen.tune0(Ytr,gama=gamm1,lambda=lamb1);S1=jpen_h0(Ytr,gamma00=jp_opt1[2],lambda00=jp_opt1[1]);  S1[abs(S1)<1e-6]=0;
  return(S1);
}

# JPEN INV estimation
jpen.inv0=function(R1i,gam,lam=NULL) # 
{
  if(length(lam)==0) {lam=2*gam/p};
  E=(R1i+gam*diag(p))/(1+gam);   E1=sign(E)*pmax(abs(E)-lam/(2*(1+gam)),0);diag(E1)=diag(R1i);E1[abs(E1)<1e-8]=0;
  return(E1); 
}

# tuning parameter selection
jpen.inv.tune0=function(Ytr, gama=NULL,lambda=NULL)
{
  K = 4;   n = dim(Ytr)[1];p=dim(Ytr)[2];Str=cov(Ytr);R1=cov2cor(jpen.inv.init(Ytr));R1i=solve(R1);R1i[abs(R1i)<1e-8]=0.0;
  est00=jpen.inv.init(Ytr);
  if(is.null(lambda)){lambda = lamvec(R1i);l2=5} else {l2=length(lambda)}
  if(is.null(gama))   { 
    L=5; errl=matrix(0,ncol=3,nrow=L*l2)  }   else   {     L = length(gama); errl = matrix(0, ncol = 3, nrow = L * l2)   }
  for (l in 1:l2) {
    st = (l - 1) * L + 1;en = st + L - 1;
    if(is.null(lambda)) { errl[st:en, 1] = 2*lambda[l];c=lambda[l]} else {errl[st:en,1]=lambda[l];c=lambda[l]/2}
    if(is.null(gama)) {c1=(seq(1:5)/5)*c/(1+pi*sqrt(log(p)/n));gamm1=(c/c1-1);     errl[st:en, 2] = c(gamm1);}
    else {errl[st:en,2]=gama}
  }
  for (r in 1:(L * l2)) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];       gamma00 = errl[r, 2];       Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,]; 
      D11=diag(sqrt(diag(cov(Ytr1))));est0=jpen.inv0(solve(cov2cor(est00)),lam=lambda00,gam=gamma00);cov0=D11%*%solve(est0)%*%D11;
      errl[r, 3] =errl[r,3] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
  }
  k = which.min(errl[, 3]);   opt = errl[k, -3];
  return(c(opt))
}


#2. GLASOO 
glasso.inv.h0=function(Ytr,lambda00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  n=nrow(Ytr);m=ncol(Ytr); S1=cov(Ytr);D11=diag(sqrt(diag(S1)));R1=cor(Ytr);
  new_cov00i=glasso(R1,rho=lambda00)$wi;new_cov0i=solve(D11)%*%new_cov00i%*%solve(D11); new_covi=0.5*(new_cov0i+t(new_cov0i));return(as.matrix(new_cov0i))
}

# tuning parameter selection
glasso.inv.tune0=function(Ytr, lambda=NULL)
{
  p = dim(Ytr)[2];Str=var(Ytr);
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];
  if(length(lambda)==0)   {   lambda=lamvec(Str)}
  L=length(lambda);  errl=matrix(0,ncol=2,nrow=L);
  errl[,1]=lambda;
  for (r in 1:L) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];      Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       est1i=glasso.inv.h0(Ytr1,lambda00);est1i=0.5*(est1i+t(est1i));cov0=solve(est1i);
      errl[r, 2] =errl[r,2] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
  }
  k = which.min(errl[, 2]);   opt = errl[k,1];
  return(c(opt))
}

#3. PDSCE
pdsce.inv.h0=function(chrom2,lambda00) # given mean, estimate sigma_h0  # part of EM algorithm
{
  new_cov0i=pdsoft(cov(chrom2),lambda00)$omega;new_covi=0.5*(new_cov0i+t(new_cov0i));return(as.matrix(new_cov0i))
}

# tuning parameter selection
pdsce.inv.tune0=function(Ytr, lambda=NULL)
{
  p = dim(Ytr)[2];Str=var(Ytr);
  K = 5;   n = dim(Ytr)[1];   p = dim(Ytr)[2];
  if(length(lambda)==0)   {   lambda=lamvec(Str)}
  L=length(lambda);  errl=matrix(0,ncol=2,nrow=L);
  errl[,1]=lambda;
  for (r in 1:L) 
  {
    id=f.k.fold(n,K);
    for(k in 1:K)
    {
      lambda00 = errl[r, 1];      Ytr1=Ytr[id[[k]]$train,];Yts1=Ytr[id[[k]]$test,];       est1i=pdsce.inv.h0(Ytr1,lambda00);cov0=solve(est1i);
      errl[r, 2] =errl[r,2] -2*sum(dmvnorm(Yts1,mean=rep(0,p),sigma=cov0,log=TRUE));
    }
  }
  k = which.min(errl[, 2]);   opt = errl[k,1];
  return(c(opt))
}

#4. CLIME 

clime.inv=function(Ytr)
{
  return(sugm.select(sugm(data=Ytr,method="clime"),criterion="cv")$opt.icov);
}


