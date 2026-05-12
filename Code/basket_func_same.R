sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}


GP.unique.grid = function(xgrid){
  d = ncol(xgrid)
  uqgrid = unique(c(xgrid))
  uqidx = apply(xgrid,1:d,function(x) return(which(uqgrid==x)))
  return(list(uqgrid=uqgrid,uqidx=uqidx))
}

GP.simplex = function(xgrid,d,n){
  simplexList = list()
  d = ncol(xgrid)
  numPoly = 0
  for(i in 0:n) {
    simplexList[[i+1]] = matrix(xsimplex(d,i),nrow=d)
    numPoly = numPoly + ncol(simplexList[[i+1]])
  }
  return(list(simplexList=simplexList,numPoly=numPoly))
}

GP.design.mat = function(xgrid,a=1,b=1,alpha=0.9){
  e.degree= GP.eigen.degree(alpha=alpha,a=a,b=b,d=ncol(xgrid))  
  Xmat = GP.eigen.funcs(xgrid,e.degree$n,a=a,b=b)
  lambda = GP.eigen.value(n=e.degree$n,a=a,b=b,d=ncol(xgrid))
  return(list(Xmat=Xmat,lambda=lambda,e.degree=e.degree))
}



GP.eigen.funcs = function(xgrid,n,a,b){       #for calculating psi
  cn = sqrt(a^2+2*a*b)
  d = ncol(xgrid)
  uqset = GP.unique.grid(xgrid)
  temp = GP.simplex(xgrid,d,n)
  return(GP.eigen.funcs.comp(uqset,xgrid,temp$simplexList,
                             temp$numPoly,n,d,cn))
}

GP.eigen.funcs.comp = function(uqset,xgrid,simplexList,
                               numPoly,n,d,cn,trace=FALSE) {
  sqrt2c = sqrt(2*cn)
  D = (sqrt2c)^(d/2)
  
  hermiteList = hermite.h.polynomials(n+1,normalized=TRUE)
  uqhermite = list()
  for(i in 1:(n+1)) {
    uqhermite[[i]] = as.function(hermiteList[[i]])(sqrt2c*uqset$uqgrid)
  }
  
  eigenfunctions = matrix(NA, nrow=nrow(xgrid),ncol=numPoly)
  i = 1
  for(k in 1:length(simplexList)) {
    for(j in 1:ncol(simplexList[[k]])) {
      hermiteGrid = matrix(NA,nrow=nrow(xgrid),ncol=d)
      for(l in 1:d)
        hermiteGrid[,l] = uqhermite[[simplexList[[k]][l,j]+1]][uqset$uqidx[,l]]      
      eigenfunctions[,i] = D*apply(exp(-cn*xgrid^2)*hermiteGrid,1,prod)
      if(trace==TRUE){
        cat(i,"/",numPoly,"\n")
        flush.console()
      }
      i = i + 1
    }
  }
  return(eigenfunctions)
}

GP.gaussian.kernel = function(x1,x2,a=1,b=1){
  return(exp(-a*sum(x1^2)-a*sum(x2^2)-b*sum(abs(x1-x2)^1.99)))
}


GP.eigen.contr = function(n=10,a=1,b=1,d=2){
  cn = sqrt(a^2+2*a*b)
  A = a+b+cn
  B = b/A  
  return((1-B)^d*sum(exp(lchoose(0:n+d-1,d-1)+0:n*log(B))))
}

GP.eigen.value = function(n=10,a=1,b=1,d=2){
  cn = sqrt(a^2+2*a*b)
  A = a+b+cn
  B = b/A
  idx = c(0, choose(0:n+d,d))
  idxlist = sapply(1:(n+1),function(i) return((idx[i]+1):idx[i+1]))
  k=GP.num.eigen.funs(n=n,d=d)
  value = rep(NA,length=k)
  dvalue = (sqrt(pi/A))^d*B^(1:(n+1))
  for(i in 1:(n+1))
    value[idxlist[[i]]] = dvalue[i]
  return(value)
}

GP.num.eigen.funs = function(n=10,d=2){
  return(choose(n+d,d))
}

GP.margin.likehood = function(ygrid,xgrid,n=30,d=2,a=1,b=1,tau2=1){
  
}

GP.eigen.degree = function(alpha,a=1,b=1,d=2,trace=FALSE){
  nleft=100
  nright = 1
  n = nright
  alpha_n = GP.eigen.contr(n,a=a,b=b,d=d)
  max_iter=100
  iter = 1
  while((abs(alpha_n-alpha)>0.01)& (iter<max_iter) & (n>=nright) & (n<=nleft)){
    n = floor((nright+nleft)*0.5)
    alpha_n = GP.eigen.contr(n,a=a,b=b,d=d)
    if(alpha_n<alpha){
      nright=n
    }
    else{
      nleft=n
    }
    if(trace){
      cat("iter ",iter,": [",nright,",",nleft,"]","alpha_n=",alpha_n,"\n")
      flush.console()
    }
    iter=iter+1
  }
  return(list(n=nleft,k=GP.num.eigen.funs(nleft,d=d), alpha_n=GP.eigen.contr(nleft,a=a,b=b,d=d)))
}


GP.simulate.direct = function(xgrid,a=1,b=1,tau2=1){
  if(is.null(nrow(xgrid)))
    xgrid = matrix(xgrid,nrow=length(xgrid),ncol=1)  
  idxset = expand.grid(1:nrow(xgrid),1:nrow(xgrid))
  covmat = sapply(1:nrow(idxset),FUN = function(i){
    return(GP.gaussian.kernel(xgrid[idxset[i,1],],xgrid[idxset[i,2],],a=a,b=b))})
  covmat = matrix(covmat,nrow=nrow(xgrid),ncol=nrow(xgrid))
  
  ygrid = t(rmvnorm(1,sigma=tau2*covmat))
  
  return(ygrid)
}

GP.fit = function(ygrid,xgrid,a=1,b=1,tau2=1,sigma2=1,alpha=0.9,calpha=0.95){
  e.degree= GP.eigen.degree(alpha=alpha,a=a,b=b,d=ncol(xgrid))
  Xmat = GP.eigen.funcs(xgrid,e.degree$n,a=a,b=b)
  lambda = tau2*GP.eigen.value(n=e.degree$n,a=a,b=b,d=ncol(xgrid))
  W = 1.0/(apply(Xmat^2,2  ,sum)+sigma2/lambda)
  thetahat = W*(t(Xmat)%*%ygrid)
  fmean = Xmat%*%thetahat
  Xmat2 = Xmat*Xmat
  fsd = sqrt(Xmat2%*%W*sigma2)
  CIwidth = qnorm(1-0.5*(1-calpha)/length(ygrid))
  fucl = fmean+CIwidth*fsd
  flcl = fmean-CIwidth*fsd
  return(list(thetahat=thetahat,Xmat=Xmat,fmean=fmean,fsd=fsd,fucl=fucl,flcl=flcl))
}


twofigs.levelplot = function(z1,z2,x,y,titles=c("fig1","fig2"),cols=my.cols){
  pdat = list(z = c(z1,z2),group=rep(titles,each=length(z1)),x=rep(x,times=2),y=rep(y,times=2))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1)
}

threefigs.levelplot = function(z1,z2,z3,x,y,titles=c("fig1","fig2","fig3"),cols=my.cols){
  pdat = list(z = c(z1,z2,z3),group=rep(titles,each=length(z1)),x=rep(x,times=3),y=rep(y,times=3))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1)
}

onefig.levelplot.points = function(z,x,y,pts,titles="fig1",cols=my.cols,pch=19,
                                   pcol="black",pcex=0.5){
  pdat = list(z = z,group=rep(titles,each=length(z)),x=rep(x,times=1),y=rep(y,times=1))
  levelplot(z~x+y|group,data=pdat,col.regions=cols,cuts=length(cols)-1,
            panel=function(...){
              panel.levelplot(...)
              lpoints(pts[,1], pts[,2], 
                      pch=pch,col=pcol,cex=pcex)
            }
  )
  
}

simul.theta = function(sgrid,grid.matrix,alphap=150,alphat=300,alpha=0.1,beta=50,sbeta=0.001,beg=520
                       ,smend=400){
  k = ncol(grid.matrix)
  size2 = nrow(sgrid)
  lcenter = c(-0.5,-0.5)
  lintfun.pop = sapply(1:size2,function(i) alphap*exp(-3*sum((sgrid[i,]-lcenter)^2)))
  rcenter=c(0.5,-0.5)
  rintfun.pop = sapply(1:size2,function(i) alphap*exp(-3*sum((sgrid[i,]-rcenter)^2)))
  intfun.pop=lintfun.pop+rintfun.pop 
  intfun.pop[order(intfun.pop)[1:smend]]=sbeta*intfun.pop[order(intfun.pop)[1:smend]]
  intfun.pop[order(intfun.pop)[beg:size2]]=beta*intfun.pop[order(intfun.pop)[beg:size2]]
  lcenter.type0=c(-0.5,-0.2)
  lintfun.type0=sapply(1:size2,function(i) alphat*exp(-4*sum((sgrid[i,]-lcenter.type0)^2)))
  rcenter.type0=c(0.5,-0.2)
  rintfun.type0=sapply(1:size2,function(i) alphat*exp(-4*sum((sgrid[i,]-rcenter.type0)^2)))
  intfun.type0=lintfun.type0+rintfun.type0
  intfun.type0[order(intfun.type0)[1:smend]]=sbeta*intfun.type0[order(intfun.type0)[1:smend]]
  intfun.type0[order(intfun.type0)[beg:size2]]=beta*intfun.type0[order(intfun.type0)[beg:size2]]
  lcenter.type1=c(-0.5,-0.7)
  lintfun.type1=sapply(1:size2,function(i) alphat*exp(-4*sum((sgrid[i,]-lcenter.type1)^2)))
  rcenter.type1=c(0.5,-0.7)
  rintfun.type1=sapply(1:size2,function(i) alphat*exp(-4*sum((sgrid[i,]-rcenter.type1)^2)))
  intfun.type1=lintfun.type1+rintfun.type1
  intfun.type1[order(intfun.type1)[1:smend]]=sbeta*intfun.type1[order(intfun.type1)[1:smend]]
  intfun.type1[order(intfun.type1)[beg:size2]]=beta*intfun.type1[order(intfun.type1)[beg:size2]]
  grid.inverse=solve(t(grid.matrix)%*%grid.matrix)
  theta.pop=grid.inverse%*%t(grid.matrix)%*%log(intfun.pop)*alpha
  theta.type0=grid.inverse%*%t(grid.matrix)%*%log(intfun.type0)*alpha
  theta.type1=grid.inverse%*%t(grid.matrix)%*%log(intfun.type1)*alpha
  true.theta=matrix(c(theta.pop, rep(theta.type0/2,2), rep(theta.type1/2,2)),nrow=k)
  return(true.theta)
}
compute.intensity  = function(true.theta,covar,sgrid,grid.matrix,num.game,num.type=2){  
  p=dim(covar)[2]-1
  log.intfun = list()
  size2 = nrow(sgrid)
  for(gameth in 1:num.game)
  {      
    z_i=matrix(c(as.numeric(covar[covar$game==gameth,-1])))
    type1=2
    type2=3
    temp0.matrix=cbind(true.theta[,1],true.theta[,type1:type2]%*%z_i)
    log.intfun[[paste(gameth,"0",sep="_")]] = apply(grid.matrix%*%temp0.matrix,1,sum)      
    mtype1=4
    mtype2=5
    temp1.matrix=cbind(true.theta[,1],true.theta[,mtype1:mtype2]%*%z_i)      
    log.intfun[[paste(gameth,"1",sep="_")]] = apply(grid.matrix%*%temp1.matrix,1,sum)
    
  }  
  return(log.intfun)
}




simul.data = function(log.intfun,sgrid,num.game){
  test.data=NULL
  dxdy = (sgrid[1,1]-sgrid[2,1])^2
  size2 = nrow(sgrid)
  for(gameth in 1:num.game){
    
    lambda.miss = sum(exp(log.intfun[[paste(gameth,"0",sep="_")]]))*dxdy
    num.miss=rpois(1,lambda.miss)
    if(num.miss>0){
      max.temp0 = max(log.intfun[[paste(gameth,"0",sep="_")]])
      intfun.temp0 = exp(log.intfun[[paste(gameth,"0",sep="_")]] - max.temp0)
      miss.prob = intfun.temp0/sum(intfun.temp0)
      id_miss=sample(1:size2,size=num.miss,prob=miss.prob,replace=TRUE)
      if(num.miss==1){
        miss.grid=matrix(sgrid[id_miss,],nrow=1)
      }
      else miss.grid=sgrid[id_miss,]
      Game=rep(gameth,num.miss)
      Missmade=c(rep(0,num.miss))
      Angle=c(miss.grid[,1])
      Distance=c(miss.grid[,2])
      miss.data.temp=data.frame(Game,Missmade,Angle,Distance)
      test.data=rbind(test.data,miss.data.temp)
    }
    
    lambda.made = sum(exp(log.intfun[[paste(gameth,"1",sep="_")]]))*dxdy
    num.made=rpois(1,lambda.made)
    if(num.made>0){
      max.temp1 = max(log.intfun[[paste(gameth,"1",sep="_")]])
      intfun.temp1 = exp(log.intfun[[paste(gameth,"1",sep="_")]] - max.temp1)
      made.prob = intfun.temp1/sum(intfun.temp1)
      id_made=sample(1:size2,size=num.made,prob=made.prob,replace=TRUE)  
      if(num.made==1){
        made.grid=matrix(sgrid[id_made,],nrow=1)
      }
      else made.grid=sgrid[id_made,]
      Game=rep(gameth,num.made)
      Missmade=c(rep(1,num.made))
      Angle=c(made.grid[,1])
      Distance=c(made.grid[,2])
      made.data.temp=data.frame(Game,Missmade,Angle,Distance)
      test.data=rbind(test.data,made.data.temp)
    }
    
  }
  
  return(test.data)
}



coef.gt=function(gameth,type,data,covar,gridmat,num.type=2){
  game=unique(data$Game)
  num.game=length(unique(data$Game))
  game.th=game[gameth]
  I.type=diag(1,num.type)
  z_i=matrix(c(as.numeric(covar[covar$game==game.th,-1])))
  inter.grid=data[data$Game==game.th&data$Missmade==(type-1),]
  type.coef=kronecker(I.type[type,],t(z_i))
  phi.coef=matrix(c(1,type.coef),nrow=1)
  mat.coef=kronecker(phi.coef,gridmat)
  return(c(nrow(inter.grid),apply(mat.coef,2,sum)))
}


esimat.ini=function(data, covar, gridmat, dxdy, num.type=2){
  num.game = length(unique(data$Game))
  p = dim(covar)[2] - 1
  dimtheta = 1 + num.type*p
  I.type = diag(1, num.type)
  idx = expand.grid(1:num.type, 1:num.game)
  temp = sapply(1:nrow(idx), function(i) coef.gt(idx[i,2], idx[i,1], data, covar, gridmat, num.type=2))
  nonreq=c()
  for(j in 1:ncol(temp))
  { 
    if (temp[1,j] == 0)
      nonreq=c(nonreq,j)
  }
  if (length(nonreq) > 0) {
    temp = temp[,-nonreq] }
  log.esti = log(temp[1,]/dxdy) 
  coef.matrix = t(temp[-1,])
  L = ncol(coef.matrix)/dimtheta
  theta.vec = Solve(coef.matrix, log.esti)
  theta.initial = matrix(theta.vec, nrow=L, ncol=dimtheta)
  return(theta.initial)
  
}


inter.coef=function(gameth,type,theta,data,covar,gridmat,dxdy){
  game<-unique(data$Game)
  num.game<-length(game)
  game.th<-game[gameth]
  p=dim(covar)[2]-1
  z_i=matrix(c(as.numeric(covar[covar$game==game.th,-1])))
  type1=2+p*(type-1)
  type2=1+p*type
  theta.choose=cbind(theta[,1],theta[,type1:type2]%*%z_i)
  theta.choose=apply(theta.choose,1,sum)
  result.matrix=gridmat%*%theta.choose
  sum.vec.0=exp(result.matrix)*dxdy
  return(sum.vec.0)
}


#theta is a matrix dimension is (1+num.game+num.quarter(p+1)+num.type*(p+1))*L
#compute the sum of X_itj*theta
comp.post.pop = function(theta, data, covar, gridmat, dxdy){
  game = unique(data$Game)
  num.game = length(game)
  idx = expand.grid(1:num.type,1:num.game)
  temp = sapply(1:nrow(idx), function(i) inter.coef(idx[i,2],idx[i,1],theta,data,covar,gridmat, dxdy))
  temp1 = sum(temp)
  return(temp1)
}




comp.post.type = function(type,theta, data, covar, gridmat, dxdy){
  game=unique(data$Game)
  num.game=length(game)
  temp=sapply(1:num.game, function(i) inter.coef(i,type,theta, data, covar, gridmat, dxdy))
  temp1=sum(temp)
  return(temp1)
}

#calculate posterior likelihood for the second part given game and type
post.gts=function(gameth,type,theta,data,covar,n,a,b,num.type=2){
  game<-unique(data$Game)
  num.game<-length(game)
  game.th<-game[gameth]
  p=dim(covar)[2]-1
  z_i=matrix(c(as.numeric(covar[covar$game==game.th,-1])))
  inter.data=data[data$Game==game.th&data$Missmade==(type-1),]
  inter.grid=matrix(c(inter.data$Xdist,inter.data$Ydist),ncol=2)
  if(nrow(inter.grid)>0){
    inter.mat=GP.eigen.funcs(inter.grid,n=n,a=a,b=b)
    type1=2+p*(type-1)
    type2=1+p*type
    theta.choose=cbind(theta[,1],theta[,type1:type2]%*%z_i)
    result.matrix=inter.mat%*%theta.choose
    sum.vec.0=apply(result.matrix,1,sum)
    return(sum(sum.vec.0))  
  }
  else return(0)
}

comp.posts=function(theta,data,num.type=2,n,a,b){
  game=unique(data$Game)
  num.game=length(game)
  idx = expand.grid(1:num.type,1:num.game)
  temp=sapply(1:nrow(idx), function(i) post.gts(idx[i,2],idx[i,1],theta,data,covar,n=n,a=a,b=b,num.type=2))
  temp1=sum(temp)
  return(temp1)
}


#calculate coefficient for update theta in pop

comp.mucoef.pop = function(theta, data, covar, gridmat, dxdy){
  game = unique(data$Game)
  num.game = length(game)
  idx = expand.grid(1:num.type,1:num.game)
  temp = sapply(1:nrow(idx), function(i) inter.coef(idx[i,2],idx[i,1],theta,data,covar, gridmat, dxdy))
  temp1=apply(temp,1,sum)
  temp2=t(gridmat)%*%temp1
  return(temp2)
}




comp.mucoef.tc=function(type,k,theta, data, covar, gridmat, dxdy){
  game=unique(data$Game)
  num.game=length(game)
  p=dim(covar)[2]-1
  temp=sapply(1:num.game, function(i) inter.coef(i,type,theta, data, covar, gridmat, dxdy))
  zm=rep(0,num.game)
  for(gameth in 1:num.game){
    game.th<-game[gameth]
    z_i=c(as.numeric(covar[covar$game==game.th,-1]))
    zm[gameth]=z_i[k]
  }
  temp1=t(gridmat)%*%temp%*%zm
  return(temp1)
}



#calculate mu for specific game i
mu.gt=function(gameth,type,data,n=6,a=1,b=5){
  game=unique(data$Game)
  gameth=game[gameth]
  game.data=data[data$Game==gameth&data$Missmade==type-1,]
  game.grid=matrix(c(game.data$Xdist,game.data$Ydist),ncol=2)
  mu.gameth=apply(GP.eigen.funcs(game.grid,n,a,b),2,sum)
  return(mu.gameth)
}

#calculate mu for specific type
mu.tc=function(type,k,data,covar,n=6,a=1,b=5){
  game=unique(data$Game)
  num.game=length(game)
  temp=sapply(1:num.game, function(i) mu.gt(i,type,data,n,a,b))
  zm=rep(0,num.game)
  for(gameth in 1:num.game){
    game.th<-game[gameth]
    z_i=c(as.numeric(covar[covar$game==game.th,-1]))
    zm[gameth]=z_i[k]
  }
  temp1 =temp%*%zm
  return(temp1)   
} 



#function to plot the picture
plot.shots = function(mydat,game=c(1,2),split="Missmade"){
  cols = c("blue","red","yellow","black")
  #cols = rev(gray(c(0,0.1,0.2,0.4)))
  game.idx = which(is.element(mydat$Game,game))
  miss = which(mydat$Missmade==0)
  made = which(mydat$Missmade==1)
  splittype = unique(mydat[[split]])
  
  split.idx = list()
  for(i in 1:length(splittype)){
    split.idx[[i]] = which(mydat[[split]]==splittype[i])  
  }
  plot(0,0,type="n",
       xlim=c(-1,1),ylim=c(-1,1),xlab="Angle",ylab="Distance",asp=1)
  for(i in 1:length(splittype)){
    
    miss.idx = intersect(intersect(game.idx,split.idx[[i]]),miss)
    made.idx = intersect(intersect(game.idx,split.idx[[i]]),made)
    points(mydat$Angle[miss.idx],
           mydat$Distance[miss.idx],
           pch=4,cex=1,col=cols[i])
    
    points(mydat$Angle[made.idx],
           mydat$Distance[made.idx],
           pch=19,cex=1,col=cols[i])
  }
}



