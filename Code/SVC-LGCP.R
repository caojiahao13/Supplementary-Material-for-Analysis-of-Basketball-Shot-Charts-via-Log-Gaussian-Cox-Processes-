den.norm=function(theta,var){
  return(0.5*t(theta)%*%solve(var)%*%theta)
}

adjust_acceptance=function(accept,sgm,target = 0.5){
  y = 1. + 1000.*(accept-target)*(accept-target)*(accept-target)
  if (y < .9)
    y = .9
  if (y > 1.1)
    y = 1.1
  sgm = sgm* y  
  return(sgm)  
}


fit_SVCLGCP = function(shot_data, covar){
    source("../basketball_data_code/curry/basket_func_same.R")
    curry = shot_data
    
    
#     p=dim(covar)[2]-1
#     courtgrid = as.matrix(expand.grid(seq(-1,1,length=50),seq(-1,1,length=40)))
#     intergrid = as.matrix(expand.grid(seq(-1,1,length=60),seq(-1,1,length=60)))
#     dxdy=(intergrid[1,1]-intergrid[2,1])^2
#     a = 0.25
#     d = 2
#     b = 1.5
#     #GP.eigen.degree(0.9,a,b=9,d=2)
#     GP.eigen.degree(0.8,a=a,b=b,d=d)
#     n = GP.eigen.degree(0.8,a=a,b=b,d=d)$n
#     eigmat = GP.eigen.value(n=n,a=a,b=b,d=d)
#     eigmat = diag(eigmat)
#     pmat = eigmat
#     L = nrow(eigmat)

#     intermat = GP.eigen.funcs(intergrid, n=n, a=a, b=b)
#     courtmat = GP.eigen.funcs(courtgrid, n=n, a=a, b=b)



    # num.type = 2
    # apop = 5
    # bpop = 5
    # sigma.pop = 1#rinvgamma(1,apop,bpop)
    # atype = 5
    # btype = 5
    # sigma.type = rep(1,num.type*p)#rinvgamma(num.type*p,atype,btype)
    # N = 15000

    game = unique(curry$Game)
    num.game = length(game)



    dimtheta = 1+num.type*p
    Theta.result = array(0, dim=c(N,L,dimtheta))
    tau.pop = 0.03/(L*dimtheta)^(1/3)
    tauk = rep(0.03/(L*dimtheta)^(1/3), 2*p)
    var.pop.star = tau.pop*diag(1,L)


    mu0.grid = matrix(c(curry$Xdist,curry$Ydist),ncol=2)
    mu_pop = apply(GP.eigen.funcs(mu0.grid, n, a, b),2,sum)
    mutc = matrix(0,nrow=L, ncol=num.type*p)
    for(type in 1:2){
      for(k in 1:p) {
        jk=p*(type-1)+k
        mutc[,jk]= mu.tc(type, k, curry, covar, n=n, a=a, b=b)
      }
    }

    #mu_type=sapply(1:num.type, function(i) mu.type.th(i,curry,covar,L,n=n,a=a,b=b))
    edxdy = (courtgrid[1,1]-courtgrid[2,1])^2
    current.theta =esimat.ini(curry, covar, courtmat, edxdy, num.type=2)
    likelihood1 = rep(0,N)
    likelihood2 = rep(0,N)
    accept.pop = 0
    accept.type = rep(0,p*num.type)
    Sigmapop = rep(0,N)
    Sigmatype = matrix(0, nrow=N, ncol=p*num.type)

    for(sim in 1:N){

      #update theta0
      temp.theta = current.theta

      var.pop = sigma.pop*pmat
      mu.pop.new = current.theta[,1] + tau.pop/2*(matrix(c(-comp.mucoef.pop(current.theta,
                                                                            curry,covar, intermat, dxdy) + mu_pop),ncol=1) - solve(var.pop)%*%current.theta[,1])
      temp.theta[,1] = mu.pop.new+rnorm(L, 0, sd=sqrt(tau.pop)) 
      mu.pop.old = temp.theta[,1] + tau.pop/2*(matrix(c(-comp.mucoef.pop(temp.theta,
                                                                         curry,covar, intermat, dxdy) + mu_pop),ncol=1) - solve(var.pop)%*%temp.theta[,1])
      current.post.pop = comp.post.pop(current.theta,curry,covar,intermat, dxdy) -  #current likelihood value
        mu_pop%*%current.theta[,1] + den.norm(current.theta[,1], var.pop)
      new.post.pop = comp.post.pop(temp.theta,curry,covar,intermat, dxdy) -
        mu_pop%*%temp.theta[,1] + den.norm(temp.theta[,1], var.pop)
      pop.diff.new=temp.theta[,1]-mu.pop.new
      pop.diff.old=current.theta[,1]-mu.pop.old
      pop.prop.new=0.5*sum(pop.diff.new^2/tau.pop)
      pop.prop.old=0.5*sum(pop.diff.old^2/tau.pop)
      temp.pop.post=current.post.pop-new.post.pop+pop.prop.new-pop.prop.old
      u=runif(1)
      if(temp.pop.post>=log(u)){
        current.theta=temp.theta
        accept.pop=accept.pop+1
      }

      a.pop = apop + L*0.5
      b.pop = bpop + den.norm(current.theta[,1], pmat)
      sigma.pop=rinvgamma(1,a.pop,b.pop)
      Sigmapop[sim]=sigma.pop

      for(type in 1:2){
        for(k in 1:p) {
          jk=p*(type-1)+k+1
          sigmajk = sigma.type[jk-1]
          jkmat = eigmat
          varjk = sigmajk*jkmat
          temp.theta = current.theta
          new.mean = current.theta[,jk]+tauk[k]/2*(-comp.mucoef.tc(type,k,current.theta,curry,covar, intermat, dxdy)+
                                                     mutc[,jk-1]-solve(varjk)%*%current.theta[,jk])
          temp.theta[,jk]= new.mean + rnorm(L,0,sd=sqrt(tauk[k])) 
          curr.mean = temp.theta[,jk]+tauk[k]/2*(-comp.mucoef.tc(type,k,temp.theta,curry,covar,intermat, dxdy)+
                                                   mutc[,jk-1] - solve(varjk)%*%temp.theta[,jk])
          new.post.tc = comp.post.type(type,temp.theta,curry,covar,intermat, dxdy) -
            sum(mutc[,jk-1]*temp.theta[,jk]) + den.norm(temp.theta[,jk],varjk)
          curr.post.tc = comp.post.type(type,current.theta,curry,covar,intermat, dxdy)-
            sum(mutc[,jk-1]*current.theta[,jk]) + den.norm(current.theta[,jk],varjk)
          jkdiff.new=temp.theta[,jk]-new.mean
          jkdiff.old=current.theta[,jk]-curr.mean
          jkprop.new=0.5*sum(jkdiff.new^2/tauk[k])
          jkprop.old=0.5*sum(jkdiff.old^2/tauk[k])
          jkpost=curr.post.tc+jkprop.new-new.post.tc-jkprop.old
          jku=runif(1)
          if(jkpost>=log(jku)){
            current.theta=temp.theta
            accept.type[jk-1]=accept.type[jk-1]+1
          }
          ajk=atype + L*0.5
          bjk=btype + den.norm(current.theta[,jk],jkmat)
          sigma.type[jk-1]=rinvgamma(1,ajk,bjk)
        }
      }
      Sigmatype[sim,]=sigma.type
      #update value of b

      current.post=-comp.post.pop(current.theta, curry,covar,intermat, dxdy)+  #current likelihood value
        mu_pop%*%current.theta[,1]+ sum(mutc*current.theta[, 2:dimtheta])

      likelihood1[sim] = current.post

      if(sim%%100==0){
        cat("accept.pop", accept.pop,"\n")
        cat("accept.type", accept.type,"\n")
      }
      if(sim<=5000&sim%%100==0){
        accept.pop=accept.pop/100
        tau.pop=adjust_acceptance(accept.pop,tau.pop,0.5)
        for(tj in 1:length(accept.type)) {
          acceptt=accept.type[tj]/100
          tauk[tj]=adjust_acceptance(acceptt,tauk[tj],0.5) 
        }
        accept.pop=0
        accept.type=rep(0,p*type)
      }

      Theta.result[sim,,]=current.theta  

      cat("Iter", sim,"\n")

      flush.console()
    }

    print(tau.pop)
    print(tauk)
    print(accept.pop)
    print(accept.type)


    # save(Theta.result,file="curry_homelevel_theta.RData")
    # save(likelihood1,file="curryb_homelevel_like1.RData")
    # save(Sigmapop,file="curryb_homelevel_pop.RData")
    # save(Sigmatype,file="curryb_homelevel_type.RData")

    fit = list(Theta.result = Theta.result, likelihood1 = likelihood1, Sigmapop = Sigmapop, Sigmatype = Sigmatype)
    return(fit)

}