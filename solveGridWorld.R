library(doMC)
library(foreach)
registerDoMC(detectCores())

tForPolicy<-function(policy,g){
  T = matrix(0,g$x*g$y,g$x*g$y)
  indX = 1
  for(x in 1:g$x){
    for(y in 1:g$y){
      indY = 1
      for(xp in 1:g$x){
        for(yp in 1:g$y){
          T[indX,indY] = transProb(c(x,y),policy[x,y],c(xp,yp),g)
          indY = indY + 1
        }
      }
      indX = indX + 1
    }
  }
  return(T)
}

rForPolicy<-function(policy,g){
  R = rep(0,g$x*g$y)
  ind = 1
  for(x in 1:g$x){
    for(y in 1:g$y){
      R[ind] = expReward(c(x,y),policy[x,y],g)
      ind = ind + 1
    }
  }
  return(R)
}

valueIter<-function(v,g,gamma = 1.0){
  vA = lapply(1:length(g$actions),
              function(a){
                policy = matrix(a,g$x,g$y)
                T = tForPolicy(policy,g)
                R = rForPolicy(policy,g)
                return(R + gamma*T%*%v)
              })
  return(do.call(pmax,args=vA))
}


valueIterFast<-function(v,g,gamma = 1.0){
  v = foreach(x=1:g$x,.combine=c)%:%
    foreach(y=1:g$y,.combine=c)%dopar%{
      vA = lapply(1:length(g$actions),
                  function(a){
                    T = rep(0,g$x*g$y)
                    ind = 1
                    for(xp in 1:g$x){
                      for(yp in 1:g$y){
                        T[ind] = transProb(c(x,y),a,c(xp,yp),g)
                        ind = ind + 1
                      }
                    }
                    return(expReward(c(x,y),a,g) + gamma*sum(T*v))
                  })
      return(max(unlist(vA)))
    }
  return(v)
}


solveValueIter<-function(g,vInit=NULL,gamma = 1.0,
                         tol=1e-8,verbose=FALSE,fast=TRUE){
  if(is.null(vInit)){
    vInit = rep(0,g$x*g$y)
  }
  converged = FALSE
  while(!converged){
    if(fast)
      v = valueIterFast(vInit,g,gamma)
    else
      v = valueIter(vInit,g,gamma)

    converged = (sum((v-vInit)^2) < tol)
    if(verbose)
      cat(paste(sprintf("% 16.8f",sum(v)),"::",
                sprintf("% 16.8f\n",sum((v-vInit)^2))))
    vInit = v
  }
  return(v)
}


polEval<-function(policy,g,gamma=1.0){
  T = tForPolicy(policy,g)
  R = rForPolicy(policy,g)

  v = solve(diag(length(R)) - gamma*T,R)
  return(v)
}


polUpdate<-function(v,g,gamma=1.0,ties.method="first"){
  vA = sapply(1:length(g$actions),
              function(a){
                policy=matrix(a,g$x,g$y)
                T = tForPolicy(policy,g)
                R = rForPolicy(policy,g)
                return(R + gamma*T%*%v)
              })
  return(matrix(max.col(vA,ties.method=ties.method),g$x,g$y,byrow=TRUE))
}


solvePolIter<-function(g,polInit=NULL,gamma = 1.0){
  if(is.null(polInit)){
    polInit = matrix(sample(1:length(g$actions),
                            g$x*g$y,replace=TRUE),g$x,g$y)
  }

  converged = FALSE
  while(!converged){
    v = polEval(polInit,g,gamma)
    policy = polUpdate(v,g,gamma)

    converged = all(policy == polInit)

    polInit = policy
  }
  return(policy)
}
