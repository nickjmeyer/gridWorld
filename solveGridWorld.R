library(doMC)
library(foreach)
library(Rcpp)
library(RcppArmadillo)
registerDoMC(detectCores())

sourceCpp("valueIt.cpp")


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

valueIterIP<-function(v,g,gamma=1.0){
  indV = 1
  for(x in 1:g$x){
    for(y in 1:g$y){
      vA = rep(0,length(g$actions))
      for(a in 1:length(g$actions)){
        T = rep(0,g$x*g$y)
        indT = 1
        for(xp in 1:g$x){
          for(yp in 1:g$y){
            T[indT] = transProb(c(x,y),a,c(xp,yp),g)
            indT = indT + 1
          }
        }
        vA[a] = expReward(c(x,y),a,g) + gamma * sum(T*v)
      }
      v[indV]= max(vA)
      indV = indV + 1
    }
  }
  return(v)
}


valueIterPar<-function(v,g,gamma = 1.0){
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


valueIterFast<-function(v,g,gamma = 1.0){
  return(solveValueIterFast(v,c(g$r),
                            as.integer(g$x),as.integer(g$y),
                            as.integer(g$goal[1]-1),as.integer(g$goal[2]-1),
                            as.double(g$noise),
                            as.integer(g$s[1]-1),as.integer(g$s[2]-1),
                            unlist(g$actions),
                            as.double(1.0),as.double(1e-8)))
}


solveValueIter<-function(g,vInit=NULL,gamma = 1.0,
                         tol=1e-8,verbose=FALSE,method=NULL){
  if(is.null(vInit)){
    vInit = rep(0,g$x*g$y)
  }

  if(is.null(method))
    iterFunc = valueIter
  else if(method == "Par")
    iterFunc = valueIterPar
  else if(method == "IP")
    iterFunc = valueIterIP
  else if(method == "Fast")
    iterFunc = valueIterFast
  else
    stop(paste("Invalid method of",method))

  converged = FALSE
  while(!converged){
    v = iterFunc(vInit,g,gamma)

    converged = (norm(v-vInit,"2") < tol)
    if(verbose)
      cat(paste(sprintf("% 16.8f",sum(v)),"::",
                sprintf("% 16.8f\n",norm(v-vInit,"2"))))
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
