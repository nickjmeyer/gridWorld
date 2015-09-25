library(Rcpp)
library(RcppArmadillo)
source("gridWorld.R")

sourceCpp("valueIt.cpp")

g = newGrid()


vOptF = solveValueIterFast(rep(0,g$x*g$y),c(g$r),
                           as.integer(g$x),as.integer(g$y),
                           as.integer(g$goal[1]-1),as.integer(g$goal[2]-1),
                           as.double(g$noise),
                           as.integer(g$s[1]-1),as.integer(g$s[2]-1),
                           unlist(g$actions),
                           as.double(1.0),as.double(1e-8))


for(y in 1:g$y){
  for(x in 1:g$x){
    print(expReward(c(x,y),1,g))
  }
}
