rm(list=ls(all=TRUE))

source("gridWorld.R")

g = newGrid()

## plot grid world
plotGrid(g)

## plot the state
plotGrid(g,showState = TRUE)

## transition to a new state and plot it
for(i in 1:10){
  Sys.sleep(1)
  g = nextStep(sample(1:4,1),g)$g
  plotGrid(g,showState = TRUE)
}


## a randomly generated policy
policy = matrix(sample(1:4,g$x*g$y,replace=TRUE),nrow=g$x,ncol=g$y)

## plot the policy
plotGrid(g,policy=policy)
