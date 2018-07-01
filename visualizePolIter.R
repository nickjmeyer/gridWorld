rm(list=ls(all=TRUE))

source("gridWorld.R")
source("solveGridWorld.R")

g = newGrid()

## a randomly generated policy
polInit = matrix(sample(1:4,g$x*g$y,replace=TRUE),nrow=g$x,ncol=g$y)

## plot the policy
p = plotGrid(g,policy=polInit)

ggsave("policy_iteration_000.png", p)

converged = FALSE

count = 1

gamma = 0.99
while(!converged){
    v = polEval(polInit,g,gamma)
    policy = polUpdate(v,g,gamma)

    ## plot the policy
    p = plotGrid(g,policy=policy)
    ggsave(sprintf("policy_iteration_%03d.png", count), p)
    count = count + 1

    converged = all(policy == polInit)

    polInit = policy
}
