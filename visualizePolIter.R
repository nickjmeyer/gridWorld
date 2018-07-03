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

g = newGrid()
p = plotGrid(g, showState = TRUE)
iter = 0
ggsave(sprintf("demo_%03d.png", iter), p)
while(!(g$s[1] == 8 && g$s[2] == 9))
{
    print(policy[g$s[1], g$s[2]])
    g = nextStep(policy[g$s[1], g$s[2]], g)$g
    p = plotGrid(g, showState = TRUE)
    print(g$s)
    iter = iter + 1
    ggsave(sprintf("demo_%03d.png", iter), p)
}
