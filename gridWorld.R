rm(list=ls(all=TRUE))

library(ggplot2)
library(grid)


newGrid<-function(start = c(6,1),step = -0.1,puddle = -10,goal = 50,
                  noise = 0.8){
  g = list()
  g$x = g$y = 10

  g$stepR = step
  g$puddleR = puddle
  g$goalR = goal

  if(start[1] < 1 || start[1] > g$x)
    stop("start[1] out of bounds")
  else if(start[2] < 1 || start[2] > g$y)
    stop("start[2] out of bounds")

  if(noise < 0.0 || noise > 1.0)
    stop("noise is not a probability")

  g$goal = c(8,9)

  g$r = matrix(step,g$x,g$y)

  g$r[g$goal[1],g$goal[2]] = goal
  g$r[g$goal[1]+1,g$goal[2]] = puddle
  g$r[g$goal[1]+1,g$goal[2]-1] = puddle
  g$r[g$goal[1],g$goal[2]-1] = puddle
  g$r[g$goal[1]-1,g$goal[2]-1] = puddle
  g$r[g$goal[1]-1,g$goal[2]] = puddle
  g$r[g$goal[1]-1,g$goal[2]+1] = puddle

  ind = expand.grid(3:8,2:4)
  g$r[ind[,1],ind[,2]] = puddle

  g$r[6,c(3,4)] = step

  ind = expand.grid(3:4,6:8)
  g$r[ind[,1],ind[,2]] = puddle

  g$r[7,c(6,7)] = puddle
  g$r[c(9,10),6] = puddle

  g$r[9,c(4,5)] = puddle

  g$r[6,6] = puddle

  g$s = start

  g$noise = noise

  g$actions = list(c(0,1),c(1,0),c(0,-1),c(-1,0))

  return(g)
}


transProb<-function(s,a,sp,g){
  if(s[1] < 1 || s[1] > g$x)
    stop("s[1] is outside of bounds")
  else if(s[2] < 1 || s[2] > g$y)
    stop("s[2] is outside of bounds")
  else if(sp[1] < 1 || sp[1] > g$x)
    stop("sp[1] is outside of bounds")
  else if(sp[2] < 1 || sp[2] > g$y)
    stop("sp[2] is outside of bounds")


  if(a < 1 || a > length(g$actions))
    stop("a is not a valid action")

  if(all(s == g$goal))
    return(ifelse(all(s == sp),1.0,0.0))

  sa = s + g$actions[[a]]
  sa[1] = min(max(sa[1],1),g$x)
  sa[2] = min(max(sa[2],1),g$y)

  prob = 0.0
  if(all(sa == sp))
    prob = prob + 1.0 - g$noise

  for(i in 1:length(g$actions)){
    sa = s + g$actions[[i]]
    sa[1] = min(max(sa[1],1),g$x)
    sa[2] = min(max(sa[2],1),g$y)

    if(all(sa == sp))
      prob = prob + g$noise/length(g$actions)
  }
  return(prob)
}



expReward<-function(s,a,g){
  if(s[1] < 1 || s[1] > g$x)
    stop("s[1] is outside of bounds")
  else if(s[2] < 1 || s[2] > g$y)
    stop("s[2] is outside of bounds")

  if(a < 1 || a > length(g$actions))
    stop("a is not a valid action")

  if(all(s == g$goal))
    return(0.0)

  sa = s + g$actions[[a]]
  sa[1] = min(max(sa[1],1),g$x)
  sa[2] = min(max(sa[2],1),g$y)

  er = g$r[sa[1],sa[2]] * (1.0 - g$noise)

  for(i in 1:length(g$actions)){
    sa = s + g$actions[[i]]
    sa[1] = min(max(sa[1],1),g$x)
    sa[2] = min(max(sa[2],1),g$y)

    er = er + g$r[sa[1],sa[2]] * g$noise / length(g$actions)
  }
  return(er)
}



nextStep<-function(a,g){
  if(a < 1 || a > length(g$actions))
    stop("a is not a valid action")



  tuple = list()
  tuple$s = g$s
  tuple$a = a

  if(!all(g$s == g$goal)){
    if(runif(1) < g$noise){
      print("noisy")
      a = sample(1:length(g$actions),1)
    }

    g$s = g$s + g$actions[[a]]
    g$s[1] = min(max(g$s[1],1),g$x)
    g$s[2] = min(max(g$s[2],1),g$y)


    tuple$r = g$r[g$s[1],g$s[2]]
    tuple$sp = g$s
  }
  else{
    tuple$r = 0
    tuple$sp = g$s
  }

  tuple$g = g

  return(tuple)
}



plotGrid<-function(g,policy=NULL,showState=FALSE){
  pData = NULL
  ind = 1
  for(x in 1:g$x){
    for(y in 1:g$y){
      pData = rbind(pData,c(x-0.5,y-0.5,g$r[x,y],ind,1))
      pData = rbind(pData,c(x+0.5,y-0.5,g$r[x,y],ind,2))
      pData = rbind(pData,c(x+0.5,y+0.5,g$r[x,y],ind,3))
      pData = rbind(pData,c(x-0.5,y+0.5,g$r[x,y],ind,4))
      ind = ind + 1
    }
  }
  pData = as.data.frame(pData)
  names(pData) = c("x","y","r","group","by")
  pData$r = factor(pData$r,levels=c(g$goalR,g$stepR,g$puddleR),
                   labels=c(paste("goal(",g$goalR,")",sep=""),
                            paste("step(",g$stepR,")",sep=""),
                            paste("puddle(",g$puddleR,")",sep="")))


  p = ggplot()
  p = p + geom_polygon(data = pData, aes(x=x,y=y,group=group,
                                         by=by,fill=r),
                       color = "black",size=0.25)

  p = p + scale_fill_manual("Rewards",
                            values=c("firebrick2",
                                     "chartreuse3",
                                     "dodgerblue3"))

  if(!is.null(policy)){
    if(nrow(policy) != g$x || ncol(policy) != g$y)
      stop("dimension of policy is invalid")
    else if(any(policy < 1 || policy > length(g$actions)))
      stop("policy exceeds bounds")

    pData = NULL
    for(x in 1:g$x){
      for(y in 1:g$y){
        if(!(x == g$goal[1] && y == g$goal[2])){
          a = g$actions[[policy[x,y]]] * 0.35

          pData = rbind(pData,c(x-a[1],y-a[2],x+a[1],y+a[2]))
        }
      }
    }
    pData = as.data.frame(pData)
    names(pData) = c("x","y","xend","yend")

    p = p + geom_segment(data = pData, aes(x=x,y=y,xend=xend,yend=yend),
                         arrow=arrow(length=unit(0.1,"in"),angle=50))
  }


  if(showState){
    sPoint = data.frame(x = g$s[1], y = g$s[2])
    p = p + geom_point(data=sPoint,aes(x=x,y=y),color="black",
                       size=6,shape=21,fill="gold1")
  }

  p = p + coord_fixed()

  p = p + theme(panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title= element_blank())

  p = p + ggtitle("Grid World")

  print(p)
}
