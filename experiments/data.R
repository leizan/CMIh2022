library(numbers)
library(hesim)

indexLabel = c('0'='A', '1'='B', '2'='C', '3'='D', '4'='E')
labelIndex = c('A'='0', 'B'='1', 'C'='2', 'D'='3', 'E'='4')

fMap = function(number){
  mapping = c('C', 'A', 'D', 'B')
  return(mapping[number+1])
}

gMap = function(number){
  mapping = c('B', 'A', 'C')
  return(mapping[number+1])
}

pMap = function(number){
  mapping = c('C', 'A', 'D', 'E', 'B')
  return(mapping[number+1])
}

qMap = function(number){
  mapping = c('D', 'E', 'B', 'A', 'C')
  return(mapping[number+1])
}

hMap = function(numberX, numberY){
  valueX = as.numeric(labelIndex[numberX])
  valueY = as.numeric(labelIndex[numberY])
  return(5*valueX+valueY)
}

mMap = function(numberX, numberY){
  mapping = c('AA'='0', 'AB'='0', 'AC'='0',
              'AD'='1', 'AE'='1', 'BA'='1', 'BB'='1', 'BC'='1',
              'BD'='2', 'BE'='2', 'CA'='2', 'CB'='2', 'CC'='2',
              'CD'='3', 'CE'='3', 'DA'='3', 'DB'='3', 'DC'='3',
              'DD'='4', 'DE'='4', 'EA'='4', 'EB'='4', 'EC'='4',
              'ED'='0', 'EE'='0')
  XY = paste(numberX, numberY, sep="")
  return(as.numeric(mapping[XY]))
}

nMap = function(numberX){
  mapping = c('A'='1', 'B'='3', 'C'='0', 'D'='2', 'E'='4')
  return(as.numeric(mapping[numberX]))
}

rMap = function(numberY){
  mapping = c('B','A','C','D')
  return(mapping[numberY+1])
}

sMap = function(numberX){
  mapping = c('A'='1', 'B'='4', 'C'='0', 'D'='2', 'E'='3')
  return(as.numeric(mapping[numberX]))
}

tMap = function(numberY){
  mapping = c('D','E','B','A','C')
  return(mapping[numberY+1])
}


Chain_a1 = function(n){
  x = runif(n, min=0, max=100)
  z = x + rnorm(n)
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rbinom(1, abs(round(z[i])), 0.5)
  }
  data = data.frame(x,y,z)
  return(data)
}

Chain_a2 = function(n){
  x = runif(n, min=0, max=100)
  z = x + rnorm(n)
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rpois(1, abs(round(z[i])))
  }
  data = data.frame(x,y,z)
  return(data)
}

Chain_b = function(n){
  x = sample(c('A','B','C','D','E'), n, replace=T, prob=c(0.6,0.1,0.1,0.1,0.1))
  z = rep(0,n)
  y = rep(0,n)
  for(i in 1:n){
    z[i] = rnorm(1,mean=nMap(x[i]),sd=2)
    y[i] = rMap(mod(round(z[i]+rnorm(1)),4))
  }
  x = as.numeric(labelIndex[x])
  y = as.numeric(labelIndex[y])
  data = data.frame(x,y,z)
  return(data)
}

Chain_c1 = function(n){
  x = runif(n, min=0, max=100)
  z = rep(0,n)
  for(i in 1:n){
    z[i] = sample(round(x[i]):(round(x[i])+2),1,replace=T)
  }
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rbinom(1,z[i],0.5)
  }
  data = data.frame(x,y,z)
  return(data)
}

Chain_c2 = function(n){
  x = runif(n, min=0, max=100)
  z = rep(0,n)
  for(i in 1:n){
    z[i] = sample(round(x[i]):(round(x[i])+2),1,replace=T)
  }
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rpois(1, z[i])
  }
  data = data.frame(x,y,z)
  return(data)
}

Chain_d = function(n){
  x = sample(c('A','B','C','D','E'), n, replace=T,prob=c(0.6,0.1,0.1,0.1,0.1))
  z = rep(0,n)
  y = rep(0,n)
  for(i in 1:n){
    probZ = rep(0.025,5)
    probY = rep(0.025,5)
    probZ[sMap(x[i])+1]=0.9
    z[i] = rcat(1,probZ) -1
    probY[as.numeric(labelIndex[tMap(z[i])])+1]=0.9
    y[i] = rcat(1,probY) -1
  }
  x = as.numeric(labelIndex[x])
  data = data.frame(x,y,z)
  return(data)
}

Fork_a1 = function(n){
  z = runif(n, min=0, max=100)
  x = z + rnorm(n) 
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rbinom(1,round(z[i]),0.5)
  }
  data = data.frame(x,y,z)
  return(data)
}

Fork_a2 = function(n){
  z = runif(n, min=0, max=100)
  x = z + rnorm(n) 
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rpois(1, round(z[i]))
  }
  data = data.frame(x,y,z)
  return(data)
}

Fork_b = function(n){
  z = rexp(n, rate=0.1)
  x = sapply(mod(round(z+rnorm(n)),4),fMap)
  y = sapply(mod(round(z+rnorm(n)),3),gMap)
  x = as.numeric(labelIndex[x])
  y = as.numeric(labelIndex[y])
  data = data.frame(x,y,z)
  return(data)
}

Fork_c1 = function(n){
  z = sample(0:100,n,replace=T)
  x = z + rnorm(n)
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rbinom(1,z[i],0.5)
  }
  data = data.frame(x,y,z)
  return(data)
}

Fork_c2 = function(n){
  z = sample(0:100,n,replace=T)
  x = z + rnorm(n)
  y = rep(0,n)
  for(i in 1:n){
    y[i] = rpois(1, z[i])
  }
  data = data.frame(x,y,z)
  return(data)
}

Fork_d = function(n){
  z = sample(c('A','B','C','D','E'), n, replace=T,prob=c(0.6,0.1,0.1,0.1,0.1))
  z = as.numeric(labelIndex[z])
  x = rep(0,n)
  y = rep(0,n)
  for(i in 1:n){
    probX = rep(0.025,5)
    probY = rep(0.025,5)
    probX[as.numeric(labelIndex[pMap(z[i])])+1]=0.9
    probY[as.numeric(labelIndex[qMap(z[i])])+1]=0.9
    x[i] = rcat(1,probX) -1
    y[i] = rcat(1,probY) -1
  }
  data = data.frame(x,y,z)
  return(data)
}

Collider_a1 = function(n){
  x = rnorm(n, mean=50, sd=25)
  y = rbinom(n, 100, 0.5)
  z = x + y + rnorm(n)
  data = data.frame(x,y,z)
  return(data)
}

Collider_a2 = function(n){
  x = rnorm(n, mean=50, sd=25)
  y = rpois(n, 100)
  z = x + y + rnorm(n)
  data = data.frame(x,y,z)
  return(data)
}

Collider_b = function(n){
  x = sample(c('A','B','C','D','E'), n, replace=T,prob=c(0.6,0.1,0.1,0.1,0.1))
  y = sample(c('A','B','C','D','E'), n, replace=T,prob=c(0.6,0.1,0.1,0.1,0.1))
  z = rep(0,n)
  for(i in 1:n){
    z[i] = rnorm(1,mean=hMap(x[i],y[i]),sd=1)
  }
  x = as.numeric(labelIndex[x])
  y = as.numeric(labelIndex[y])
  data = data.frame(x,y,z)
  return(data)
}

Collider_c1 = function(n){
  x = rnorm(n, mean=50, sd=25) # sd=50 
  y = rbinom(n, 100, 0.5)
  z = rep(0,n)
  eleList = abs(round(x + y + rnorm(n)))
  for(i in 1:n){
    z[i] = rbinom(1, eleList[i], 0.5)
  }
  data = data.frame(x,y,z)
  return(data)
}

Collider_c2 = function(n){
  x = rnorm(n, mean=50, sd=25) # sd=50
  y = rpois(n, 100)
  z = rep(0,n)
  eleList = abs(round(x + y + rnorm(n)))
  for(i in 1:n){
    z[i] = rbinom(1, eleList[i], 0.5)
  }
  data = data.frame(x,y,z)
  return(data)
}

Collider_d = function(n){
  x = sample(c('A','B','C','D','E'), n, replace=T,prob=c(0.6,0.1,0.1,0.1,0.1))
  y = sample(c('A','B','C','D','E'), n, replace=T,prob=c(0.6,0.1,0.1,0.1,0.1))
  z = rep(0,n)
  for(i in 1:n){
    probZ = rep(0.025,5)
    probZ[mMap(x[i],y[i])+1]=0.9
    z[i] = rcat(1,probZ) -1
  }
  x = as.numeric(labelIndex[x])
  y = as.numeric(labelIndex[y])
  data = data.frame(x,y,z)
  return(data)
}