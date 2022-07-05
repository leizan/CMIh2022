test2 = function(n,m=3){
  x = sample(0:(m-1),n,replace=T)
  y = rep(0,n)
  for(i in 1:n){
    y[i] = runif(1, min=x[i], max=x[i]+2)
  }
  return(data.frame(x,y))
}

test3 = function(n,p=0.15){
  x = rexp(n, rate=1)
  y = rep(0,n)
  for(i in 1:n){
    coin = runif(1, min=0, max=1)
    if(coin < p){
      y[i] = 0
    }else{
      y[i] = rpois(1, lambda=x[i])
    }
  }
  return(data.frame(x,y))
}

test5 = function(n,m=3){
  x = sample(0:(m-1),n,replace=T)
  y = rep(0,n)
  z = rbinom(n, p=0.5, size=3)
  for(i in 1:n){
    y[i] = runif(1, min=x[i], max=x[i]+2)
  }
  return(data.frame(x,y,z))
}

test6 = function(n,p=0.15){
  x = rexp(n, rate=1)
  y = rep(0,n)
  z = rbinom(n, p=0.5, size=3)
  for(i in 1:n){
    coin = runif(1, min=0, max=1)
    if(coin < p){
      y[i] = 0
    }else{
      y[i] = rpois(1, lambda=x[i])
    }
  }
  return(data.frame(x,y,z))
}

test7 = function(n){
  z = rbinom(n, p=0.5, size=9)
  x = rep(0,n)
  y = rep(0,n)
  for(i in 1:n){
    x[i] = rnorm(1, mean=z[i], sd=1)
    y[i] = rnorm(1, mean=z[i], sd=1)
  }
  return(data.frame(x,y,z))
}

test8 = function(n){
  x = sample(0:4,n,replace=T) 
  y = rep(0,n)
  z = rep(0,n)
  for(i in 1:n){
    y[i] = runif(1, min=x[i], max=x[i]+2)
    z[i] = rbinom(1, p=0.5, size=x[i])
  }
  return(data.frame(x,y,z))
}

test9 = function(n){
  # X -> Z -> Y  ----> I(X,Y | Z) = 0
  x = rexp(n, rate=10) 
  z = rep(0,n)
  y = rep(0,n)
  for(i in 1:n){
    z[i] = rpois(1, lambda=x[i])
    y[i] = rbinom(1, p=0.5, size=z[i]+5) 
  }
  return(data.frame(x,y,z))
}
