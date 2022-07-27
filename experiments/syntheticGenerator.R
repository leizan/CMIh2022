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

test10 = function(n,m=3){
  z_1 = sample(0:(m-1),n,replace=T)
  z_2 = rbinom(n, p=0.5, size=3)
  z_3 = rexp(n, rate=1) 
  z_4 = rexp(n, rate=10) 
  x_1 = rep(0,n)
  x_2 = rep(0,n)
  x_3 = rep(0,n)
  y = rep(0,n)
  covM = matrix(c(1,0,0,1), nrow=2, ncol=2)
  for(i in 1:n){
    x1x2 = mvrnorm(n=1, mu=c(z_3[i],z_4[i]), Sigma=covM)
    x_1[i] = x1x2[1]
    x_2[i] = x1x2[2]
    x_3[i] = rbinom(1, p=0.5, size=z_1[i]+z_2[i]) 
    y[i] = rbinom(1, p=0.5, size=z_1[i]+z_2[i])
  }
  return(data.frame(x_1,x_2,x_3,y,z_1,z_2,z_3,z_4))
}

test11 = function(n,tCov=0.6,m=3,p=0.15){
  covM = matrix(c(1,tCov,tCov,1), nrow=2, ncol=2)
  xy = mvrnorm(n=n, Sigma=covM, mu=c(0,0))
  x_1 = xy[,1]
  y_1 = xy[,2]
  x_2 = sample(0:(m-1),n,replace=T)
  y_2 = rep(0,n)
  for(i in 1:n){
    y_2[i] = runif(1, min=x_2[i], max=x_2[i]+2)
  }
  x_3 = rexp(n, rate=1)
  y_3 = rep(0,n)
  for(i in 1:n){
    coin = runif(1, min=0, max=1)
    if(coin < p){
      y_3[i] = 0
    }else{
      y_3[i] = rpois(1, lambda=x_3[i])
    }
  }
  return(data.frame(x_1,x_2,x_3,y_1,y_2,y_3))
}
