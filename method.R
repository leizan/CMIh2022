library(parallelDist)
library(qlcMatrix)
library(doParallel)
library(foreach)

# MS_local
localPermuTestMS_local = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  
  # rang transformation applied to no nominal dimension
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = KNN.CMI_estimates_MS(data=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, k=k)
  
  neighbourList = list()
  
  zz = as.matrix(dd[,zinds])
  zCatAbs = isCat[which((isCat-2)>=1)]
  zCat = rep(0, length(zCatAbs))
  indexCat = 1 
  for(i in zCatAbs){
    zCat[indexCat] = which(zinds==i)
    indexCat = indexCat + 1 
  }
 
  ndimZ = length(zinds)
  
  for(p in 1:total){
    disArray = zz
    pointZ = zz[p,]
    for(d in 1:ndimZ){
      if(d %in% zCat){
        for(n in 1:total){
          if(disArray[n,][[d]]== pointZ[[d]]){
            disArray[n,][[d]] = 0
          }else{
            disArray[n,][[d]] = Inf
          }
        }
      }else{
        disArray[,d] = abs(disArray[,d]-pointZ[[d]])
      }
    }
    disFinal = as.numeric(rowMax(disArray, ignore.zero=F))
    disKPer = sort(disFinal)[kPer]
    neighbourList[[p]] = which(disFinal<=disKPer)
  }  
  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("MS-R.R")
    usedIndice = c()
    xNew = rep(0,total)
    for(p in 1:total){
      # Shuffle neighbourList for each p 
      neighbourList[[p]] = sample(neighbourList[[p]])
    }
    # Create random permutation
    pi=sample(1:total)
    for(n in pi){
      j = neighbourList[[n]][1]
      nPer = length(neighbourList[[n]])
      m = 1
      while((j %in% usedIndice)&(m<nPer)){
        m=m+1
        j = neighbourList[[n]][m]
      }
      xNew[n]= dd[j,][[1]]
      usedIndice=append(usedIndice, j)
    }
    data = dd
    data[[1]] = xNew
    KNN.CMI_estimates_MS(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat, k=k)
    
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}

# CMIh_local
localPermuTestCMIh_local = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  # rang transformation applied to no nominal dimension
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = mixedEstimator(data=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat)

  neighbourList = list()
  zz = as.matrix(dd[,zinds])
  zCatAbs = isCat[which((isCat-2)>=1)]
  zCat = rep(0, length(zCatAbs))
  indexCat = 1 
  for(i in zCatAbs){
    zCat[indexCat] = which(zinds==i)
    indexCat = indexCat + 1 
  }
  
  ndimZ = length(zinds)
  
  for(p in 1:total){
    disArray = zz
    pointZ = zz[p,]
    for(d in 1:ndimZ){
      if(d %in% zCat){
        for(n in 1:total){
          if(disArray[n,][[d]]== pointZ[[d]]){
            disArray[n,][[d]] = 0
          }else{
            disArray[n,][[d]] = Inf
          }
        }
      }else{
        disArray[,d] = abs(disArray[,d]-pointZ[[d]])
      }
    }
    disFinal = as.numeric(rowMax(disArray, ignore.zero=F))
    disKPer = sort(disFinal)[kPer]
    neighbourList[[p]] = which(disFinal<=disKPer)
  }  
  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("mixedCmiIEstimator.R")
    usedIndice = c()
    xNew = rep(0,total)
    for(p in 1:total){
      # Shuffle neighbourList for each p 
      neighbourList[[p]] = sample(neighbourList[[p]])
    }
    # Create random permutation
    pi=sample(1:total)
    for(n in pi){
      j = neighbourList[[n]][1]
      nPer = length(neighbourList[[n]])
      m = 1
      while((j %in% usedIndice)&(m<nPer)){
        m=m+1
        j = neighbourList[[n]][m]
      }
      xNew[n]= dd[j,][[1]]
      usedIndice=append(usedIndice, j)
    }
    data = dd
    data[[1]] = xNew
    mixedEstimator(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
    
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}

#############################################################################
#############################################################################

#MS_local_const
localPermuTestMS_local_const = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  # rang transformation applied to no nominal dimension
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = KNN.CMI_estimates_MS(data=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, k=k)
  
  # k class min of Z to update k for local permutation
  zDis = intersect(zinds, isCat)
  
  if(length(zDis)!=0){
    zEleTab = table(dd[,zDis])
  }
  
  neighbourList = list()
  zz = as.matrix(dd[,zinds])
  zCatAbs = isCat[which((isCat-2)>=1)]
  zCat = rep(0, length(zCatAbs))
  indexCat = 1 
  for(i in zCatAbs){
    zCat[indexCat] = which(zinds==i)
    indexCat = indexCat + 1 
  }
  
  ndimZ = length(zinds)
  
  for(p in 1:total){
    disArray = zz
    pointZ = zz[p,]
    for(d in 1:ndimZ){
      if(d %in% zCat){
        for(n in 1:total){
          if(disArray[n,][[d]]== pointZ[[d]]){
            disArray[n,][[d]] = 0
          }else{
            disArray[n,][[d]] = Inf
          }
        }
      }else{
        disArray[,d] = abs(disArray[,d]-pointZ[[d]])
      }
    } 
    if(length(zCat)!=0){
      eleNumb = do.call(`[`, c(list(zEleTab), as.character(zz[p,zCat])))
      kPer = min(eleNumb,kPer)
    }
    disFinal = as.numeric(rowMax(disArray, ignore.zero=F))
    disKPer = sort(disFinal)[kPer]
    neighbourList[[p]] = which(disFinal<=disKPer)
  }  
  
  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("MS-R.R")
    usedIndice = c()
    newIndice = rep(0,total)
    for(p in 1:total){
      # Shuffle neighbourList for each p 
      neighbourList[[p]] = sample(neighbourList[[p]])
    }
    # Create random permutation
    pi=sample(1:total)
    for(n in pi){
      j = neighbourList[[n]][1]
      nPer = length(neighbourList[[n]])
      m = 1
      while((j %in% usedIndice)&(m<nPer)){
        m=m+1
        j = neighbourList[[n]][m]
      }
      newIndice[n]= dd[j,][[1]]
      usedIndice=append(usedIndice, j)
    }
    
    data = dd
    data[[1]] = newIndice
    
    KNN.CMI_estimates_MS(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat, k=k)
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}

#CMIh_local_const
localPermuTestCMIh_local_const = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  # rang transformation applied to no nominal dimension
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = mixedEstimator(data=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  
  # k class min of Z to update k for local permutation
  zDis = intersect(zinds, isCat)
  
  if(length(zDis)!=0){
    zEleTab = table(dd[,zDis])
  }
  
  neighbourList = list()
  zz = as.matrix(dd[,zinds])
  zCatAbs = isCat[which((isCat-2)>=1)]
  zCat = rep(0, length(zCatAbs))
  indexCat = 1 
  for(i in zCatAbs){
    zCat[indexCat] = which(zinds==i)
    indexCat = indexCat + 1 
  }
  
  ndimZ = length(zinds)
  
  for(p in 1:total){
    disArray = zz
    pointZ = zz[p,]
    for(d in 1:ndimZ){
      if(d %in% zCat){
        for(n in 1:total){
          if(disArray[n,][[d]]== pointZ[[d]]){
            disArray[n,][[d]] = 0
          }else{
            disArray[n,][[d]] = Inf
          }
        }
      }else{
        disArray[,d] = abs(disArray[,d]-pointZ[[d]])
      }
    } 
    if(length(zCat)!=0){
      eleNumb = do.call(`[`, c(list(zEleTab), as.character(zz[p,zCat])))
      kPer = min(eleNumb,kPer)
    }
    disFinal = as.numeric(rowMax(disArray, ignore.zero=F))
    disKPer = sort(disFinal)[kPer]
    neighbourList[[p]] = which(disFinal<=disKPer)
  }  
  
  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("mixedCmiIEstimator.R")
    usedIndice = c()
    newIndice = rep(0,total)
    for(p in 1:total){
      # Shuffle neighbourList for each p 
      neighbourList[[p]] = sample(neighbourList[[p]])
    }
    # Create random permutation
    pi=sample(1:total)
    for(n in pi){
      j = neighbourList[[n]][1]
      nPer = length(neighbourList[[n]])
      m = 1
      while((j %in% usedIndice)&(m<nPer)){
        m=m+1
        j = neighbourList[[n]][m]
      }
      newIndice[n]= dd[j,][[1]]
      usedIndice=append(usedIndice, j)
    }
    
    data = dd
    data[[1]] = newIndice
    
    mixedEstimator(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}

###############################################################################
###############################################################################

#MS 
globalPermuTest_MS = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = KNN.CMI_estimates_MS(data=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, k=k)

  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("MS-R.R")
    data = dd
    data[[1]] = sample(dd[[1]])
    KNN.CMI_estimates_MS(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat, k=k)
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}

# MS without condition
globalPermuTest_MS_without_condition = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = KNN.MI_estimates_MS(data=dd, xind=xind, yind=yind, isCat=isCat, k=k)
  
  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("MS-R.R")
    data = dd
    data[[1]] = sample(dd[[1]])
    KNN.MI_estimates_MS(data=data, xind=xind, yind=yind, isCat=isCat, k=k)
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}

#CMIh
globalPermuTest_CMIh = function(dd, xind, yind, zinds, isCat, B=2, kPer=7, k=k){
  
  total = dim(dd)[1]
  dimNum = dim(dd)[2]
  for(d in 1:dimNum){
    if(!(d %in% isCat)){
      dd[[d]] = rank(dd[[d]],ties.method = c("first"))
    }
  }
  
  I_original = mixedEstimator(data=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  
  estimatesKnn = foreach(l = 1:B, .combine=rbind, .packages = c('parallelDist','qlcMatrix','dplyr','MASS')) %dopar% {
    source("mixedCmiIEstimator.R")
    data = dd
    data[[1]] = sample(dd[[1]])
    mixedEstimator(data=data, xind=xind, yind=yind, zinds=zinds, isCat=isCat)
  }
  pValue = length(which(estimatesKnn[,1] >= I_original))/B
  return(pValue)
}
