
# The input data is purely continuous and well tailored, the output is the relative distance array
getDistArray = function(data){
  N = dim(data)[1]
  nDim = dim(data)[2]
  inds = 1:nDim
  disArray = array(rep(NaN, N*N*nDim), dim =c(N,N,nDim))
  for(m in inds){
    dataDim = as.matrix(data[,m])
    disArray[,,m] = as.matrix(parDist(dataDim, method="manhattan"))
  }
  return(disArray)
}

#Get epsilon for each point 
getEpsilonDistance = function(k, disArray){
  N = dim(disArray)[1]
  epsilonDisArray = rep(0,N)
  for(point_i in 1:N){
    coord_dists =  as.matrix(disArray[,point_i,])
    # print(coord_dists)
    dists = as.numeric(rowMax(coord_dists, ignore.zero=F))
    ordered_dists = sort(dists)
    epsilonDisArray[point_i]= 2*ordered_dists[k+1]
  }
  return(epsilonDisArray)
}

#Find the interset among serveral sets 
findInterCluster = function(eleInEachClass){
  interCluster = eleInEachClass[[1]]
  for(m in 2:length(eleInEachClass)){
    interCluster = intersect(interCluster, eleInEachClass[[m]])
  }
  return(interCluster)
}

#Estimator of entropy for purely continuous data
conEntroEstimator = function(data,k,dN){
  N = dim(data)[1]
  if(N==1){
    return(0)
  }
  distArray = getDistArray(data)
  epsilonDis = getEpsilonDistance(k, distArray)
  if(0 %in% epsilonDis){
    epsilonDis = epsilonDis[epsilonDis!=0]
    N = length(epsilonDis)
    if(N==0){
      return(0)
    }
    entropy = -digamma(k) + digamma(N) + (dN*sum(log(epsilonDis)))/N
    return(entropy)
  }
  # The maximum norm is used, so Cd=1, log(Cd)=0
  entropy = -digamma(k) + digamma(N) + (dN*sum(log(epsilonDis)))/N
  return(entropy)
}


#Calculate the information entropy for mixed data
mixedEntroEstimator=function(data, dimCon, dimDis){
  # Input data should be as matrix
  dN = length(dimCon)
  N = dim(data)[1]
  entroCon = 0
  entroDis = 0
  if(length(dimCon)!=0){
    dataCon = as.matrix(data[,dimCon])
  }
  if(length(dimDis)!=0){
    dataDis = as.matrix(data[,dimDis])
  }
  if(length(dimCon)==0 && length(dimDis)==0){
    # print('Input data is NULL!!!')
  }
  # If the data is purely continuous!
  if(length(dimDis)==0 && length(dimCon)!=0){
    entroCon = conEntroEstimator(dataCon, max(1,round(0.1*N)), dN)
  }
  if(length(dimDis)!=0){
    #classByDimList = rep(c(), length(dimDis))
    classByDimList = vector(mode="list", length=length(dimDis))
    for(i in 1:length(dimDis)){
      classByDimList[[i]]=unique(dataDis[,i])
    }
    classList = as.matrix(expand.grid(classByDimList))
    #indexInClass = rep(c(), dim(classList)[1])
    indexInClass = vector(mode="list", length=dim(classList)[1])
    for(i in 1:dim(classList)[1]){
      classElement = as.numeric(classList[i,])
      if(length(classElement)==1){
        indexInClass[[i]] = which(dataDis==classElement)
      }else{
        #eleInEachClass=rep(c(),length(classElement))
        eleInEachClass=vector(mode="list", length=length(classElement))
        for(m in 1:length(classElement)){
          eleInEachClass[[m]] = which(dataDis[,m]==classElement[m])
        }
        indexInClass[[i]] = findInterCluster(eleInEachClass)
      }
    }
    #Remove the empty bins 
    for(i in length(indexInClass):1){
      if(length(indexInClass[[i]])==0){
        indexInClass[[i]] = NULL
      }
    }
    probBins=rep(0,length(indexInClass))
    for(i in 1:length(indexInClass)){
      probBins[i] = length(indexInClass[[i]])/N
    }
    for(i in probBins){
      entroDis = entroDis - i*log(i)
    }
  }
  if(length(dimDis)!=0 && length(dimCon)!=0){
    for(i in 1:length(probBins)){
      entroCon = entroCon + conEntroEstimator(as.matrix(dataCon[indexInClass[[i]],]), max(1,round(0.1*length(indexInClass[[i]]))), dN)*probBins[i]
    }
  }
  res = entroCon + entroDis    
  return(res)
}

mixedEstimator = function(data, xind, yind, zinds, isCat){
  xDimCon = setdiff(xind,isCat)
  xDimDis = setdiff(xind, xDimCon)
  yDimCon = setdiff(yind,isCat)
  yDimDis = setdiff(yind, yDimCon)
  zDimCon = setdiff(zinds,isCat)
  zDimDis = setdiff(zinds, zDimCon)
  
  conXYZ = c(xDimCon,yDimCon,zDimCon)
  disXYZ = c(xDimDis,yDimDis,zDimDis)
  hXYZ = mixedEntroEstimator(data, conXYZ, disXYZ)
  conXZ = c(xDimCon,zDimCon)
  disXZ = c(xDimDis,zDimDis)
  hXZ = mixedEntroEstimator(data, conXZ, disXZ)
  conYZ = c(yDimCon,zDimCon)
  disYZ = c(yDimDis,zDimDis)
  hYZ = mixedEntroEstimator(data, conYZ, disYZ)
  conZ = zDimCon
  disZ = zDimDis
  hZ = mixedEntroEstimator(data, conZ, disZ)
  cmi = hXZ + hYZ - hXYZ - hZ
  return(cmi)
}

