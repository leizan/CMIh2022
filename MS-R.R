# Mise a jour 
getPairwiseDistArray_IO = function(data, inds, isCat=c()){
  N = dim(data)[1]
  nDim = dim(data)[2]
  if(length(inds)==0){
    inds = 1:nDim
  }
  disArray = array(rep(NaN, N*N*nDim), dim =c(N,N,nDim))
  for(m in inds){
    dataDim = as.matrix(data[,m])
    if(m %in% isCat){
      for(i in 1:N){
        for(j in 1:N){
          if(dataDim[i]==dataDim[j]){
            disArray[i,j,m] = 0
          }else{
            disArray[i,j,m] = 1
          }
        }
      }
    }else{
      dataDim = as.matrix(dataDim, nrow=N)
      disArray[,,m] = as.matrix(parDist(dataDim, method="manhattan"))
    }
  }
  return(disArray)
}

getPairwiseDistArray_II = function(data, inds, isCat=c()){
  N = dim(data)[1]
  nDim = dim(data)[2]
  if(length(inds)==0){
    inds = 1:nDim
  }
  disArray = array(rep(NaN, N*N*nDim), dim =c(N,N,nDim))
  for(m in inds){
    dataDim = as.matrix(data[,m])
    if(m %in% isCat){
      for(i in 1:N){
        for(j in 1:N){
          if(dataDim[i]==dataDim[j]){
            disArray[i,j,m] = 0
          }else{
            disArray[i,j,m] = Inf
          }
        }
      }
    }else{
      dataDim = as.matrix(dataDim, nrow=N)
      disArray[,,m] = as.matrix(parDist(dataDim, method="manhattan"))
    }
  }
  return(disArray)
}

getPointCoordDists = function(distArray, point_i, inds=c()){
  if(length(inds)==0){
    inds = 1:dim(distArray)[3]
  }
  obsDists =  distArray[,point_i, inds]
  return(obsDists)
}

getKnnDist = function(coord_dists, k){
  dists = as.numeric(rowMax(coord_dists, ignore.zero=F))
  ordered_dists = sort(dists)
  k_tilde = length(which(dists <= ordered_dists[k+1]))-1
  return(c(k_tilde, ordered_dists[k+1]))
}

getKnnDistCont = function(coord_dists, k){
  dists = as.numeric(rowMax(coord_dists, ignore.zero=F))
  ordered_dists = sort(dists)
  return(c(k, ordered_dists[k+1]))
}

getKnnDistMin = function(coord_dists, k, isCat){
  if(length(isCat)!=0){
    distsCat = as.numeric(rowMax(coord_dists[,isCat], ignore.zero=F))
    k = length(which(distsCat==0))-1
  }
  dists = as.numeric(rowMax(coord_dists, ignore.zero=F))
  ordered_dists = sort(dists)
  k_tilde = max(length(which(dists <= ordered_dists[k+1]))-1,1)
  return(c(k_tilde, ordered_dists[k+1]))
}

getKnnDistRavk = function(coord_dists, k){
  dists = as.numeric(rowMax(coord_dists, ignore.zero=F))
  ordered_dists = sort(dists)
  if(ordered_dists[k+1]==0){
    k_tilde = length(which(dists == 0)) - 1
  }else{
    k_tilde = k
    
  }
  return(c(k_tilde, ordered_dists[k+1]))
}

countNeighborsCont = function(coord_dists, rho, coords=c()){
  # if(length(coords)==0){
  #   coords = 1:dim(coord_dists)[2]
  # }
  if(length(coords)>=1){
    
    dists = rowMax(coord_dists[,coords], ignore.zero=F)
  }else{
    return(1)
  }
  count = max(length(which(dists < rho)) - 1, 1)
  return(count)
}

countNeighbors= function(coord_dists, rho, coords=c()){
  if(length(coords)==0){
    coords = 1:dim(coord_dists)[2]
  }
  if(length(coords)>=1){
    dists = rowMax(coord_dists[,coords], ignore.zero=F)
  }else{
    return(1)
  }
  count = max(length(which(dists <= rho)) - 1, 1)
  return(count)
}

miEstMS = function(coord_dists, k, N, x_coords, y_coords){
  resKR = getKnnDist(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nx = countNeighbors(coord_dists, rho, x_coords)
  ny = countNeighbors(coord_dists, rho, y_coords)
  xiMS = digamma(k_tilde) - digamma(nx) - digamma(ny) + digamma(N)
  return(xiMS)
}

miEstFP = function(coord_dists, k, N, x_coords, y_coords){
  resKR =  getKnnDistCont(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nx = countNeighborsCont(coord_dists, rho, x_coords)
  ny = countNeighborsCont(coord_dists, rho, y_coords)
  xiFP = digamma(k_tilde) - digamma(nx) - digamma(ny) + digamma(N)
  return(xiFP)
}

miEstRAVK1 = function(coord_dists, k, N, x_coords, y_coords){
  resKR = getKnnDistRavk(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nx = countNeighbors(coord_dists, rho, x_coords)
  ny = countNeighbors(coord_dists, rho, y_coords)
  xiRAVK1 = digamma(k_tilde) - log(nx + 1) - log(ny + 1) + log(N+1)
  return(xiRAVK1)  
}

miEstRAVK2 = function(coord_dists, k, N, x_coords, y_coords){
  resKR = getKnnDist(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nx = countNeighbors(coord_dists, rho, x_coords)
  ny = countNeighbors(coord_dists, rho, y_coords)
  xiRAVK2 = digamma(k_tilde) - log(nx + 1) - log(ny + 1) + log(N+1)
  return(xiRAVK2)
}

cmiEstFP = function(coord_dists, k,  x_coords, y_coords, z_coords){
  resKR =  getKnnDistCont(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nxz = countNeighborsCont(coord_dists, rho, c(x_coords,z_coords))
  nyz = countNeighborsCont(coord_dists, rho, c(y_coords,z_coords))
  nz = countNeighborsCont(coord_dists, rho, z_coords)
  xiFP = digamma(k_tilde) - digamma(nxz) - digamma(nyz) + digamma(nz)
  return(xiFP)
}

cmiEstRAVK1 = function(coord_dists, k,  x_coords, y_coords, z_coords){
  resKR = getKnnDistRavk(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nxz = countNeighbors(coord_dists, rho, c(x_coords,z_coords))
  nyz = countNeighbors(coord_dists, rho, c(y_coords,z_coords))
  nz = countNeighbors(coord_dists, rho, z_coords)
  xiRAVK1 = digamma(k_tilde) - log(nxz + 1) - log(nyz + 1) + log(nz + 1)
  return(xiRAVK1)  
}

cmiEstRAVK2 = function(coord_dists, k,  x_coords, y_coords, z_coords){
  resKR = getKnnDist(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nxz = countNeighbors(coord_dists, rho, c(x_coords,z_coords))
  nyz = countNeighbors(coord_dists, rho, c(y_coords,z_coords))
  nz = countNeighbors(coord_dists, rho, z_coords)
  xiRAVK2 = digamma(k_tilde) - log(nxz + 1) - log(nyz + 1) + log(nz + 1)
  return(xiRAVK2)
}

cmiEstMS = function(coord_dists, k,  x_coords, y_coords, z_coords){
  resKR = getKnnDist(coord_dists=coord_dists, k=k)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nxz = countNeighbors(coord_dists, rho, c(x_coords,z_coords))
  nyz = countNeighbors(coord_dists, rho, c(y_coords,z_coords))
  nz = countNeighbors(coord_dists, rho, z_coords)
  xiMS = digamma(k_tilde) - digamma(nxz) - digamma(nyz) + digamma(nz)
  return(xiMS)
}

cmiEstMin = function(coord_dists, k,  x_coords, y_coords, z_coords, isCat){
  resKR = getKnnDistMin(coord_dists=coord_dists, k=k, isCat=isCat)
  k_tilde = as.integer(resKR[1])
  rho = resKR[2]
  nxz = countNeighbors(coord_dists, rho, c(x_coords,z_coords))
  nyz = countNeighbors(coord_dists, rho, c(y_coords,z_coords))
  nz = countNeighbors(coord_dists, rho, z_coords)
  xiMS = digamma(k_tilde) - digamma(nxz) - digamma(nyz) + digamma(nz)
  return(xiMS)
}

KNN.MI_estimates_FP = function(data, xind, yind, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  miFP = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    miFP = miFP + miEstFP(coord_dists=coord_dists, k=k,  N=N, x_coords=x_coords, y_coords=y_coords)
    
  }
  return(miFP/N)
}

KNN.MI_estimates_RAVK1 = function(data, xind, yind, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  miRAVK1 = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    miRAVK1 = miRAVK1 + miEstRAVK1(coord_dists=coord_dists, k=k, N=N, x_coords=x_coords, y_coords=y_coords)
  }
  return(miRAVK1/N)
}

KNN.MI_estimates_RAVK2 = function(data, xind, yind, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  miRAVK2 = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    miRAVK2 = miRAVK2 + miEstRAVK2(coord_dists=coord_dists, k=k, N=N, x_coords=x_coords, y_coords=y_coords)
    
  }
  return(miRAVK2/N)
}

KNN.MI_estimates_Parall = function(data, xind, yind, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  miFP = 0
  miRAVK2 = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
  
    miFP = miFP + miEstFP(coord_dists=coord_dists, k=k,  N=N, x_coords=x_coords, y_coords=y_coords)
    miRAVK2 = miRAVK2 + miEstRAVK2(coord_dists=coord_dists, k=k, N=N, x_coords=x_coords, y_coords=y_coords)
  }
  return(c(miFP/N, miRAVK2/N))
}

KNN.MI_estimates_MS = function(data, xind, yind, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  miMS = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  for(i in 1:N){
    MI_coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    miMS = miMS + miEstMS(coord_dists=MI_coord_dists, k=k,  N=N, x_coords=x_coords, y_coords=y_coords)
  }
  return(miMS/N)
} 

KNN.CMI_estimates_FP = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiFP = 0
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    x_coords = 1:length(xind)
    y_coords = (length(xind)+1):length(c(xind, yind))
    z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
    
    cmiFP = cmiFP + cmiEstFP(coord_dists=coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
    
  }
  return(cmiFP/N)
}

KNN.CMI_estimates_RAVK1 = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiRAVK1 = 0
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    x_coords = 1:length(xind)
    y_coords = (length(xind)+1):length(c(xind, yind))
    z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
    
    cmiRAVK1 = cmiRAVK1 + cmiEstRAVK1(coord_dists=coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
  }
  return(cmiRAVK1/N)
}

KNN.CMI_estimates_RAVK2 = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiRAVK2 = 0
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    x_coords = 1:length(xind)
    y_coords = (length(xind)+1):length(c(xind, yind))
    z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
    
    cmiRAVK2 = cmiRAVK2 + cmiEstRAVK2(coord_dists=coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
    
  }
  return(cmiRAVK2/N)
}

KNN.CMI_estimates_Parall = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiFP = 0
  cmiRAVK2 = 0
  for(i in 1:N){
    coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    x_coords = 1:length(xind)
    y_coords = (length(xind)+1):length(c(xind, yind))
    z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
    
    cmiFP = cmiFP + cmiEstFP(coord_dists=coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
    cmiRAVK2 = cmiRAVK2 + cmiEstRAVK2(coord_dists=coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
  }
  return(c(cmiFP/N,cmiRAVK2/N))
}

KNN.MI_estimates_MS_INF = function(data, xind, yind, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind)
  distArray = getPairwiseDistArray_II(data=as.matrix(data), inds=inds, isCat=isCat)
  miMS = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  for(i in 1:N){
    MI_coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    miMS = miMS + miEstMS(coord_dists=MI_coord_dists, k=k,  N=N, x_coords=x_coords, y_coords=y_coords)
  }
  return(miMS/N)
} 

KNN.CMI_estimates_MS = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiMS = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
  for(i in 1:N){
    CMI_coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    cmiMS = cmiMS + cmiEstMS(coord_dists=CMI_coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
  }
  return(cmiMS/N)
}

KNN.CMI_estimates_MS_INF = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_II(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiMS = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
  for(i in 1:N){
    CMI_coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    cmiMS = cmiMS + cmiEstMS(coord_dists=CMI_coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords)
  }
  return(cmiMS/N)
}

KNN.CMI_estimates_MIN = function(data, xind, yind, zinds, k= 7, isCat=c(), logE=T){
  N = dim(data)[1]
  inds = c(xind, yind, zinds)
  distArray = getPairwiseDistArray_IO(data=as.matrix(data), inds=inds, isCat=isCat)
  cmiMIN = 0
  x_coords = 1:length(xind)
  y_coords = (length(xind)+1):length(c(xind, yind))
  z_coords = (length(c(xind,yind))+1):length(c(xind, yind, zinds))
  if(length(isCat) != 0){
    for(m in 1:length(isCat)){
      if(isCat[m] %in% xind){
        isCat[m] = which(xind == isCat[m])
      }else if(isCat[m] %in% yind){
        isCat[m] = which(yind == isCat[m]) + length(xind)
      }else{
        isCat[m] = which(zinds == isCat[m]) + length(c(xind,yind))
      }
    }
  }
  for(i in 1:N){
    CMI_coord_dists = getPointCoordDists(distArray=distArray, point_i=i, inds=inds)
    cmiMIN = cmiMIN + cmiEstMin(coord_dists=CMI_coord_dists, k=k,  x_coords=x_coords, y_coords=y_coords, z_coords=z_coords, isCat=isCat)
  }
  return(cmiMIN/N)
}

