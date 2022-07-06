rm(list = ls())

source("source.R")

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


getCat = function(name){
  flag = strsplit(name,"_")[[1]][2]
  if(flag== 'c1'){
    return(c(2,3))
  }else if(flag == 'c2'){
    return(c(3))
  }else if(flag== 'a1'){
    return(c(2))
  }else if(flag== 'b'){
    return(c(1,2))
  }else if(flag== 'd'){
    return(c(1,2,3))
  }else{
    return(c())
  }
}


functionList = c(Chain_a1, Chain_a2, Chain_b, Chain_c1, Chain_c2, Chain_d,

                Fork_a1, Fork_a2, Fork_b, Fork_c1, Fork_c2, Fork_d,

                 Collider_a1, Collider_a2, Collider_b, Collider_c1, Collider_c2, Collider_d)

nameExper = c("Chain_a1", "Chain_a2", "Chain_b", "Chain_c1", "Chain_c2", "Chain_d",

              "Fork_a1", "Fork_a2", "Fork_b", "Fork_c1", "Fork_c2", "Fork_d",

              "Collider_a1", "Collider_a2", "Collider_b", "Collider_c1", "Collider_c2", "Collider_d")


indice  = 1
xind = c(1)
yind = c(2)
zinds = c(3)
inds = c(xind, yind, zinds)

module_accept = list()
for(name in nameExper){
  module_accept[[name]] = c(0,0)
}

acceptMS_local = data.frame(module_accept)
acceptCMIh_local = data.frame(module_accept)

acceptMS_local_const = data.frame(module_accept)
acceptCMIh_local_const = data.frame(module_accept)

acceptGlobal_MS = data.frame(module_accept)
acceptGlobal_CMIh = data.frame(module_accept)

total = length(functionList)
pb = txtProgressBar(min = 0, max = total, style = 3)
set.seed(3344)
repeatTime = 10 # 10
n = 500 # 500
k = 50 # 50
B = 1000 # 1000

module_list = list()
for(name in nameExper){
  module_list[[name]] = rep(0,repeatTime)
}

pValueMS_local = data.frame(module_list)
pValueCMIh_local = data.frame(module_list)

pValueMS_local_const = data.frame(module_list)
pValueCMIh_local_const = data.frame(module_list)

pValueGlobal_MS = data.frame(module_list)
pValueGlobal_CMIh = data.frame(module_list)


for(exper in functionList){
  name = nameExper[indice]
  print(name)
  resMS_local = rep(0,repeatTime)
  resCMIh_local = rep(0,repeatTime)
  
  resMS_local_const = rep(0,repeatTime)
  resCMIh_local_const = rep(0,repeatTime)
  
  resGlobal_MS = rep(0,repeatTime)
  resGlobal_CMIh = rep(0,repeatTime)
  for(m in 1:repeatTime){
    dd = exper(n=n)
    isCat = getCat(name)
    
    resMS_local[m] = localPermuTestMS_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local[m] = localPermuTestCMIh_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    
    resMS_local_const[m] = localPermuTestMS_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local_const[m] = localPermuTestCMIh_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    
    resGlobal_MS[m] = globalPermuTest_MS(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resGlobal_CMIh[m] = globalPermuTest_CMIh(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
  }
  
  pValueMS_local[name]  = resMS_local
  acceptMS_local[name][1,] = length(which(resMS_local>=0.01))/repeatTime 
  acceptMS_local[name][2,] = length(which(resMS_local>=0.05))/repeatTime
  
  pValueCMIh_local[name]  = resCMIh_local
  acceptCMIh_local[name][1,] = length(which(resCMIh_local>=0.01))/repeatTime 
  acceptCMIh_local[name][2,] = length(which(resCMIh_local>=0.05))/repeatTime
  
  pValueMS_local_const[name]  = resMS_local_const
  acceptMS_local_const[name][1,] = length(which(resMS_local_const>=0.01))/repeatTime 
  acceptMS_local_const[name][2,] = length(which(resMS_local_const>=0.05))/repeatTime
  
  pValueCMIh_local_const[name] = resCMIh_local_const
  acceptCMIh_local_const[name][1,] = length(which(resCMIh_local_const>=0.01))/repeatTime 
  acceptCMIh_local_const[name][2,] = length(which(resCMIh_local_const>=0.05))/repeatTime
  
  pValueGlobal_MS[name]  = resGlobal_MS
  acceptGlobal_MS[name][1,] = length(which(resGlobal_MS>=0.01))/repeatTime 
  acceptGlobal_MS[name][2,] = length(which(resGlobal_MS>=0.05))/repeatTime
  
  pValueGlobal_CMIh[name]  = resGlobal_CMIh
  acceptGlobal_CMIh[name][1,] = length(which(resGlobal_CMIh>=0.01))/repeatTime 
  acceptGlobal_CMIh[name][2,] = length(which(resGlobal_CMIh>=0.05))/repeatTime
  
  if(strsplit(name,"_")[[1]][1] == "Collider"){

    acceptMS_local[name][1,] = 1 - length(which(resMS_local>=0.01))/repeatTime 
    acceptMS_local[name][2,] = 1 - length(which(resMS_local>=0.05))/repeatTime
    
    acceptCMIh_local[name][1,] = 1 - length(which(resCMIh_local>=0.01))/repeatTime 
    acceptCMIh_local[name][2,] = 1 - length(which(resCMIh_local>=0.05))/repeatTime
    
    acceptMS_local_const[name][1,] = 1 - length(which(resMS_local_const>=0.01))/repeatTime 
    acceptMS_local_const[name][2,] = 1 - length(which(resMS_local_const>=0.05))/repeatTime
    
    acceptCMIh_local_const[name][1,] = 1 - length(which(resCMIh_local_const>=0.01))/repeatTime 
    acceptCMIh_local_const[name][2,] = 1 - length(which(resCMIh_local_const>=0.05))/repeatTime
    
    acceptGlobal_MS[name][1,] = 1 - length(which(resGlobal_MS>=0.01))/repeatTime 
    acceptGlobal_MS[name][2,] = 1 - length(which(resGlobal_MS>=0.05))/repeatTime
    
    acceptGlobal_CMIh[name][1,] = 1 - length(which(resGlobal_CMIh>=0.01))/repeatTime 
    acceptGlobal_CMIh[name][2,] = 1 - length(which(resGlobal_CMIh>=0.05))/repeatTime
  }
  
  indice = indice + 1
  setTxtProgressBar(pb, indice)
  
  

  write.table(pValueMS_local, file="pValue/MS_local.tab", sep="\t", row.names=T, col.names=T, quote=F)
  write.table(pValueCMIh_local, file="pValue/CMIh_local.tab", sep="\t", row.names=T, col.names=T, quote=F)
  
  write.table(pValueMS_local_const, file="pValue/MS_local_const.tab", sep="\t", row.names=T, col.names=T, quote=F)
  write.table(pValueCMIh_local_const, file="pValue/CMIh_local_const.tab", sep="\t", row.names=T, col.names=T, quote=F)
  
  write.table(pValueGlobal_MS, file="pValue/Global_MS.tab", sep="\t", row.names=T, col.names=T, quote=F)
  write.table(pValueGlobal_CMIh, file="pValue/Global_CMIh.tab", sep="\t", row.names=T, col.names=T, quote=F)
  
  resAccept = data.frame(CMIh_local_001 = as.numeric(acceptCMIh_local[1,]),
                         CMIh_local_005 = as.numeric(acceptCMIh_local[2,]),
                         
                         CMIh_local_const_001 = as.numeric(acceptCMIh_local_const[1,]),
                         CMIh_local_const_005 = as.numeric(acceptCMIh_local_const[2,]),
                         
                         CMIh_global_001 = as.numeric(acceptGlobal_CMIh[1,]),
                         CMIh_global_005 = as.numeric(acceptGlobal_CMIh[2,]),
                         
                         MS_local_001 = as.numeric(acceptMS_local[1,]),
                         MS_local_005 = as.numeric(acceptMS_local[2,]),
                         
                         MS_local_const_001 = as.numeric(acceptMS_local_const[1,]),
                         MS_local_const_005 = as.numeric(acceptMS_local_const[2,]),
                         
                         MS_global_001 = as.numeric(acceptGlobal_MS[1,]),
                         MS_global_005 = as.numeric(acceptGlobal_MS[2,]),
                         
                         row.names = nameExper)
  write.table(resAccept, file="results/resAccept.tab", sep="\t", row.names=T, col.names=T, quote=F)
}
