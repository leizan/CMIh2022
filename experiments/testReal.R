rm(list = ls())

source("source.R")

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


getCat = function(name){
  flag = strsplit(name,"_")[[1]][1]
  if( flag == 'conDependence'){
    return(c(1,3))
  }else{
    return(c(3))
  }
}

getZinds = function(name){
  flag = strsplit(name,"_")[[1]][1]
  if( flag == 'conDependence'){
    return(c(3))
  }else{
    return(c(3,4,5))
  }
}

nameExper = c("chain_2021-11-15_3000.csv", "chain_2021-11-15_8000.csv", "chain_2021-11-15_10000.csv", "chain_2021-12-13_4000.csv", "chain_2021-12-13_15000.csv",
              "chain_2021-12-13_23000.csv", "chain_2022-01-18_1000.csv", "chain_2022-01-18_4000.csv", "chain_2022-01-18_6000.csv", "chain_2022-02-23_20000.csv",
              "chain_2022-02-23_30000.csv", "chain_2022-02-23_40000.csv", 
              
              "fork_2021-11-15_3000.csv", "fork_2021-11-15_8000.csv", "fork_2021-11-15_10000.csv", "fork_2021-12-13_4000.csv", "fork_2021-12-13_15000.csv",
              "fork_2021-12-13_23000.csv", "fork_2022-01-18_1000.csv", "fork_2022-01-18_4000.csv", "fork_2022-01-18_6000.csv", "fork_2022-02-23_20000.csv",
              "fork_2022-02-23_30000.csv", "fork_2022-02-23_40000.csv",
              
              "conDependence_2021-11-15_3000.csv", "conDependence_2021-11-15_8000.csv", "conDependence_2021-11-15_10000.csv", "conDependence_2021-12-13_4000.csv",
              "conDependence_2021-12-13_15000.csv", "conDependence_2021-12-13_23000.csv", "conDependence_2022-01-18_1000.csv", "conDependence_2022-01-18_4000.csv",
              "conDependence_2022-01-18_6000.csv", "conDependence_2022-02-23_20000.csv", "conDependence_2022-02-23_30000.csv", "conDependence_2022-02-23_40000.csv")

indice  = 1
xind = c(1)
yind = c(2)

pValueMS_local = list()
for(name in nameExper){
  pValueMS_local[[name]] = c(0,0)
}

pValueMS_local = data.frame(pValueMS_local)
pValueCMIh_local = data.frame(pValueMS_local)

pValueMS_local_const = data.frame(pValueMS_local)
pValueCMIh_local_const = data.frame(pValueMS_local)

pValueGlobal_MS = data.frame(pValueMS_local)
pValueGlobal_CMIh = data.frame(pValueMS_local)


total = length(nameExper)
pb = txtProgressBar(min = 1, max = total, style = 3)
repeatTime = 1
k = 100 # 100
B = 1000 # 1000

for(name in nameExper){
  print(name)
  isCat = getCat(name)
  zinds = getZinds(name)
  
  resMS_local = rep(0,repeatTime)
  resCMIh_local = rep(0,repeatTime)

  resMS_local_const = rep(0,repeatTime)
  resCMIh_local_const = rep(0,repeatTime)

  resGlobal_MS = rep(0,repeatTime)
  resGlobal_CMIh = rep(0,repeatTime)

  for(m in 1:repeatTime){
    dd = read.csv(paste('ordinal_data/', name, sep=""))
    resMS_local[m] = localPermuTestMS_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local[m] = localPermuTestCMIh_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)
    
    resMS_local_const[m] = localPermuTestMS_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local_const[m] = localPermuTestCMIh_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)
  
    resGlobal_MS[m] = globalPermuTest_MS(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resGlobal_CMIh[m] = globalPermuTest_CMIh(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)
  }

  name = gsub('-', '.', name)
  pValueMS_local[name][1,]  = mean(resMS_local)
  pValueCMIh_local[name][1,]  = mean(resCMIh_local)
  
  pValueMS_local_const[name][1,]  = mean(resMS_local_const)
  pValueCMIh_local_const[name][1,]  = mean(resCMIh_local_const)

  pValueGlobal_MS[name][1,]  = mean(resGlobal_MS)
  pValueGlobal_CMIh[name][1,]  = mean(resGlobal_CMIh)

  indice = indice + 1

  setTxtProgressBar(pb, indice)

  res = data.frame(CMIh_local_mean = as.numeric(pValueCMIh_local[1,]),
                   CMIh_local_const_mean = as.numeric(pValueCMIh_local_const[1,]),
                   Global_CMIh_mean = as.numeric(pValueGlobal_CMIh[1,]),

                   ##################################################
                   ##################################################
                   
                   MS_local_mean = as.numeric(pValueMS_local[1,]),
                   MS_local_const_mean = as.numeric(pValueMS_local_const[1,]),
                   Global_MS_mean = as.numeric(pValueGlobal_MS[1,]),
                    
                   row.names = nameExper)
  
  write.table(res, file="results/res_real.tab", sep="\t", row.names=T, col.names=T, quote=F)
}

