rm(list = ls())

source("source.R")

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)



nameExper = c("Temparture")

indice  = 1

pValueMS_local = list()
for(name in nameExper){
  pValueMS_local[[name]] = c(0,0)
}

pValueMS_local = data.frame(pValueMS_local)
pValueCMIh_local = data.frame(pValueMS_local)

pValueMS_local_const = data.frame(pValueMS_local)
pValueCMIh_local_const = data.frame(pValueMS_local)

pValueMS_global = data.frame(pValueMS_local)
pValueCMIh_global = data.frame(pValueMS_local)


total = length(nameExper)
pb = txtProgressBar(min = 1, max = total, style = 3)
repeatTime = 1
k = 34 # 100
B = 1000 # 1000

xind = c(1)
yind = c(2)
zinds = c(3)
isCat = c(3)


for(name in nameExper){
  print(name)
  
  resMS_local = rep(0,repeatTime)
  resCMIh_local = rep(0,repeatTime)
  
  resMS_local_const = rep(0,repeatTime)
  resCMIh_local_const = rep(0,repeatTime)
  
  resMS_global = rep(0,repeatTime)
  resCMIh_global = rep(0,repeatTime)
  
  for(m in 1:repeatTime){
    raw_data_1 = read.table("ordinal_data/pair0003.txt",sep=" ",header=F)
    raw_data_2 = read.table("ordinal_data/pair0020.txt",sep=" ",header=F)
    conTem = raw_data_1$V2
    orderTem = conTem[order(conTem)]
    ordiTem = rep(0, length(conTem))
    for(i in 1:length(conTem)){
      if(conTem[i] <= orderTem[116]){
        ordiTem[i] = 0
      }else if(conTem[i] > orderTem[232]){
        ordiTem[i] = 2
      }else{
        ordiTem[i] = 1
      }
    }
    dd = data.frame(S1 = raw_data_1[,1],  
                    S2 = raw_data_2[,1],  
                    S5 = ordiTem)

    resMS_local[m] = localPermuTestMS_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local[m] = localPermuTestCMIh_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)
     
    resMS_local_const[m] = localPermuTestMS_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local_const[m] = localPermuTestCMIh_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)
    
    resMS_global[m] = globalPermuTest_MS_without_condition(dd=dd, xind=xind, yind=yind, zinds=c(), isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_global[m] = globalPermuTest_CMIh(dd=dd, xind=xind, yind=yind, zinds=c(), isCat=isCat, B=B, kPer=5)
    
  }
  
  pValueMS_local[name][1,]  = mean(resMS_local)
  pValueCMIh_local[name][1,]  = mean(resCMIh_local)

  pValueMS_local_const[name][1,]  = mean(resMS_local_const)
  pValueCMIh_local_const[name][1,]  = mean(resCMIh_local_const)

  pValueMS_global[name][1,]  = mean(resMS_global)
  pValueCMIh_global[name][1,]  = mean(resCMIh_global)


  indice = indice + 1

  setTxtProgressBar(pb, indice)

  res = data.frame(CMIh_local_mean = as.numeric(pValueCMIh_local[1,]),
                   CMIh_local_const_mean = as.numeric(pValueCMIh_local_const[1,]),
                   CMIh_global = as.numeric(pValueCMIh_global[1,]),

                   ##################################################
                   ##################################################

                   MS_local_mean = as.numeric(pValueMS_local[1,]),
                   MS_local_const_mean = as.numeric(pValueMS_local_const[1,]),
                   MS_global = as.numeric(pValueMS_global[1,]),

                   row.names = nameExper)

  write.table(res, file="results/res_Temperature.tab", sep="\t", row.names=T, col.names=T, quote=F)
}

