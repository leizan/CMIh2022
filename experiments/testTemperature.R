rm(list = ls())

source("source.R")

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)



nameExper = c("Temparture")


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
    raw_data = read.table("ordinal_data/pair0003_pair0020.txt",sep=" ",header=F)
    conTem = raw_data$V3
    orderTem = conTem[order(conTem)]
    ordiTem = rep(0, length(conTem))
    for(i in 1:length(conTem)){
      if(conTem[i] <= orderTem[as.integer(length(conTem)/3)]){
        ordiTem[i] = 0
      }else if(conTem[i] > orderTem[as.integer(length(conTem)/3*2)]){
        ordiTem[i] = 2
      }else{
        ordiTem[i] = 1
      }
    }
    dd = data.frame(longitude = raw_data[,1],  
                    latitude = raw_data[,2],  
                    temperature = ordiTem)

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


  res = data.frame(CMIh_local_mean = as.numeric(pValueCMIh_local[1,]),
                   CMIh_local_const_mean = as.numeric(pValueCMIh_local_const[1,]),
                   CMIh_global = as.numeric(pValueCMIh_global[1,]),

                   ##################################################
                   ##################################################

                   MS_local_mean = as.numeric(pValueMS_local[1,]),
                   MS_local_const_mean = as.numeric(pValueMS_local_const[1,]),
                   MS_global = as.numeric(pValueMS_global[1,]),

                   row.names = nameExper)

  write.table(res, file="results/res_temperature.tab", sep="\t", row.names=T, col.names=T, quote=F)
}

