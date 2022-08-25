rm(list = ls())

source("source.R")

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)



nameExper = c("ADHD_GAH", "ADHD_HAM")

indice  = 0

pValueMS_local = list()
for(name in nameExper){
  pValueMS_local[[name]] = c(0,0)
}

pValueMS_local = data.frame(pValueMS_local)
pValueCMIh_local = data.frame(pValueMS_local)

pValueMS_local_const = data.frame(pValueMS_local)
pValueCMIh_local_const = data.frame(pValueMS_local)

total = length(nameExper)
pb = txtProgressBar(min = 0, max = total, style = 3)
repeatTime = 1
k = 42 # 35
B = 1000 # 1000

xind = c(1)
yind = c(2)
zinds = c(3)
isCat = c(1)


listCSV = c('Brown_TestRelease_phenotypic.csv', 'KKI_phenotypic.csv', 'NYU_phenotypic.csv', 'OHSU_phenotypic.csv',
            'OHSU_TestRelease_phenotypic.csv', 'Peking_1_phenotypic.csv', 'Peking_1_TestRelease_phenotypic.csv', 'Pittsburgh_phenotypic.csv')

raw_data = data.frame()
for(file in listCSV){
  csvFile = read.csv(paste("ordinal_data/nicholsn-adhd-200/", file, sep=""))
  csvFile = dplyr::select(csvFile, Gender, Hyper.Impulsive, Med.Status, Inattentive)%>%
    filter(Med.Status != 'pending' & Med.Status > 0 & Inattentive > 0 & Inattentive != 'N/A')
  raw_data = rbind(raw_data, csvFile)
}




for(name in nameExper){
  print(name)
  if(name == "ADHD_GAH"){
    dd = data.frame(G=raw_data[,1],
                    HI=raw_data[,2],
                    IN=raw_data[,4])
  }else if(name == "ADHD_HAM"){
    dd = data.frame(MED=raw_data[,3],
                    HI=raw_data[,2],
                    IN=raw_data[,4])
  }

  resMS_local = rep(0,repeatTime)
  resCMIh_local = rep(0,repeatTime)

  resMS_local_const = rep(0,repeatTime)
  resCMIh_local_const = rep(0,repeatTime)

  for(m in 1:repeatTime){

    resMS_local[m] = localPermuTestMS_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local[m] = localPermuTestCMIh_local(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)

    resMS_local_const[m] = localPermuTestMS_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5, k=k)
    resCMIh_local_const[m] = localPermuTestCMIh_local_const(dd=dd, xind=xind, yind=yind, zinds=zinds, isCat=isCat, B=B, kPer=5)

  }


  pValueMS_local[name][1,]  = mean(resMS_local)
  pValueCMIh_local[name][1,]  = mean(resCMIh_local)

  pValueMS_local_const[name][1,]  = mean(resMS_local_const)
  pValueCMIh_local_const[name][1,]  = mean(resCMIh_local_const)
  
  indice = indice + 1
  
  setTxtProgressBar(pb, indice)


  res = data.frame(CMIh_local_mean = as.numeric(pValueCMIh_local[1,]),
                   CMIh_local_const_mean = as.numeric(pValueCMIh_local_const[1,]),

                   ##################################################
                   ##################################################

                   MS_local_mean = as.numeric(pValueMS_local[1,]),
                   MS_local_const_mean = as.numeric(pValueMS_local_const[1,]),

                   row.names = nameExper)

  write.table(res, file="results/res_adhd.tab", sep="\t", row.names=T, col.names=T, quote=F)
}

