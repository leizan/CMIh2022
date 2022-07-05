#### Experiment X in the paper

rm(list = ls())
source("source.R")

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

loops = 1:100 #1:100
n = 2000 #2000
nz = 0:4 # 0:4
coeff = 0.1
k = round(coeff*n) 

## Test 1: uniform dependence
mean_vec = data.frame(NZ=nz, LH=rep(0, length(nz)), CMIh=rep(0, length(nz)), MS=rep(0, length(nz)), FP=rep(0, length(nz)),RAVK=rep(0, length(nz)))
sd_vec = data.frame(NZ=nz, LH=rep(0, length(nz)), CMIh=rep(0, length(nz)), MS=rep(0, length(nz)), FP=rep(0, length(nz)), RAVK=rep(0, length(nz)))
meanse_vec = data.frame(NZ=nz, LH=rep(0, length(nz)), CMIh=rep(0, length(nz)), MS=rep(0, length(nz)), FP=rep(0, length(nz)), RAVK=rep(0, length(nz)))
mean_time = data.frame(NZ=nz, LH=rep(0, length(nz)), CMIh=rep(0, length(nz)), MS=rep(0, length(nz)), FP=rep(0, length(nz)), RAVK=rep(0, length(nz)))
sd_time = data.frame(NZ=nz, LH=rep(0, length(nz)), CMIh=rep(0, length(nz)), MS=rep(0, length(nz)), FP=rep(0, length(nz)), RAVK=rep(0, length(nz)))

m = 5
trueI = log(m) - (m-1)*log(2)/m
index = 1
set.seed(1)
for(s in nz){
    print(s)
    estimates = data.frame(LH=rep(0, length(loops)), CMIh=rep(0, length(loops)), MS=rep(0, length(loops)), FP=rep(0, length(loops)), RAVK=rep(0, length(loops)))
    times = data.frame(LH=rep(0, length(loops)), CMIh=rep(0, length(loops)), MS=rep(0, length(loops)), FP=rep(0, length(loops)), RAVK=rep(0, length(loops)))
    for(l in loops){
        dd = test2(n,m=m)
        if(s > 0){
          for(gz in 1:s){
            dd = data.frame(dd, z=rbinom(n, p=0.5, size=3))
          }
        }
        estimate = NULL
        if(s == 0){
          isCat=c(1)
          start1 = Sys.time()
          estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(1),logE=T,number_cores_use = cores[1]-1)
          times[l,1] = difftime(Sys.time(), start1, units = "secs")
          
          start2 = Sys.time()
          estimates[l,2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(),isCat=c(1))
          times[l,2] = difftime(Sys.time(), start2, units = "secs")
          
          start3 = Sys.time()
          estimates[l,3] = KNN.MI_estimates_MS(data=dd, xind=1,yind=2,k=k,isCat=c(1))
          times[l,3] = difftime(Sys.time(), start3, units = "secs")
          
          start4 = Sys.time()
          estimates[l,4] = KNN.MI_estimates_FP(data=dd, xind=1,yind=2,k=k,isCat=c(1))
          times[l,4] = difftime(Sys.time(), start4, units = "secs")
          
          start5 = Sys.time()
          estimates[l,5] = KNN.MI_estimates_RAVK2(data=dd, xind=1,yind=2,k=k,isCat=c(1))
          times[l,5] = difftime(Sys.time(), start5, units = "secs")
        }else{
          isCat=c(1,1:s + 2)
          start1 = Sys.time()
          estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(1:s + 2),isCat=c(1,1:s + 2),logE=T,number_cores_use = cores[1]-1)
          times[l,1] = difftime(Sys.time(), start1, units = "secs")
          
          start2 = Sys.time()
          estimates[l,2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(1:s + 2),isCat=c(1,1:s + 2))
          times[l,2] = difftime(Sys.time(), start2, units = "secs")
          
          start3 = Sys.time()
          estimates[l,3] = KNN.CMI_estimates_MS(data=dd, xind=1,yind=2,zinds=c(1:s + 2),k=k,isCat=c(1,1:s + 2))
          times[l,3] = difftime(Sys.time(), start3, units = "secs")
          
          start4 = Sys.time()
          estimates[l,4] = KNN.CMI_estimates_FP(data=dd, xind=1,yind=2,zinds=c(1:s + 2),k=k,isCat=c(1,1:s + 2))
          times[l,4] = difftime(Sys.time(), start4, units = "secs")
          
          start5 = Sys.time()
          estimates[l,5] = KNN.CMI_estimates_RAVK2(data=dd, xind=1,yind=2,zinds=c(1:s + 2),k=k,isCat=c(1,1:s + 2))
          times[l,5] = difftime(Sys.time(), start5, units = "secs")        
        }
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
       
    }
    for(i in 1:5){
      mean_vec[index,i+1] = mean(estimates[,i])
      sd_vec[index,i+1] = sd(estimates[,i])
      meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / length(loops)
      mean_time[index,i+1] = mean(times[,i])
      sd_time[index,i+1] = sd(times[,i])
    }
    print(mean_vec[index,])
    print(sd_vec[index,])
    print(meanse_vec[index,])
    print(mean_time[index,])
    print(sd_time[index,])
    index = index + 1
}
write.table(mean_vec, file="results/testmd_mean_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/testmd_sd_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/testmd_meanse_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(mean_time, file="results/testmd_mean_time_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_time, file="results/testmd_sd_time_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)