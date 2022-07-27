#### Experiments I-XI in the paper

rm(list = ls())
source("source.R")


cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

steps = 1:20 # 1:20
loops = 100 # 1:100
fact = 100 # 100
coeff = 0.1 

# test 1: two dim Gaussian with covariance tCov

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

tCov = 0.6
trueI = -0.5 * log(1-tCov^2)
covM = matrix(c(1,tCov,tCov,1), nrow=2, ncol=2)
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = mvrnorm(n=n, Sigma=covM, mu=c(0,0))
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(),isCat=c())
    estimation[3] = KNN.MI_estimates_MS(data=dd, xind=1,yind=2,k=k,isCat=c())
    res = KNN.MI_estimates_Parall(data=dd, xind=1,yind=2,k=k,isCat=c())
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test1_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test1_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test1_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 2: uniform dependence

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

m = 5
trueI = log(m) - (m-1)*log(2)/m
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test2(n,m=m)
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(1),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(),isCat=c(1))
    estimation[3] = KNN.MI_estimates_MS(data=dd, xind=1,yind=2,k=k,isCat=c(1))
    res = KNN.MI_estimates_Parall(data=dd, xind=1,yind=2,k=k,isCat=c(1))
    estimation[4] = res[1]
    estimation[5] = res[2]

    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test2_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test2_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test2_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 3: zero-inflated poisson dependence

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

p = 0.15
trueI = (1-p) * 0.3012
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test3(n,p=p)
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(2),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(),isCat=c(2))
    estimation[3] = KNN.MI_estimates_MS(data=dd, xind=1,yind=2,k=k,isCat=c(2))
    res = KNN.MI_estimates_Parall(data=dd, xind=1,yind=2,k=k,isCat=c(2))
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test3_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test3_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test3_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 4: extension of test 1, but not affects the true value

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

tCov = 0.6
trueI = -0.5 * log(1-tCov^2)
covM = matrix(c(1,tCov,tCov,1), nrow=2, ncol=2)
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = mvrnorm(n=n, Sigma=covM, mu=c(0,0))
    dd = data.frame(dd, z=rbinom(n, p=0.5, size=3))
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(3),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(3))
    estimation[3] = KNN.CMI_estimates_MS(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(3))
    res = KNN.CMI_estimates_Parall(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(3))
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test4_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test4_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test4_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 5: extension of test 2, but not affects the true value

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

m = 5
trueI = log(m) - (m-1)*log(2)/m
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test5(n,m=m)
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(1,3),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(1,3))
    estimation[3] = KNN.CMI_estimates_MS(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(1,3))
    res = KNN.CMI_estimates_Parall(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(1,3))
    estimation[4] = res[1]
    estimation[5] = res[2]

    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test5_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test5_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test5_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 6: extension of test 3, but not affects the true value.

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

p = 0.15
trueI = (1-p) * 0.3012
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test6(n,p=p)
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(2,3),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(2,3))
    estimation[3] = KNN.CMI_estimates_MS(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(2,3))
    res = KNN.CMI_estimates_Parall(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(2,3))
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test6_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test6_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test6_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 7: conditional independence

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

trueI = 0
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test7(n)
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(3),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(3))
    estimation[3] = KNN.CMI_estimates_MS(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(3))
    res = KNN.CMI_estimates_Parall(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(3))
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test7_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test7_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test7_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 8: conditional independence

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

trueI = 0
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test8(n)
    estimate = CMI.estimates(data=dd, xind=2,yind=3,zinds=c(1),isCat=c(1,3),logE=T,number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=2,yind=3,zinds=c(1),isCat=c(1,3))
    estimation[3] = KNN.CMI_estimates_MS(data=dd, xind=2,yind=3,zinds=c(1),k=k,isCat=c(1,3))
    res = KNN.CMI_estimates_Parall(data=dd, xind=2,yind=3,zinds=c(1),k=k,isCat=c(1,3))
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test8_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test8_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test8_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 9: conditional independence

mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

trueI = 0
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,5)
    dd = test9(n)
    estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(2,3), logE=T, number_cores_use=cores[1]-1)
    estimation[1] = estimate[[1]] / n ## pure likelihood
    estimation[2] = mixedEstimator(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(2,3))
    estimation[3] = KNN.CMI_estimates_MS(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(2,3))
    res = KNN.CMI_estimates_Parall(data=dd, xind=1,yind=2,zinds=c(3),k=k,isCat=c(2,3))
    estimation[4] = res[1]
    estimation[5] = res[2]
    estimation
  }
  for(i in 1:5){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test9_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test9_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test9_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 10: multidimensional conditional independent 

mean_vec = data.frame(Samples=steps*fact, CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

m = 5
trueI = 0
index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,4)
    dd = test10(n,m)
    estimation[1] = mixedEstimator(data=dd, xind=c(1,2,3),yind=c(4),zinds=c(5,6,7,8),isCat=c(3,4,5,6))
    estimation[2] = KNN.CMI_estimates_MS(data=dd, xind=c(1,2,3),yind=c(4),zinds=c(5,6,7,8),k=k,isCat=c(3,4,5,6))
    res = KNN.CMI_estimates_Parall(data=dd, xind=c(1,2,3),yind=c(4),zinds=c(5,6,7,8),k=k,isCat=c(3,4,5,6))
    estimation[3] = res[1]
    estimation[4] = res[2]
    estimation
  }
  for(i in 1:4){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test10_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test10_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test10_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)



# test 11: multidimensional conditional dependent 

mean_vec = data.frame(Samples=steps*fact, CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, CMIh=rep(0, length(steps)), MS=rep(0, length(steps)), FP=rep(0, length(steps)), RAVK=rep(0, length(steps)))

tCov = 0.6
m = 5
p = 0.15
trueI = -0.5 * log(1-tCov^2) + log(m) - (m-1)*log(2)/m + (1-p) * 0.3012

index = 1
set.seed(1)
for(s in steps){
  print(s)
  n = s * fact
  k = round(coeff*n)
  estimates = foreach(l = 1:loops, .combine=rbind) %dopar% {
    source("source.R")
    estimation = rep(0,4)
    dd = test11(n=n,tCov=tCov,m=m,p=p)
    estimation[1] = mixedEstimator(data=dd, xind=c(1,2,3),yind=c(4,5,6),zinds=c(),isCat=c(2,6))
    estimation[2] = KNN.MI_estimates_MS(data=dd, xind=c(1,2,3),yind=c(4,5,6),k=k,isCat=c(2,6))
    res = KNN.MI_estimates_Parall(data=dd, xind=c(1,2,3),yind=c(4,5,6),k=k,isCat=c(2,6))
    estimation[3] = res[1]
    estimation[4] = res[2]
    estimation
  }
  for(i in 1:4){
    mean_vec[index,i+1] = mean(estimates[,i])
    sd_vec[index,i+1] = sd(estimates[,i])
    meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / loops
  }
  print(mean_vec[index,])
  print(sd_vec[index,])
  print(meanse_vec[index,])
  index = index + 1
}
write.table(mean_vec, file='results/test11_mean.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file='results/test11_sd.tab', sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file='results/test11_meanse.tab', sep="\t", row.names=F, col.names=T, quote=F)
