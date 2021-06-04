#Calculation of Respiration Rate and CPI
#setwd("Set Working directory to folder containing data")

mod.out.fun = function(max.time,min.time,dat,unique.id.use) {
  
  mod = dat$DO_ppm[which(dat$time_min <= max.time & dat$time_min >= min.time  & dat$unique.id == unique.id.use)] ~ dat$time_min[which(dat$time_min <= max.time & dat$time_min >= min.time & dat$unique.id == unique.id.use)]
  mod.sum = summary(lm(mod))
  mod.rate = -mod.sum$coefficients[2,1]
  mod.rsq = mod.sum$r.squared
  return(c(mod.rate,mod.rsq))
  
}


fac.to.char.fun = function(matrix.in) {
  
  i = sapply(matrix.in,is.factor)
  matrix.in[i] = lapply(matrix.in[i],as.character)
  
  return(matrix.in)
  
}
path.in = "//Data Package/01_Rates/"
path.out = "//Data Package/01_Rates/"

dat = fac.to.char.fun(read.csv(file = paste(path.in,"Stegen_EC_Raw_DO.csv",sep="")))
head(dat)

dat$cycles = dat$Sample
dat$cycles = gsub(pattern = " cycles",replacement = "",x = dat$cycles)
dat$cycles = as.numeric(gsub(pattern = c(" cycle"),replacement = "",x = dat$cycles))

dat$unique.id = paste(dat$ID,dat$Sample,sep="_")

head(dat)
str(dat)

id.and.sample = unique(dat$unique.id)

dat.rates = numeric()

for (i in id.and.sample) {
  
  four.pt.mod = mod.out.fun(max.time = 90,min.time = 0,dat = dat,unique.id.use = i)
  dat.rates = rbind(dat.rates,c(i,unique(dat$cycles[dat$unique.id == i]),four.pt.mod))
  
}

colnames(dat.rates) = c("unique.id","cycles","four.pt.rate","four.pt.rsq")
dat.rates = as.data.frame(dat.rates)
dat.rates = fac.to.char.fun(dat.rates)
dat.rates[2:ncol(dat.rates)] = lapply(dat.rates[2:ncol(dat.rates)],as.numeric)
str(dat.rates)
head(dat.rates)

## save dat.rates
write.csv(dat.rates, file = paste(path.out,"Stegen_EC_DO_Rates.csv",sep=""),row.names = F,quote = F)
##

var.levels = 0:5 # number of cycles of experimental wet/dry

cont.point.comp = numeric()

for (i in var.levels) {
  
  temp.rates = dat.rates$four.pt.rate[which(dat.rates$cycles == i)]
  cont.point.comp = rbind(cont.point.comp,c(i,sum(temp.rates[which(temp.rates > median(temp.rates))]) / sum(temp.rates)))
  
}

colnames(cont.point.comp) = c("cycles","bgc.above.med")
cont.point.comp = as.data.frame(cont.point.comp)
cont.point.comp

write.csv(cont.point.comp,paste(path.out,"Stegen_EC_CPIs.csv",sep=""),row.names = F,quote = F)

