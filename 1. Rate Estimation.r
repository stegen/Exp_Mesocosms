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
path.in = "//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/01_Rates/"
path.out = "//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/01_Rates/"

dat = fac.to.char.fun(read.csv(file = paste(path.in,"Stegen_EC_Raw_DO.csv",sep="")))
dat$unique.id = paste(dat$Sample.ID,dat$Cycles,sep="_")
head(dat)

id.and.sample = unique(dat$unique.id)

dat.rates = numeric()

for (i in id.and.sample) {
  
  four.pt.mod = mod.out.fun(max.time = 90,min.time = 0,dat = dat,unique.id.use = i)
  dat.rates = rbind(dat.rates,c(i,four.pt.mod))
  
}

colnames(dat.rates) = c("unique.id","four.pt.rate","four.pt.rsq")
dat.rates = as.data.frame(dat.rates)
dat.rates = fac.to.char.fun(dat.rates)
dat.rates[2:ncol(dat.rates)] = lapply(dat.rates[2:ncol(dat.rates)],as.numeric)
str(dat.rates)
head(dat.rates)

## save dat.rates
write.csv(dat.rates, file = paste(path.out,"Stegen_EC_DO_Rates.csv",sep=""),row.names = F,quote = F)
##

var.levels = c('0 cycles','1 cycle','2 cycles','3 cycles','4 cycles','5 cycles') # number of cycles of experimental wet/dry

cont.point.comp = numeric()

for (i in var.levels) {
  
  temp.rates = dat.rates$four.pt.rate[grep(pattern = i,x = dat.rates$unique.id)]
  cont.point.comp = rbind(cont.point.comp,c(i,sum(temp.rates[which(temp.rates > median(temp.rates))]) / sum(temp.rates)))
  
}

colnames(cont.point.comp) = c("Cycles","CPI")
cont.point.comp = as.data.frame(cont.point.comp)
cont.point.comp

write.csv(cont.point.comp,paste(path.out,"Stegen_EC_CPIs.csv",sep=""),row.names = F,quote = F)

