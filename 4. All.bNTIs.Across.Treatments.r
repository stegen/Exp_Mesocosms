# This is for making box plots of bNTI from amplicon data and from cDNA data, across treatments, and plotting with respect to rates and CPI

rm(list=ls())

############
# change these things

# define some types and names
bnti.type = 'cDNA'  # 'cDNA' or 'gDNA' 
bnti.name =  'cDNA_bNTI_weighted_rare_27227' #'cDNA_bNTI_weighted_rare_27227'#  or 'gDNA_bNTI_weighted_rare_15106'

#############

# path for output
path.out = "D:/Experimental Mesocosms/Data Package/02_bNTIcalculations"

###### load functions

fac.to.char.fun = function(matrix.in) {
  
  i = sapply(matrix.in,is.factor)
  matrix.in[i] = lapply(matrix.in[i],as.character)
  
  return(matrix.in)
  
}

### read in metadata 
#change the cDNA/gDNA sample ID to the information provided in column J of the metadata file.
meta = fac.to.char.fun(read.csv("D:/Experimental Mesocosms/Data Package/02_bNTIcalculations/Metadata.csv"))
str(meta)
head(meta)

### read in rate data
rate.in = read.csv("D:/Experimental Mesocosms/Data Package/01_Rates/Stegen_EC_DO_Rates.csv")
head(rate.in)

cpi.in = read.csv("D:/Experimental Mesocosms/Data Package/01_Rates/Stegen_EC_CPIs.csv")
cpi.in$cycles = paste(cpi.in$cycles,"cycles",sep=""); cpi.in$cycles[which(cpi.in$cycles == '1cycles')] = '1cycle' # changing names to match the other data
cpi.in$num.cycles = as.numeric(gsub("cycle.*","",cpi.in$cycles))
cpi.in

### read in the control point contribution data
cpc.in = read.csv("D:/Experimental Mesocosms/Data Package/01_Rates/Stegen_EC_DO_Rates.csv")
cpc.in$cycles = paste(cpc.in$cycles,"cycles",sep=""); cpc.in$cycles[which(cpc.in$cycles == '1cycles')] = '1cycle' # changing names to match the other data
cpc.in$sample.id = gsub("_.*","",cpc.in$unique.id)
length(unique(cpc.in$sample.id)) - nrow(cpc.in) # should be 0
cpc.in$cDNA = NA
for (i in 1:nrow(cpc.in)) {
  
  cpc.in$cDNA[i] = meta$cDNA[which(meta$Sample.ID == cpc.in$sample.id[i])]
  
}
head(cpc.in)

pdf(paste(path.out,"Rate.Box.CPI.Across.Treats_.pdf",sep=""),width = 8,height = 7)
  par(pty="s")
  boxplot(cpc.in$four.pt.rate ~ cpc.in$cycles,xlab="Experimental Treatment",ylab="Rate of O2 Consumption",cex.lab=1.5,cex.axis=1); abline(h=c(2,-2),col=2,lwd=2)
  par(new=T,pty="s")
  plot(cpi.in$bgc.above.med ~ cpi.in$num.cycles,ylim=c(0.5,max(cpi.in$bgc.above.med)),pch=19,typ="b",axes=F,ylab="",xlab="",col=4) # the order of cpi.in must have increasing cycles
  axis(side = 4,at = seq(0.5,0.9,by=0.1))
  mtext("Control Point Influence", side = 4, line = 3,cex = 1.5)
dev.off()

cpc.out = cpc.in

# paths to bNTI output

amplicon.bnti.path.in = "D://Experimental Mesocosms/Data Package/02_bNTIcalculations/"

# read in bNTI data

if (bnti.type %in% c('cDNA')) {
  
  bnti.in = read.csv(paste(amplicon.bnti.path.in,bnti.name,".csv",sep=""),row.names = 1)
  dim(bnti.in)
  bnti.in[1:5,1:5]
  
  for (i in 1:nrow(bnti.in)) {
    
    print(rownames(bnti.in)[i])
    rownames(bnti.in)[i] = meta$cDNA[which( meta[,bnti.type] == rownames(bnti.in)[i]  )]
    colnames(bnti.in)[i] = meta$cDNA[which( meta[,bnti.type] == colnames(bnti.in)[i]  )]
    print(rownames(bnti.in)[i])
    
  }
  
}
  
  # cpc analyses
  cpc.out$median.bnti.full = NA
  cpc.out$median.bnti.treat = NA
  cpc.out$pt.color = NA
  
  # put in codes for colors based on treatments, for plotting
  col.temp = 1
  for (i in unique(cpc.out$cycles)) {
    
    cpc.out$pt.color[grep(i,cpc.out$cycles)] = col.temp
    col.temp = col.temp + 1
  
  }
  head(cpc.out)
  
  bnti.for.cpc = bnti.in[which(rownames(bnti.in) %in% cpc.out$cDNA),which(rownames(bnti.in) %in% cpc.out$cDNA)]
  bnti.for.cpc[1:5,1:5]
  rownames(bnti.for.cpc)
  
  cpc.out = cpc.out[which(cpc.out$cDNA %in% rownames(bnti.for.cpc)),]
  dim(cpc.out)
  
  for (i in 1:nrow(cpc.out)) {
    
    cpc.out$median.bnti.full[i] = median(c(as.vector(bnti.for.cpc[,cpc.out$cDNA[i]]),as.vector(t(bnti.for.cpc[cpc.out$cDNA[i],]))),na.rm = T)
    cpc.out$median.bnti.treat[i] = median(c(as.vector(bnti.for.cpc[grep(cpc.out$cycles[i],rownames(bnti.for.cpc)),cpc.out$cDNA[i]]),as.vector(t(bnti.for.cpc[cpc.out$cDNA[i],grep(cpc.out$cycles[i],colnames(bnti.for.cpc))]))),na.rm = T)
    
  }
  head(cpc.out)
#save cpcout_cDNA
write.csv(cpc.out, file="cpcout_cdNA.csv", row.names=FALSE)  
  
plot(cpc.out$four.pt.rate ~ cpc.out$median.bnti.full,pch=19,col=cpc.out$pt.color)
  plot(cpc.out$four.pt.CPC.within.treat ~ cpc.out$median.bnti.treat,pch=19,col=cpc.out$pt.color)
  plot(cpc.out$four.pt.CPC.percent.dev.treat ~ cpc.out$median.bnti.treat,pch=19,col=cpc.out$pt.color)
  
  pdf(paste(path.out,"Rate.vs.bNTI.Treat_",bnti.name,".pdf",sep=""))
    mod.to.plot = log(cpc.out$four.pt.rate+ I(min(cpc.out$four.pt.CPC[cpc.out$four.pt.rate > 0])/2)) ~ (abs(cpc.out$median.bnti.treat))
    mod = summary(lm(mod.to.plot)); mod
    plot(mod.to.plot,ylab="Ln(Rate of O2 Consumption)",xlab=paste('Abs(',bnti.name,")",sep=""),cex.lab=1.5,pch=19,cex=1.3)
    abline(mod,lwd=2,col=2)
    mtext(text = paste("R2 = ",round(mod$r.squared,digits = 2)," ",sep=""),side=3,line = -1.2,las=1,adj=1)
    mtext(text = paste("p = ",round(mod$coefficients[2,4],digits = 6)," ",sep=""),side=3,line = -2,las=1,adj=1)
  dev.off()
  
  # box plots
  bnti.for.box = numeric()
  
  for (i in c('0cycles','1cycles','2cycles','3cycles','4cycles','5cycles')) {
    
    bnti.temp = bnti.in[grep(pattern = i,x = rownames(bnti.in)),grep(pattern = i,x = colnames(bnti.in))]
    print(c(i,identical(rownames(bnti.temp),colnames(bnti.temp))))
    print(dim(bnti.temp))
    print(rownames(bnti.temp))
    
    bnti.for.box = rbind(bnti.for.box,cbind(as.vector(as.dist(bnti.temp)),i))
    
    rm('bnti.temp')
    
    
  }
  
  bnti.for.box = fac.to.char.fun(as.data.frame(bnti.for.box))
  colnames(bnti.for.box) = c('bNTI','Treatment')
  bnti.for.box$bNTI = as.numeric(bnti.for.box$bNTI)
  str(bnti.for.box)
  head(bnti.for.box)

##save bnti.for.box to get gDNA bNTI and cDNA bNTI values
write.csv(bnti.for.box, file="bNTIbox_cdNA.csv")
  
  pdf(paste(path.out,"Box.Across.Treats_",bnti.name,".pdf",sep=""),width = 8,height = 7)
    par(pty="s")
    boxplot(bnti.for.box$bNTI ~ bnti.for.box$Treatment,xlab="Experimental Treatment",ylab=bnti.name,cex.lab=1.5,cex.axis=1); abline(h=c(2,-2),col=2,lwd=2)
    par(new=T,pty="s")
    plot(cpi.in$bgc.above.med ~ cpi.in$num.cycles,ylim=c(0.5,max(cpi.in$bgc.above.med)),pch=19,typ="b",axes=F,ylab="",xlab="",col=4) # the order of cpi.in must have increasing cycles
    axis(side = 4,at = seq(0.5,0.9,by=0.1))
    mtext("Control Point Influence", side = 4, line = 3,cex = 1.5)
  dev.off()
  
  mean.bnti = (tapply(bnti.for.box$bNTI,INDEX = bnti.for.box$Treatment,FUN = 'mean')); mean.bnti = as.data.frame(mean.bnti)
  median.bnti = (tapply(bnti.for.box$bNTI,INDEX = bnti.for.box$Treatment,FUN = 'median')); median.bnti = as.data.frame(median.bnti)
  var.bnti = (tapply(bnti.for.box$bNTI,INDEX = bnti.for.box$Treatment,FUN = 'var')); var.bnti = as.data.frame(var.bnti)
  
  treat.metrics = merge(mean.bnti,median.bnti,by=0)
  treat.metrics = merge(treat.metrics,var.bnti,by.x='Row.names',by.y=0)
  treat.metrics = merge(treat.metrics,cpi.in,by.x='Row.names',by.y='cycles')
  treat.metrics
  
  mod.to.plot = treat.metrics$bgc.above.med ~ treat.metrics$median.bnti
  xlab = paste("Median ",bnti.name,sep="")
  ylab = "Control Point Influence"
  mod = summary(lm(mod.to.plot)); mod
  
  pdf(paste(path.out,"CPI.v.bNTI_",bnti.name,".pdf",sep=""))
    par(pty="s")  
    plot(mod.to.plot,xlab=xlab,ylab=ylab,cex.lab=1.5,cex=1.5,pch=19)
    abline(mod,lwd=2,col=2)
    mtext(text = paste("R2 = ",round(mod$r.squared,digits = 2)," ",sep=""),side=3,line = -1.2,las=1,adj=1)
    mtext(text = paste("p = ",round(mod$coefficients[2,4],digits = 2)," ",sep=""),side=3,line = -2,las=1,adj=1)
  dev.off()

