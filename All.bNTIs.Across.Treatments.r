# run the script once for gDNA and once again for cDNA
rm(list=ls())

############
# change these things

# define some types and names
bnti.type = 'Sample_ID_OTU_Table_gDNA' # Sample_ID_OTU_Table_cDNA or Sample_ID_OTU_Table_gDNA
if (bnti.type == 'Sample_ID_OTU_Table_gDNA') {bnti.name = 'gDNA_bNTI_weighted_rare_15106'} # cDNA_bNTI_weighted_rare_27227 or gDNA_bNTI_weighted_rare_15106
if (bnti.type == 'Sample_ID_OTU_Table_cDNA') {bnti.name = 'cDNA_bNTI_weighted_rare_27227'} # cDNA_bNTI_weighted_rare_27227 or gDNA_bNTI_weighted_rare_15106

#############

# path for output
path.out = "//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/07_Analyses/"

###### load functions

fac.to.char.fun = function(matrix.in) {
  
  i = sapply(matrix.in,is.factor)
  matrix.in[i] = lapply(matrix.in[i],as.character)
  
  return(matrix.in)
  
}

### read in metadata 
meta = fac.to.char.fun(read.csv("//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/Metadata_4Jun2021.csv"))
str(meta)
head(meta)

### read in rate data
rate.in = read.csv("//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/01_Rates/Stegen_EC_DO_Rates.csv")
rate.in$four.pt.rate[rate.in$four.pt.rate < 0] = 0
head(rate.in)

cpi.in = read.csv("//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/01_Rates/Stegen_EC_CPIs.csv")
cpi.in

### merge metadata with rates and CPI
meta.rates = merge(meta,rate.in,by="unique.id",all=T)
meta.rates = merge(meta.rates,cpi.in,by="Cycles",all=T)
dim(meta.rates)
rm('meta')

pdf(paste(path.out,"Fig4a_Rate.Box.CPI.Across.Dry.Days.pdf",sep=""),width = 8,height = 7)
  par(pty="s")
  boxplot(meta.rates$four.pt.rate ~ meta.rates$Daysdry,xlab="Cumulative Days Dry",ylab=expression(Respiration ~ Rate ~ (mg ~ O[2] ~ L^-1 ~ min^-1)),cex.lab=2,cex.axis=1.5); abline(h=c(2,-2),col=2,lwd=2)
  par(new=T,pty="s")
  plot(unique(meta.rates$CPI[order(meta.rates$Daysdry)]) ~ c(1:6),xlim=c(0.5,6.5),ylim=c(0.5,max(meta.rates$CPI)),pch=19,typ="b",axes=F,ylab="",xlab="",col=rgb(red = 34,green = 151,blue = 230,maxColorValue = 255),cex=1.5) # the order of cpi.in must have increasing cycles
  axis(side = 4,at = seq(0.5,0.9,by=0.1),cex.axis=1.5)
  mtext("Control Point Influence", side = 4, line = 3,cex = 2)
  mtext(text = ' a',side = 1,line = -1.25,cex=2,adj = 0)
dev.off()

# paths to bNTI output
amplicon.bnti.path.in = "//pnl/projects/PREMIS/Experimental Mesocosms/Data Package/02_bNTIcalculations/"

# read in bNTI data

if (bnti.type %in% c('Sample_ID_OTU_Table_cDNA','Sample_ID_OTU_Table_gDNA')) {
  
  bnti.in = read.csv(paste(amplicon.bnti.path.in,bnti.name,".csv",sep=""),row.names = 1)
  dim(bnti.in)
  bnti.in[1:5,1:5]
  
  rownames(bnti.in) = gsub(pattern = 'cycles',replacement = " cycles",rownames(bnti.in))
  colnames(bnti.in) = gsub(pattern = 'cycles',replacement = " cycles",colnames(bnti.in))
  rownames(bnti.in) = gsub(pattern = '1 cycles',replacement = "1 cycle",rownames(bnti.in))
  colnames(bnti.in) = gsub(pattern = '1 cycles',replacement = "1 cycle",colnames(bnti.in))
  
}
  
  # bNTI analyses
  meta.rates$median.bnti.treat = NA
  meta.rates$pt.color = NA
  meta.rates$pt.color[which(meta.rates$Cumulative.Treatment == 'Inundated')] = '#8073ac'
  meta.rates$pt.color[which(meta.rates$Cumulative.Treatment == 'Dry')] = '#e08214'
  
  meta.rates$SampleID_bNTI = gsub(pattern = 'cycles',replacement = " cycles",meta.rates$SampleID_bNTI)
  meta.rates$SampleID_bNTI = gsub(pattern = '1 cycles',replacement = "1 cycle",meta.rates$SampleID_bNTI)
  
  for (i in 1:nrow(meta.rates)) {
    
    bnti.temp.vals = c(as.vector(bnti.in[grep(meta.rates$Cycles[i],rownames(bnti.in)),meta.rates$SampleID_bNTI[i]]),as.vector(t(bnti.in[meta.rates$SampleID_bNTI[i],grep(meta.rates$Cycles[i],rownames(bnti.in))])))
    meta.rates$median.bnti.treat[i] = median(bnti.temp.vals,na.rm = T)
    print(c(i,meta.rates$SampleID_bNTI[i],length(bnti.temp.vals)))
    rm('bnti.temp.vals')
    
  }
  head(meta.rates)
  
  pdf(paste(path.out,"Fig3_LnRate.vs.AbsbNTI_",bnti.name,".pdf",sep=""))
    par(pty="s")  
    if (bnti.name == "cDNA_bNTI_weighted_rare_27227" ) { x.lab.temp = expression(paste("Putatively Active Community |",beta,"NTI|"))}  
    if (bnti.name == "gDNA_bNTI_weighted_rare_15106" ) { x.lab.temp = expression(paste("Whole Community |",beta,"NTI|"))}  
    mod.to.plot = log(meta.rates$four.pt.rate + I(min(meta.rates$four.pt.rate[meta.rates$four.pt.rate > 0],na.rm = T)/2)) ~ abs(meta.rates$median.bnti.treat)
    mod = summary(lm(mod.to.plot)); mod
    plot(mod.to.plot,ylab=expression(Ln(Respiration ~ Rate ~ (mg ~ O[2] ~ L^-1 ~ min^-1))),xlab=x.lab.temp,cex.lab=2,cex.axis=1.5,pch=19,cex=1.5,col=meta.rates$pt.color)
    if (round(mod$coefficients[2,4],digits = 6) < 0.05) { abline(mod,lwd=2,col=2) }
    rsq.temp = round(mod$r.squared,digits = 2)
    mtext(text = substitute(paste("R"^2," = ", rsq ," ",sep=""),list(rsq=rsq.temp)),side=3,line = -1.2,las=1,adj=1)
    if (mod$coefficients[2,4] > 0.001) {    mtext(text = paste("p = ",round(mod$coefficients[2,4],digits = 2)," ",sep=""),side=3,line = -2,las=1,adj=1)}
    if (mod$coefficients[2,4] <= 0.001) {    mtext(text = paste("p = ",round(mod$coefficients[2,4],digits = 6)," ",sep=""),side=3,line = -2,las=1,adj=1)} 
    if (bnti.name == "cDNA_bNTI_weighted_rare_27227" ) {mtext(text = ' b',side = 1,line = -1.25,cex=2,adj = 0)}
    if (bnti.name == "gDNA_bNTI_weighted_rare_15106" ) {mtext(text = ' a',side = 1,line = -1.25,cex=2,adj = 0)}
  dev.off()
  

  # box plots
  bnti.for.box = numeric()
  
  for (i in c('0 cycles','1 cycle','2 cycles','3 cycles','4 cycles','5 cycles')) {
    
    bnti.temp = bnti.in[grep(pattern = i,x = rownames(bnti.in)),grep(pattern = i,x = colnames(bnti.in))]
    print(c(i,identical(rownames(bnti.temp),colnames(bnti.temp))))
    print(dim(bnti.temp))
    print(rownames(bnti.temp))
    
    dry.days.temp = unique(meta.rates$Daysdry[grep(pattern = i,x = meta.rates$Cycles)])
    
    bnti.for.box = rbind(bnti.for.box,cbind(as.vector(as.dist(bnti.temp)),i,dry.days.temp))
    
    rm('bnti.temp','dry.days.temp')
    
    
  }
  
  bnti.for.box = fac.to.char.fun(as.data.frame(bnti.for.box))
  colnames(bnti.for.box) = c('bNTI','Cycles',"Daysdry")
  bnti.for.box$bNTI = as.numeric(bnti.for.box$bNTI)
  bnti.for.box$Daysdry = as.numeric(bnti.for.box$Daysdry)
  str(bnti.for.box)
  head(bnti.for.box)
  
  if (bnti.name == "cDNA_bNTI_weighted_rare_27227" ) {
      pdf(paste(path.out,"Fig4b_bNTI.Box.vs.Daysdry_",bnti.name,".pdf",sep=""),width = 8,height = 7)
      y.lab.temp = expression(paste("Putatively Active Community |",beta,"NTI|"))
    
    }
  if (bnti.name == "gDNA_bNTI_weighted_rare_15106" ) {
      pdf(paste(path.out,"FigS3_bNTI.Box.vs.Daysdry_",bnti.name,".pdf",sep=""),width = 8,height = 7)
      y.lab.temp = expression(paste("Whole Community |",beta,"NTI|"))
    
    }
  
    par(pty="s")
    boxplot(bnti.for.box$bNTI ~ bnti.for.box$Daysdry,xlab="Cumulative Days Dry",ylab=y.lab.temp,cex.lab=2,cex.axis=1.5); abline(h=c(2,-2),col=2,lwd=2)
    par(new=T,pty="s")
    plot(unique(meta.rates$CPI[order(meta.rates$Daysdry)]) ~ c(1:6),xlim=c(0.5,6.5),ylim=c(0.5,max(meta.rates$CPI)),pch=19,typ="b",axes=F,ylab="",xlab="",col=rgb(red = 34,green = 151,blue = 230,maxColorValue = 255),cex=1.5) # the order of cpi.in must have increasing cycles
    axis(side = 4,at = seq(0.5,0.9,by=0.1),cex.axis=1.5)
    mtext("Control Point Influence", side = 4, line = 3,cex = 2)
    
    if (bnti.name == "cDNA_bNTI_weighted_rare_27227" ) {mtext(text = ' b',side = 1,line = -1.25,cex=2,adj = 0)}
  
  dev.off()



