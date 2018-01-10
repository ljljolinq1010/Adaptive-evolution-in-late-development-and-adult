library("ggplot2")
library("reshape2")
library("RColorBrewer")
######* module analysis *#####
 moduletest<-function(module,selectome,orgName) {
   if(orgName=="D.melanogaster") {
     earlyM<-selectome[selectome$Ensembl.Gene.ID%in%module$Early.embryo,]
     middleM<-selectome[selectome$Ensembl.Gene.ID%in%module$Middle.embryo,]
     lateM<-selectome[selectome$Ensembl.Gene.ID%in%module$Late.embryo,]
     larvaM<-selectome[selectome$Ensembl.Gene.ID%in%module$Larva,]
     pupAdultM<-selectome[selectome$Ensembl.Gene.ID%in%module$Pupae.Adult,]
     
     moduleList<-list(earlyM,middleM,lateM,larvaM,pupAdultM)
     moduleDF<-rbind(earlyM,middleM,lateM,larvaM,pupAdultM)
     moduleNames<-names(module)
     legendName<-"Melanogaster group branch"
     
   }
   if (orgName=="D.rerio") {
     clevBlasM<-selectome[selectome$Ensembl.Gene.ID%in%module$Cleav.Blastula,]
     gastrulaM<-selectome[selectome$Ensembl.Gene.ID%in%module$Gastrula,]
     segmentationM<-selectome[selectome$Ensembl.Gene.ID%in%module$Segmentation,]
     pharyngulaM<-selectome[selectome$Ensembl.Gene.ID%in%module$Pharyngula,]
     larvaM<-selectome[selectome$Ensembl.Gene.ID%in%module$Larva,]
     juvenileM<-selectome[selectome$Ensembl.Gene.ID%in%module$Juvenile,]
     adultM<-selectome[selectome$Ensembl.Gene.ID%in%module$Adult,]
     
     moduleList<-list(clevBlasM,gastrulaM,segmentationM,pharyngulaM,larvaM,juvenileM,adultM)
     moduleDF<-rbind(clevBlasM,gastrulaM,segmentationM,pharyngulaM,larvaM,juvenileM,adultM)
     moduleNames<-names(module)
     legendName<-"Clupeocephala branch"
     
   }
   if (orgName=="M.musculus") {
     earlyM<-selectome[selectome$Ensembl.Gene.ID%in%module$Early.embryo,]
     middleM<-selectome[selectome$Ensembl.Gene.ID%in%module$Middle.embryo,]
     lateM<-selectome[selectome$Ensembl.Gene.ID%in%module$Late.embryo,]
     
     moduleList<-list(earlyM,middleM,lateM)
     moduleDF<-rbind(earlyM,middleM,lateM)
     moduleNames<-names(module)
     legendName<-"Murinae branch"
     
   }
   
   ## retrieve lrt  
   # weak evidence of positive selection 
   lrt_weak <- lapply(moduleList, function(x) sqrt(sqrt(x$lrt[x$lrt>0])))
   # strong evidence of positive selection 
   lrt_strong <- lapply(moduleList, function(x) sqrt(sqrt(x$lrt[x$lrt>0&x$qvalue<0.2])))
   
   #####* test proportion of positive selection *#####
   ## proportion
   allGeneNum<-unlist(lapply(moduleList, nrow))
   lrt_weak_GeneNum<-unlist(lapply(lrt_weak, length))
   lrt_strong_GeneNum<-unlist(lapply(lrt_strong, length))
   lrt_weak_prop<-lrt_weak_GeneNum/allGeneNum
   lrt_strong_prop<-lrt_strong_GeneNum/allGeneNum
   ## Chi-square Test of Goodness-of-Fit
   expected<-allGeneNum/sum(allGeneNum)
   observed_weak<-lrt_weak_GeneNum
   observed_strong<-lrt_strong_GeneNum
   chq_weak<-chisq.test(x = observed_weak,p = expected)
   chq_strong<-chisq.test(x = observed_strong,p = expected)
   ## plot
   if (orgName=="D.melanogaster") {
     ylim_lrt_weak=c(0,1)
     ylim_lrt_strong=c(0,0.1)
     text_y<-c(-0.003,0.003)
     ylim_boxplot<-c(-0.2,3)
     boxColor<-c("grey","grey","grey","grey","red")
   }
   if (orgName=="D.rerio") {
     ylim_lrt_weak<-c(0,1)
     ylim_lrt_strong<-c(0,0.4)
     text_y<-c(-0.015,0.015)
     ylim_boxplot<-c(-0.2,3.5)
     boxColor<-c("grey","grey","grey","blue","grey","red","red")
   }
   if (orgName=="M.musculus") {
     ylim_lrt_weak=c(0,1)
     ylim_lrt_strong=c(0,0.02)
     text_y<-c(-0.0006,0.0006)
     ylim_boxplot<-c(-0.2,2.5)
     boxColor<-c("grey","blue","red")
   }
   # proportion of weak lrt
   prop<-barplot(lrt_weak_prop,ylim=ylim_lrt_weak,col=3,ylab=expression(paste("Proportion of genes with", " ", paste(Delta,"lnL"), " > 0")), main=legendName)
   legend("topleft",legend=paste0("p=",signif((chq_weak$p.value),2)), bty = 'n',cex = 1)
   text(x =prop, y = -0.04, srt = 45,cex.lab=1, adj = 1,  labels = moduleNames, xpd = TRUE)
   text(x=prop, y=0.04, labels=paste0("n=", allGeneNum))
   # proportion of strong lrt
   prop<-barplot(lrt_strong_prop,ylim=ylim_lrt_strong,col=3,ylab=expression(paste("Proportion of genes with", " ", paste(Delta,"lnL"), " > 0")), main=legendName)
   legend("topleft",legend=paste0("p=",signif((chq_strong$p.value),2)), bty = 'n',cex = 1)
   text(x =prop, y = text_y[1], srt = 45,cex.lab=1, adj = 1,  labels = moduleNames, xpd = TRUE)
   text(x=prop, y=text_y[2], labels=paste0("n=", allGeneNum))

   #####* test strength of positive selection *#####
   ## significant test for mean value 
   meanDistri<-c()
   p_temp<-c()
   pValue<-c()
   for (i in 1:length(lrt_weak)) {
     randomModule<-replicate(10000,sample(unlist(lrt_weak),length(lrt_weak[[i]]),replace=F))
     meanDistri<-apply(randomModule,2,mean)
     meanLrt<-mean(meanDistri)
     sdLrt<-sd(meanDistri)
     if (meanLrt > mean(lrt_weak[[i]]) ) {
       p_temp<-pnorm(mean(lrt_weak[[i]]),meanLrt,sdLrt,lower.tail = T)
       } else {
       p_temp<-pnorm(mean(lrt_weak[[i]]),meanLrt,sdLrt,lower.tail = F)
     }
     pValue<-cbind(pValue,p_temp)
   }
   ## ggplot 
   moduleDF<-subset(moduleDF,lrt>0)
   moduleDF$lrt<-sqrt(sqrt(moduleDF$lrt))
   moduleDF$moduleNames<-rep(moduleNames,unlist(lapply(lrt_weak, length))) 
   moduleDF$moduleNames<-as.character(moduleDF$moduleNames)
   moduleDF$moduleNames<-factor(moduleDF$moduleNames,levels=unique(moduleDF$moduleNames))
   
   if(orgName=="D.melanogaster") {
     color<-c("blue","black","black","black","red")
     y1=-0.5
     y2=3
   }
   if (orgName=="M.musculus") {
     color<-c("black","blue","red")
     y1=-0.25
     y2=2.5
   }
   if (orgName=="D.rerio") {
     color<-c("black","black","black","blue","black","red","black")
     y1=-0.5
     y2=3.5
   }
 
   print(
     ggplot(data=moduleDF,aes(x=moduleNames,y=lrt,color=moduleNames))+geom_jitter(position=position_jitter(0.1),pch=16,cex=2,alpha = 0.5) +
       ylim(y1,y2)+stat_summary(fun.y = mean,fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5,color="black")+
       scale_color_manual(values=color) +
       labs(title=legendName,x="", y = expression(paste("Fourth root"," ",(paste(Delta,"lnL")))),size=20)+
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(),axis.line = element_line(colour = "black"),
             axis.text.y = element_text(size=18,color="black"),
             axis.title.y = element_text(color="black", size=18, face="bold"),
             legend.position="none",axis.text.x=element_text(angle = 45, hjust = 1,size=18), plot.title = element_text(hjust = 0.5,size=20,face="bold"))+
       geom_hline(yintercept = mean(moduleDF$lrt),colour="green",linetype=2,size=1)
   )

 } 
 
## fly
selectome <- read.table("data/selectome/fly_melanogaster_group_selectome.txt",sep="\t",h=T) 
module<-read.table("results/modules/fly_module_genes.txt",sep="\t",h=T) 
moduletest(module,selectome,"D.melanogaster")
 
## zf
selectome <- read.table("data/selectome/zf_clupeocephala_selectome.txt",sep="\t",h=T) 
module<-read.table("results/modules/zf_module_genes.txt",sep="\t",h=T) 
moduletest(module,selectome,"D.rerio")
 
## mouse
selectome <- read.table("data/selectome/mouse_murine_selectome.txt",sep="\t",h=T) 
module<-read.table("results/modules/mouse_module_genes.txt",sep="\t",h=T) 
moduletest(module,selectome,"M.musculus")
 

####* transcriptome index analysis *####
library(RColorBrewer)
myPalette <- brewer.pal(9, "Set1")

transIndex(transcriptome,selectome,orgName)  {
  if (orgName=="D.melanogaster") {
    timePoint<-c(2:28) 
    devTime1<-log2(c(seq(2,24,by=2),44,70,96,108,120,132,144,156,168,192,216,240,264,360,960 ))
    
    devTime2<-c("2h","4h","6h","8h","10h","12h","14h","","18h","","22h","","2d","3d","4d","","5d",
                          "","6d","","7d","8d","9d","","11d","15d","40d")
    devTimeColor<-c(rep(myPalette[2],3),rep(myPalette[1],2),rep(myPalette[3],19),rep(myPalette[4],3))
    modules = list(Early.development = 1:3, Middle.devlopment = 4:5, Late.development = 6:24,Adult=25:27)
    legendName<-"Melanogaster group branch"
    ylimBoxplot<-c( 0.64, 0.76)
    ytext=0.631
    ylimRegression<-c(0.67,0.73)
  }
  if (orgName=="M.musculus") {
    timePoint<-c(2:18) 
    devTime1<-c(0.5,1.5,2,3.5,7.5,8.5,9,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5)
    
    devTime2<-c("0.5d","1.5d","2d","3.5d","7.5d","8.5d","9d","9.5d","10.5d","11.5d","12.5d","13.5d","14.5d","15.5d","16.5d","17.5d","18.5d")
    devTimeColor<-c(myPalette[5],rep(myPalette[2],4),rep(myPalette[1],6),rep(myPalette[3],6))
    modules = list(Maternal.stage=1,Early.development = 2:5, Middle.devlopment = 6:11, Late.development = 12:17)
    legendName<-"Murinae branch"
    ylimBoxplot<-c(0.25,0.33)
    ylimRegression<-c(0.273,0.295)
    ytext<-0.2435
  }
  if (orgName=="D.rerio") {
    timePoint<-c(3:61) ## remove unfertilized egg
    devTime1<-log2(c(0.25,0.75,1.25,1.75,2.25,2.75,3.33,4,4.66,5.33,6,7,8,9,10,10.33,11,11.66,
                       12,13,14,15,16,17,18,19,20,21,22,23,25,27,30,34,38,42,48,60,72,96,144,192,240,
                       336,432,720,960,1080,1320,1560,1920,2160,2520,2880,5040,6480,10080,12960,15120))
    devTime2<-c("0.25h","0.75h","1.25h","1.75h","2.25h","","3.5h","","","5.5h","","","8h","","","","","11h",
                "","","","","","","18h","","","","","","25h","","","","38h","","2d","","3d","4d","6d","","10d",
                "","18d","30d","40d","","55d","","80d","","3.5m","","7m","9m","1y2m","","1y9m")
    devTimeColor<-c(rep(myPalette[5],4),rep(myPalette[2],11),rep(myPalette[1],21),rep(myPalette[3],16),rep(myPalette[4],7))
    modules = list(Maternal.stage=1:4,Early.development = 5:15, Middle.devlopment = 16:36, Late.development = 37:52,Adult=53:59)
    legendName<-"Clupeocephala branch"
    ylimBoxplot<-c(1.03,1.16)
    ylimRegression<-c(1.075,1.10)
    ytext=1.02
  }
  ## calculate index
  transSelec<- merge(transcriptome, selectome, by="Ensembl.Gene.ID")
  transSelec$lrt <- ifelse(transSelec$lrt<0,0, transSelec$lrt)
  transSelec$lrt <-sqrt(sqrt(transSelec$lrt))
  transIndex<-c()
  transIndex<-apply(transSelec[timePoint], 2, function(x) sum(x*(transSelec[,"lrt"]))/sum(x))
  ## bootstrap analysis
  cat("\nbootstrap analysis...")
  bootGeneID<-replicate(10000,sample(transSelec$Ensembl.Gene.ID,replace=T))
  transIndexBoot<-c()
  transIndexBoot1<-c()
  for (i in 1:10000) {
    tempID<-data.frame(bootGeneID[,i])
    names(tempID)<-"Ensembl.Gene.ID"
    tempTransSelec<-merge(tempID,transSelec,by="Ensembl.Gene.ID")
    
    transIndexBoot1<-apply(tempTransSelec[timePoint],2, function(x) sum(x*(tempTransSelec[,"lrt"]))/sum(x))
    transIndexBoot<-rbind(transIndexBoot,transIndexBoot1)
  }
  ## calculate  mean index of each module, and compare the mean with wilcox test.
  meanIndex<-c()
  meanIndex1<-c()
  for (i in 1:10000) {
    meanIndex1 <- lapply( modules,function(x) mean(transIndexBoot[i,][x]) )
    meanIndex <- rbind(meanIndex,meanIndex1)
  }
  ## pairwise test
  pt<-c()
  pt$meanIndex<-unlist(meanIndex)
  pt<-data.frame(pt)
  pt$group=rep(c(1:length(modules)),each=10000)
  pairwiseWT<-pairwise.wilcox.test( pt$meanIndex,pt$group,p.adjust.method = "fdr")
  pValue<-c()
  if (orgName=="D.rerio") {
    pValue$one<-c(1,pairwiseWT$p.value[,1])
    pValue$second<-c(pairwiseWT$p.value[,1][1],1,pairwiseWT$p.value[,2][c(2:4)])
    pValue$third<-c(pairwiseWT$p.value[,1][2],pairwiseWT$p.value[,2][2],1,pairwiseWT$p.value[,3][c(3:4)])
    pValue$forth<-c(pairwiseWT$p.value[,1][3],pairwiseWT$p.value[,2][3],pairwiseWT$p.value[,3][3],1,pairwiseWT$p.value[,4][4])
    pValue$fifth<-c(pairwiseWT$p.value[4,],1)
    pValue<-data.frame(pValue)
    colnames(pValue)<-names(modules)
    rownames(pValue)<-names(modules)
   
  } else {
    pValue$one<-c(1,pairwiseWT$p.value[,1])
    pValue$second<-c(pairwiseWT$p.value[,1][1],1,pairwiseWT$p.value[,2][c(2:3)])
    pValue$third<-c(pairwiseWT$p.value[,1][2],pairwiseWT$p.value[,2][2],1,pairwiseWT$p.value[,3][3])
    pValue$forth<-c(pairwiseWT$p.value[3,],1)
    pValue<-data.frame(pValue)
    colnames(pValue)<-names(modules)
    rownames(pValue)<-names(modules)
  }
  

  ## boxplot
  if (orgName=="D.rerio"){
    boxplotData<-matrix(unlist(meanIndex),ncol = 5,nrow = 10000)
    
  } else {
    boxplotData<-matrix(unlist(meanIndex),ncol = 4,nrow = 10000)
    
  }
  boxplot(boxplotData,las=2,ylim=ylimBoxplot,pch=16,outcex=0.5,boxwex=0.7, xaxt = "n",main=legendName,cex.lab=1.2,cex.main=1.2,
          col=c(unique(devTimeColor)),ylab =expression(paste("Transcriptoem index of fourth root"," ",(paste(Delta,"lnL")))))
  text(x =  seq_along(names(modules)), y = ytext, srt = 45,cex.lab=1.2,adj = 1,  labels = names(modules), xpd = TRUE)
  ## regression plot
  if (orgName=="D.rerio") {
    lmd4 <- lm(unlist(transIndex) ~ poly(devTime1, 4, raw=TRUE))
    a<-summary(lmd4)$coef[,1][[1]]
    b<-summary(lmd4)$coef[,1][[2]]
    c<-summary(lmd4)$coef[,1][[3]]
    d<-summary(lmd4)$coef[,1][[4]]
    e<-summary(lmd4)$coef[,1][[5]]
    r2<-signif(summary(lmd4)$adj.r.squared, 2)
    f<-summary(lmd4)$fstatistic
    polyModel <- function(x) { eval(a) + eval(b)*x + eval(c)*x^2+ eval(d)*x^3+ eval(e)*x^4}
  } else {
    lmd2 <- lm(unlist(transIndex) ~ poly(devTime1, 2, raw=TRUE))
    a<-summary(lmd2)$coef[,1][[1]]
    b<-summary(lmd2)$coef[,1][[2]]
    c<-summary(lmd2)$coef[,1][[3]]
    r2<-signif(summary(lmd2)$adj.r.squared, 2)
    f<-summary(lmd2)$fstatistic
    polyModel <- function(x) { eval(a) + eval(b)*x + eval(c)*x^2}
  }
  
  curve(polyModel, min(devTime1), max(devTime1), col="black",xlab="Time",
        ylab=expression(paste("Transcriptoem index of fourth root"," ",(paste(Delta,"lnL")))),
        ylim=ylimRegression, main=legendName,xaxt="n",lwd=6,lty=1,cex.lab=1.2,cex.axis=1.2,cex.main=1.2)
  points(devTime1, unlist(transIndex), pch=16, lwd=6)
  for (j in 1:length(timePoint)) {
    axis(side=1, at=devTime1[j], col.axis=devTimeColor[j], labels=devTime2[j], las=2,cex.axis=1.2) # Add development stages as labels, each color represents one meta development stage 
  } 
  myP<-signif(pf(f[1],f[2],f[3],lower.tail=F), 2)
  rp = vector('expression',2)
  rp[1] = substitute(expression(R^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(p == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(myP, digits = 2)))[2]
  legend("topleft",legend=rp, bty = 'n',cex = 1.2,col=c("black","white"),lty=c(1,1),lwd=c(2,2))

}
   
## fly
selectome <- read.table("data/selectome/fly_melanogaster_group_selectome.txt",sep="\t",h=T) 
transcriptome<-read.table("data/expression/fly_RNAseq_Log2.txt",sep="\t",h=T) 
moduletest(transcriptome,selectome,"D.melanogaster")

## zf
selectome <- read.table("data/selectome/zf_clupeocephala_selectome.txt",sep="\t",h=T) 
transcriptome<-read.table("data/expression/zf_Microarray_Log2.txt",sep="\t",h=T) 
moduletest(transcriptome,selectome,"D.rerio")

## mouse
selectome <- read.table("data/selectome/mouse_murine_selectome.txt",sep="\t",h=T) 
transcriptome<-read.table("data/expression/mouse_RNAseq_Log2.txt",sep="\t",h=T) 
moduletest(transcriptome,selectome,"M.musculus")

#####* pathway analysis *#####
pathwaySelAna<-function(orgName,pathwaySel,pathwayID_geneID,ensemblID_entrezID,transcriptome) {
  ## retrieve positive selected pathways
  positiveSelID<-subset(pathwaySel,setQ<0.2)$setID.orig
  ## retrieve genes of positive selected pathways
  positiveGeneID<-lapply(positiveSelID, function(x) subset(pathwayID_geneID,setID==x)$geneID)
  positiveGeneID<-data.frame(unique(unlist(positiveGeneID)))
  names(positiveGeneID)<-"geneID"
  positiveGeneID<-merge(positiveGeneID,ensemblID_entrezID,by="geneID")
  
  ## retrieve genes of non-positive selected pathways
  nonPositiveGeneID<-lapply(positiveSelID, function(x) subset(pathwayID_geneID,setID!=x)$geneID)
  nonPositiveGeneID<-data.frame(unique(unlist(nonPositiveGeneID)))
  names(nonPositiveGeneID)<-"geneID"
  nonPositiveGeneID<-merge(nonPositiveGeneID,ensemblID_entrezID,by="geneID")
  nonPosGeneID<-nonPositiveGeneID

  ## retrieve gene expression
  posExpression<-merge(positiveGeneID,transcriptome, by="Ensembl.Gene.ID")
  medianPosExp<-apply(posExpression[3:length(posExpression)],2, function(x) median(x))
  nonPosExpression<-merge(nonPosGeneID,transcriptome, by="Ensembl.Gene.ID")
  medianNonPosExp<-apply(nonPosExpression[3:length(nonPosExpression)],2, function(x) median(x))
  ratio<-medianPosExp/medianNonPosExp
 
  if (orgName=="D.melanogaster") {
    devTime1<-log2(c(seq(2,24,by=2),44,70,96,108,120,132,144,156,168,192,216,240,264,360,960 ))
    
    devTime2<-c("2h","4h","6h","8h","10h","12h","14h","","18h","","22h","","2d","3d","4d","","5d",
                "","6d","","7d","8d","9d","","11d","15d","40d")
    devTimeColor<-c(rep(myPalette[2],3),rep(myPalette[1],2),rep(myPalette[3],19),rep(myPalette[4],3))
    legendName<-"Melanogaster group branch"
    ylimRegression<-c(0.8,1.2)
  }
  if (orgName=="M.musculus") {
    devTime1<-c(1.5,2,3.5,7.5,8.5,9,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5)
    devTime2<-c("1.5d","2d","3.5d","7.5d","8.5d","9d","9.5d","10.5d","11.5d","12.5d","13.5d","14.5d","15.5d","16.5d","17.5d","18.5d")
    devTimeColor<-c(rep(myPalette[2],4),rep(myPalette[1],6),rep(myPalette[3],6))
    legendName<-"Murinae branch"
    ylimRegression<-c(0.1,0.9)
  }
  if (orgName=="D.rerio") {
    devTime1<-log2(c(0.25,0.75,1.25,1.75,2.25,2.75,3.33,4,4.66,5.33,6,7,8,9,10,10.33,11,11.66,
                     12,13,14,15,16,17,18,19,20,21,22,23,25,27,30,34,38,42,48,60,72,96,144,192,240,
                     336,432,720,960,1080,1320,1560,1920,2160,2520,2880,5040,6480,10080,12960,15120))
    devTime2<-c("0.25h","0.75h","1.25h","1.75h","2.25h","","3.5h","","","5.5h","","","8h","","","","","11h",
                "","","","","","","18h","","","","","","25h","","","","38h","","2d","","3d","4d","6d","","10d",
                "","18d","30d","40d","","55d","","80d","","3.5m","","7m","9m","1y2m","","1y9m")
    devTimeColor<-c(rep(myPalette[5],4),rep(myPalette[2],11),rep(myPalette[1],21),rep(myPalette[3],16),rep(myPalette[4],7))
    legendName<-"Clupeocephala branch"
    ylimRegression<-c(0.65,1.1)
  }
  ## regression plot (polynomial model start with degree 3, progressively increase degree until the improvement is not significant (anova test > 0.05))
  if (orgName=="M.musculus") {
    lmd5 <- lm(unlist(ratio) ~ poly(devTime1, 5, raw=TRUE))
    a<-summary(lmd5)$coef[,1][[1]]
    b<-summary(lmd5)$coef[,1][[2]]
    c<-summary(lmd5)$coef[,1][[3]]
    d<-summary(lmd5)$coef[,1][[4]]
    e<-summary(lmd5)$coef[,1][[5]]
    h<-summary(lmd5)$coef[,1][[6]]
    
    r2<-signif(summary(lmd5)$adj.r.squared, 2)
    f<-summary(lmd5)$fstatistic
    polyModel <- function(x) { eval(a) + eval(b)*x + eval(c)*x^2+ eval(d)*x^3+ eval(e)*x^4+ eval(h)*x^5}
  } else {
    lmd4 <- lm(unlist(ratio) ~ poly(devTime1, 4, raw=TRUE))
    a<-summary(lmd4)$coef[,1][[1]]
    b<-summary(lmd4)$coef[,1][[2]]
    c<-summary(lmd4)$coef[,1][[3]]
    d<-summary(lmd4)$coef[,1][[4]]
    e<-summary(lmd4)$coef[,1][[5]]
    r2<-signif(summary(lmd4)$adj.r.squared, 2)
    f<-summary(lmd4)$fstatistic
    polyModel <- function(x) { eval(a) + eval(b)*x + eval(c)*x^2+ eval(d)*x^3+ eval(e)*x^4}
  }

  curve(polyModel, min(devTime1), max(devTime1), col="cyan3",xlab="Time",
        ylab="Ratio of median expression",
        ylim=ylimRegression, main=legendName,xaxt="n",lwd=6,lty=1,cex.lab=1.2,cex.axis=1.2,cex.main=1.2)
  points(devTime1, unlist(ratio), pch=16, lwd=6,col="cyan3")
  for (j in 1:length(devTime1)) {
    axis(side=1, at=devTime1[j], col.axis=devTimeColor[j], labels=devTime2[j], las=2,cex.axis=1.2) # Add development stages as labels, each color represents one meta development stage 
  } 
  myP<-signif(pf(f[1],f[2],f[3],lower.tail=F), 2)
  rp = vector('expression',2)
  rp[1] = substitute(expression(R^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(p == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(myP, digits = 2)))[2]
  legend("topleft",legend=rp, bty = 'n',cex = 1.2,col=c("black","white"),lty=c(1,1),lwd=c(2,2))
}


## zf
pathwaySel<-read.table("data/pathway_selection/zf_clupeocephala_setscores_postpruning.txt",sep="\t",h=T,quote = "")
pathwayID_geneID<-read.table("data/pathway_selection/zf_pathwayID_geneID.txt",sep="\t",h=F,quote = "")
names(pathwayID_geneID)<-c("setID","setName","geneID")
ensemblID_entrezID<-read.table("data/pathway_selection/zf_ensembl_entrez_ID_one2one.txt",sep="\t",h=F)
names(ensemblID_entrezID)<-c("Ensembl.Gene.ID","geneID")
transcriptome<-read.table("data/expression/zf_Microarray_Log2.txt",sep="\t",h=T)
transcriptome$egg.0min<-NULL
pathwaySelAna("D.rerio",pathwaySel,pathwayID_geneID,ensemblID_entrezID,transcriptome)

## mouse
pathwaySel<-read.table("data/pathway_selection/mouse_murinae_setscores_postpruning.txt",sep="\t",h=T,quote = "")
pathwayID_geneID<-read.table("data/pathway_selection/mouse_pathwayID_geneID",sep="\t",h=F,quote = "")
names(pathwayID_geneID)<-c("setID","setName","geneID")
ensemblID_entrezID<-read.table("data/pathway_selection/mouse_ensembl_entrez_ID_one2one.txt",sep="\t",h=F)
names(ensemblID_entrezID)<-c("Ensembl.Gene.ID","geneID")
transcriptome<-read.table("data/expression/mouse_RNAseq_Log2.txt",sep="\t",h=T)
transcriptome$Mean.2cell<-NULL ## because the median expression of this time point is 0, we need to remove it
pathwaySelAna("M.musculus",pathwaySel,pathwayID_geneID,ensemblID_entrezID,transcriptome)

## fly
pathwaySel<-read.table("data/pathway_selection/fly_setscores_postpruning.txt",sep="\t",h=T,quote = "")
pathwayID_geneID<-read.table("data/pathway_selection/fly_pathwayID_geneID",sep="\t",h=F,quote = "")
names(pathwayID_geneID)<-c("setID","setName","geneID")
ensemblID_entrezID<-read.table("data/pathway_selection/fly_ensembl_entrez_ID_one2one.txt",sep="\t",h=F)
names(ensemblID_entrezID)<-c("Ensembl.Gene.ID","geneID")
transcriptome<-read.table("data/expression/fly_RNAseq_Log2.txt",sep="\t",h=T)
pathwaySelAna("D.melanogaster",pathwaySel,pathwayID_geneID,ensemblID_entrezID,transcriptome)






















 