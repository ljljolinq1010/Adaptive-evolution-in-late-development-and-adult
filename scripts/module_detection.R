library("ggplot2")
library("matrixStats")
library("RColorBrewer")
library("pheatmap")
library("gplots")
myPalette <- brewer.pal(9, "Set1")

########*fly*########
transcriptome  <- read.table("data/expression/fly_RNAseq_Log2.txt",sep="\t",h=T) 
transcriptome$mean<-rowMeans(transcriptome[-1])
transcriptome$sd<-rowSds(as.matrix(transcriptome[,c(2:28)]))
transcriptome[,c(2:28)]<-apply(transcriptome[,c(2:28)],2,function(x) (x-transcriptome$mean)/transcriptome$sd)
transcriptome<-na.omit(transcriptome)

## PCA
pca<-prcomp((transcriptome[,c(2:28)]))
## atan2
transcriptome<-cbind(transcriptome,pca$x[,1],pca$x[,2])
names(transcriptome)[c(31,32)]<-c("PCA1","PCA2")
transcriptome$atan2<-atan2(transcriptome$PCA1,transcriptome$PCA2)
## find the first gene 
transcriptome<-transcriptome[order(transcriptome$atan2,decreasing = T),]

## heatmap
devName<-c("2h","4h","6h","8h","10h","12h","14h","16h","18h","20h","22h","24h","2d","3d","4d","4.5d","5d",
           "5.5d","6d","6.5d","7d","8d","9d","10d","11d","15d","40d")

transMatrix<-as.matrix(transcriptome[,c(2:28)],ncol=27,nrow=nrow(transcriptome))
pheatmap(transMatrix, cluster_cols=F,cluster_rows=F,scale='none',main="D.melanogaster",show_rownames=F,xlab="Time",
       labels_col = devName,color=colorRampPalette(c("blue","white","red"))(256))

## define modules 
x<-c(1:27)
A=1;K=4;M1=7;B=1;
y_E = K - (K-A) / ( 1 + exp(-B*(x-M1)) ); #early !!!
plot(x,y_E)

A=1;K=4;M1=3;M2=7;B=1;
y1 = A + (K-A) / ( 1 + exp(-B*(x-M1)) ); 
y2 = K - (K-A) / ( 1 + exp(-B*(x-M2)) ); 
y_M=c(y1[x<0.5*(M1+M2)],y2[x>=0.5*(M1+M2)]); #middle!!!
plot(x,y_M)

A=1;K=4;M1=7;M2=11;B=1;
y1 = A + (K-A) / ( 1 + exp(-B*(x-M1)) ); 
y2 = K - (K-A) / ( 1 + exp(-B*(x-M2)) ); 
y_Late=c(y1[x<0.5*(M1+M2)],y2[x>=0.5*(M1+M2)]); #Late !!!
plot(x,y_Late)


A=1;K=4;M1=15;M2=19;B=1;
y1 = A + (K-A) / ( 1 + exp(-B*(x-M1)) ); 
y2 = K - (K-A) / ( 1 + exp(-B*(x-M2)) ); 
y_Larva=c(y1[x<0.5*(M1+M2)],y2[x>=0.5*(M1+M2)]); #Larva !!!
plot(x,y_Larva)

A=1;K=4;M2=19;B=1;
y_A = A + (K-A) / ( 1 + exp(-B*(x-M2)) ); #Adult !!!
plot(x,y_A)

## correlation
transcriptome$corTest_E<-apply(transcriptome[,c(2:28)],1,function(x) cor(x,y_E))
transcriptome$corTest_M<-apply(transcriptome[,c(2:28)],1,function(x) cor(x,y_M))
transcriptome$corTest_Late<-apply(transcriptome[,c(2:28)],1,function(x) cor(x,y_Late))
transcriptome$corTest_Larva<-apply(transcriptome[,c(2:28)],1,function(x) cor(x,y_Larva))
transcriptome$corTest_A<-apply(transcriptome[,c(2:28)],1,function(x) cor(x,y_A))
transcriptome$maxCor<-(colnames(transcriptome[,c(34:38)])[apply(transcriptome[,c(34:38)],1,which.max)])

## plot
transcriptome$module<-rep(6, nrow(transcriptome))
transcriptome$module<-ifelse (((transcriptome$corTest_E > quantile (transcriptome$corTest_E,probs=0.95))&(transcriptome$maxCor=="corTest_E")),1,transcriptome$module)
transcriptome$module<-ifelse (((transcriptome$corTest_M > quantile (transcriptome$corTest_M,probs=0.95))&(transcriptome$maxCor=="corTest_M")),2,transcriptome$module)
transcriptome$module<-ifelse (((transcriptome$corTest_Late > quantile (transcriptome$corTest_Late,probs=0.95))&(transcriptome$maxCor=="corTest_Late")),3,transcriptome$module)
transcriptome$module<-ifelse (((transcriptome$corTest_Larva > quantile (transcriptome$corTest_Larva,probs=0.95))&(transcriptome$maxCor=="corTest_Larva")),4,transcriptome$module)
transcriptome$module<-ifelse (((transcriptome$corTest_A > quantile (transcriptome$corTest_A,probs=0.95))&(transcriptome$maxCor=="corTest_A")),5,transcriptome$module)

ggplot(transcriptome, aes(x=PCA1, y=PCA2,color=as.factor(module))) +geom_point()+
  labs( color = "Modules\n",size=30) +theme_update(plot.title = element_text(hjust = 0.5))+ggtitle("D.melanogaster")+
  scale_color_manual(labels = c("Early embryo","Middle embryo","Late embryo","larva","Pupae.Adult","Non-module"), values=c( myPalette[2], myPalette[1],myPalette[3],myPalette[8],myPalette[4],"darkgrey"))+
  theme(
        axis.text = element_text(size=18,color="black"),
        axis.title = element_text(color="black", size=18, face="bold"),
        plot.title = element_text(hjust = 0.5,size=20,face="bold"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=18))

## for each moduel, we keep genes with correlation coefficient rank in top 5%  (if use 10%, some modules were mixed)
transcriptome_E<-transcriptome[(transcriptome$corTest_E > quantile (transcriptome$corTest_E,probs=0.95))&(transcriptome$maxCor=="corTest_E"),]
transcriptome_M<-transcriptome[(transcriptome$corTest_M> quantile (transcriptome$corTest_M,probs=0.95))&(transcriptome$maxCor=="corTest_M"),]
transcriptome_Late<-transcriptome[(transcriptome$corTest_Late > quantile (transcriptome$corTest_Late,probs=0.95))&(transcriptome$maxCor=="corTest_Late"),]
transcriptome_Larva<-transcriptome[(transcriptome$corTest_Larva > quantile (transcriptome$corTest_Larva,probs=0.95))&(transcriptome$maxCor=="corTest_Larva"),]
transcriptome_A<-transcriptome[(transcriptome$corTest_A > quantile (transcriptome$corTest_A,probs=0.95))&(transcriptome$maxCor=="corTest_A"),]

## modules
max.len = max(nrow(transcriptome_E),nrow(transcriptome_M),nrow(transcriptome_Late),nrow(transcriptome_Larva),nrow(transcriptome_A))

gene_E<-data.frame(transcriptome_E$Ensembl.Gene.ID)
gene_E[c((nrow(gene_E)+1):max.len),]<-NA

gene_M<-data.frame(transcriptome_M$Ensembl.Gene.ID)
gene_M[c((nrow(gene_M)+1):max.len),]<-NA

gene_Late<-data.frame(transcriptome_Late$Ensembl.Gene.ID)
gene_Late[c((nrow(gene_Late)+1):max.len),]<-NA

gene_Larva<-data.frame(transcriptome_Larva$Ensembl.Gene.ID)
gene_Larva[c((nrow(gene_Larva)+1):max.len),]<-NA



moduleGenes<-data.frame(gene_E$transcriptome_E.Ensembl.Gene.ID,gene_M$transcriptome_M.Ensembl.Gene.ID,
                        gene_Late$transcriptome_Late.Ensembl.Gene.ID,gene_Larva$transcriptome_Larva.Ensembl.Gene.ID,
                        transcriptome_A$Ensembl.Gene.ID)

names(moduleGenes)<-c("Early embryo", "Middle embryo","Late embryo","Larva","Pupae.Adult")

########*mouse*########
transcriptome  <- read.table("data/expression/mouse_RNAseq_Log2.txt",sep="\t",h=T) 
transcriptome$mean<-rowMeans(transcriptome[-1])
transcriptome$sd<-rowSds(as.matrix(transcriptome[,c(2:18)]))
transcriptome[,c(2:18)]<-apply(transcriptome[,c(2:18)],2,function(x) (x-transcriptome$mean)/transcriptome$sd)
transcriptome<-na.omit(transcriptome)
## PCA
pca<-prcomp((transcriptome[,c(2:18)]))
## anta2
transcriptome<-cbind(transcriptome,pca$x[,1],pca$x[,2])
names(transcriptome)[c(21,22)]<-c("PCA1","PCA2")
transcriptome$atan2<-atan2(transcriptome$PCA1,transcriptome$PCA2)
## find the first gene
transcriptome<-transcriptome[order(transcriptome$atan2,decreasing = T),]

## heatmap
devName<-c("0.5d","1.5d","2d","3.5d","7.5d","8.5d","9d","9.5d","10.5d","11.5d","12.5d","13.5d","14.5d","15.5d","16.5d","17.5d","18.5d")
transMatrix <- as.matrix(transcriptome[,c(2:18)],ncol=17,nrow=nrow(transcriptome))
pheatmap(transMatrix, cluster_cols=F,cluster_rows=F,scale='none',main="M.musculus",show_rownames=F,xlab="Time",
         labels_col = devName,color=colorRampPalette(c("blue","white","red"))(256))

## define modules 
x<-c(1:17)
A=1;K=8;M1=1;B=1;
y_E = K - (K-A) / ( 1 + exp(-B*(x-M1)) ); #early !!!
plot(x,y_E)
A=1;K=4;M1=5;M2=11;B=1; 
y1 = A + (K-A) / ( 1 + exp(-B*(x-M1)) ); 
y2 = K - (K-A) / ( 1 + exp(-B*(x-M2)) ); 
y_M=c(y1[x<0.5*(M1+M2)],y2[x>=0.5*(M1+M2)]); #middle !!!
plot(x,y_M)

A=1;K=4;M2=14;B=1;
y_L = A + (K-A) / ( 1 + exp(-B*(x-M2)) ); #late !!!
plot(x,y_L)


## correlation
transcriptome$corTest_E<-apply(transcriptome[,c(2:18)],1,function(x) cor(x,y_E))
transcriptome$corTest_M<-apply(transcriptome[,c(2:18)],1,function(x) cor(x,y_M))
transcriptome$corTest_L<-apply(transcriptome[,c(2:18)],1,function(x) cor(x,y_L))
transcriptome$maxCor<-(colnames(transcriptome[,c(24:26)])[apply(transcriptome[,c(24:26)],1,which.max)])
## plot
transcriptome$module<-rep(4, nrow(transcriptome))
transcriptome$module<-ifelse (((transcriptome$corTest_E > quantile (transcriptome$corTest_E,probs=0.9))&(transcriptome$maxCor=="corTest_E")),1,transcriptome$module)
transcriptome$module<-ifelse (((transcriptome$corTest_M > quantile (transcriptome$corTest_M,probs=0.9))&(transcriptome$maxCor=="corTest_M")),2,transcriptome$module)
transcriptome$module<-ifelse (((transcriptome$corTest_L > quantile (transcriptome$corTest_L,probs=0.9))&(transcriptome$maxCor=="corTest_L")),3,transcriptome$module)

ggplot(transcriptome, aes(x=PCA1, y=PCA2,color=as.factor(module))) +geom_point()+
  labs( color = "Modules\n",size=30) +theme_update(plot.title = element_text(hjust = 0.5))+ggtitle("M.musculus")+
  scale_color_manual(labels = c("Early embryo","Middle embryo","Late embryo","Non-module"), values=c( myPalette[2], myPalette[1],myPalette[3],"darkgrey"))+
  theme(
    axis.text = element_text(size=18,color="black"),
    axis.title = element_text(color="black", size=18, face="bold"),
    plot.title = element_text(hjust = 0.5,size=20,face="bold"),
    legend.text=element_text(size=15),
    legend.title=element_text(size=18))

## modules
transcriptome_E<-transcriptome[(transcriptome$corTest_E > quantile (transcriptome$corTest_E,probs=0.9))&(transcriptome$maxCor=="corTest_E"),]
transcriptome_M<-transcriptome[(transcriptome$corTest_M > quantile (transcriptome$corTest_M,probs=0.9))&(transcriptome$maxCor=="corTest_M"),]
transcriptome_L<-transcriptome[(transcriptome$corTest_L > quantile (transcriptome$corTest_L,probs=0.9))&(transcriptome$maxCor=="corTest_L"),]

moduleGenes<-data.frame(transcriptome_E$Ensembl.Gene.ID,transcriptome_M$Ensembl.Gene.ID,transcriptome_L$Ensembl.Gene.ID)
names(moduleGenes)<-c("Early embryo", "Middle embryo","Late embryo")

























