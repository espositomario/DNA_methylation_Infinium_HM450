rm(list=ls())
library(minfi)
setwd('/Users/marius/prj/unibo/DRD/R/prj/')
baseDir <- ("Input/")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)

Red <- data.frame(getRed(RGset))
head(Red)
Green <- data.frame(getGreen(RGset))
head(Green)

probe_red <- Red[rownames(Red)=="10737353",]
probe_green <- Green[rownames(Green)=="10737353",]
load('Illumina450Manifest_clean.RData')
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="10737353", 'Infinium_Design_Type']
address_df=data.frame(Sample=colnames(probe_green), Red=unlist(probe_red, use.names = FALSE), Green=unlist(probe_green, use.names = FALSE),Type='II',Color='both')
address_df

MSet.raw <- preprocessRaw(RGset)
dim(MSet.raw)

qc <- getQC(MSet.raw)
plotQC(qc)

controlStripPlot(RGset, controls="NEGATIVE")


detP <- detectionP(RGset) 
failed <- detP>0.05
n_failed=colSums(failed)
summary_df <- data.frame(Sample=colnames(failed),Num_failed_position=n_failed,percentage=colMeans(failed)*100,row.names = NULL)
summary_df


wt <- targets[targets$Group=="WT", "Basename"]
mt <- targets[targets$Group=="MUT", "Basename"]
wt <- gsub(baseDir, "", wt) 
mt <- gsub(baseDir, "", mt) 

wt_subset <- MSet.raw[,colnames(MSet.raw) %in% wt]
mt_subset <- MSet.raw[,colnames(MSet.raw) %in% mt]

wtBeta <- getBeta(wt_subset)
wtM <- getM(wt_subset)
mtBeta <- getBeta(mt_subset)
mtM <- getM(mt_subset)

density_wtBeta <- density(apply(wtBeta,MARGIN=1,mean,na.rm=T),na.rm=T)
density_mtBeta <- density(apply(mtBeta,MARGIN=1,mean,na.rm=T),na.rm=T)
density_wtM <- density(apply(wtM,MARGIN=1,mean,na.rm=T),na.rm=T)
density_mtM <- density(apply(mtM,MARGIN=1,mean,na.rm=T),na.rm=T)


par(mfrow=c(1,2))
#beta density MT vs WT

plot(density_wtBeta,main="Density of Beta Values",col="#03B2C9",lwd=2.5,xlab='beta')
lines(density_mtBeta,main="Density of Beta Values",col="#DE6E4B",lwd=2.5)
legend('topright', legend=c("WT","MT"), fill = c("#03B2C9","#DE6E4B"),cex=0.7)
#M-vlaue density MT vs WT
plot(density_wtM,main="Density of M Values",col="#03B2C9",lwd=2.5,xlab='M-val')
lines(density_mtM,main="Density of M Values",col="#DE6E4B",lwd=2.5)
legend('topright', legend=c("WT","MT"), fill = c("#03B2C9","#DE6E4B"),cex = 0.7)

dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)
#get beta for infinium I and II
beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
#raw mean
density_mean_beta_I <- density(apply(beta_I,1,mean,na.rm=T),na.rm=T)
density_mean_beta_II <- density(apply(beta_II,1,mean,na.rm=T),na.rm=T)
#raw sd
density_sd_of_beta_I <- density(apply(beta_I,1,sd,na.rm=T),na.rm=T)
density_sd_of_beta_II <- density(apply(beta_II,1,sd,na.rm=T),na.rm=T)

#apply Noob normalization
RGSet_Noob <- preprocessNoob(RGset)
#get beta for Infinium I and II
beta_Noob <- getBeta(RGSet_Noob)
beta_I_Noob <- beta[rownames(beta_Noob) %in% dfI$IlmnID,]
beta_II_Noob <- beta[rownames(beta_Noob) %in% dfII$IlmnID,]
#Noob mean
density_mean_beta_I_Noob <- density(apply(beta_I_Noob,1,mean,na.rm=T),na.rm=T)
density_mean_beta_II_Noob <- density(apply(beta_II_Noob,1,mean,na.rm=T),na.rm=T)
#Noob sd
density_sd_of_beta_I_Noob <- density(apply(beta_I_Noob,1,sd,na.rm=T),na.rm=T)
density_sd_of_beta_II_Noob <- density(apply(beta_II_Noob,1,sd,na.rm=T),na.rm=T)

#plot
par(mfrow=c(2,3),xpd=T)
plot(density_mean_beta_I,col="#003554",main="raw beta",lwd=3.0,xlab='beta')
lines(density_mean_beta_II,col="#D95CFF",lwd=3.0)

legend('bottomright',cex=0.8,inset=c(0, -0.35), legend=c("Type I","Type II"), fill = c("#003554","#D95CFF"))

plot(density_sd_of_beta_I,col="#003554",main="raw sd",lwd=3.0,xlab='beta sd')
lines(density_sd_of_beta_II,col="#D95CFF",lwd=3.0)
###

targets$Group<- as.factor(targets$Group)
palette(c("#DE6E4B","#03B2C9"))
boxplot(beta,col=targets$Group, border="#343434",names=NA,main="raw beta by Sample")
legend("bottomright", inset=c(0, -0.35), legend=levels(targets$Group),col=c(1,2),pch = 19,cex=0.9)


plot(density_mean_beta_I_Noob,col="#003554",main="preprocessNoob beta",lwd=3.0,xlab='beta')
lines(density_mean_beta_II_Noob,col="#D95CFF",lwd=3.0)

plot(density_sd_of_beta_I_Noob,col="#003554",main="preprocessNoob sd",lwd=3.0,xlab='beta sd')
lines(density_sd_of_beta_II_Noob,col="#D95CFF",lwd=3.0)

boxplot(beta_Noob,col=targets$Group,border="#343434",names=NA,main="preprocessNoob beta by Sample")


pca_results <- prcomp(t(beta_Noob),scale=T)
#Scree plot
library(factoextra)
fviz_eig(pca_results, addlabels = T,xlab='PC number',ylab='% of var', barfill = "#369AD3", barcolor = "#369AD3")


par(mfrow=c(1,1))
palette(c("#DE6E4B","#03B2C9"))
# Set shapes for sexes
sex_shapes <- c("M" = -0x2642L, "F" = -0x2640L)
# Create the plot
plot(pca_results$x[,1],pca_results$x[,2], col = targets$Group, pch = sex_shapes[targets$Sex],cex=1.5,xlab = "PC1(33%)", ylab = "PC2(22.3%)",main='PCA (Sex/Group)',xlim=c(-750,750),ylim=c(-750,750))
text(pca_results$x[,1],pca_results$x[,2],labels=rownames(pca_results$x),cex=0.4,pos=2,srt=-30)
legend("topright", legend=levels(targets$Group),col=c(1,2),pch = 19)


#speed up with parallelization
library(future.apply)
plan(multisession)
My_t_test <- function(x) {
  t_test <- t.test(x ~ targets$Group)
  return(t_test$p.value)} 
p_values <- future_apply(beta_Noob, 1, My_t_test)
final_ttest <- data.frame(beta_Noob,t_test_p_val = p_values)
final_ttest <- final_ttest[order(final_ttest$t_test_p_val),]
hist(final_ttest$t_test_p_val, main="P-value distribution (t-test)",xlab='p-val')
abline(v=0.05,col="#ef233c")

corr_pValues_BH <- p.adjust(final_ttest$t_test_p_val,"BH")
corr_pValues_bonferroni <- p.adjust(final_ttest$t_test_p_val,"bonferroni")
final_ttest_corr <- data.frame(final_ttest,corr_pValues_BH,corr_pValues_bonferroni)
colMeans(final_ttest_corr[,9:11]<0.05)*nrow(final_ttest_corr)

boxplot(final_ttest_corr[,9:11], ylim = c(-0.1, 1.1), col = c("#1D2322", "#CD523C", "#87DFD6"),names=NA,main='p-val before/after corrections')
legend("topright", legend=c("raw", "BH", "Bonferroni"),col=c("#1D2322", "#CD523C", "#87DFD6"),pch=19, cex=0.5, xpd=TRUE)


# WT group mean
WT_group_mean <- apply(final_ttest_corr[,targets$Group=="WT"], 1, mean)
# MUT group mean
MUT_group_mean <- apply(final_ttest_corr[,targets$Group=="MUT"], 1, mean)
# Compute delta between the means
delta <- MUT_group_mean - WT_group_mean
BH_sig<-final_ttest_corr[10]<0.05
toVolcPlot <- data.frame(delta,'minus_log10_p_val'= -log10(final_ttest_corr$t_test_p_val),BH_sig)

plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.4,col='#1D2322',xlab='Δμ_beta(MUT-WT)',ylab='-log10 p',xlim=c(-0.8,0.8))
abline(h=-log10(0.05))
nominal_sig <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.05)),]
BH_sig <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,3]==T,]
points(nominal_sig[,1], nominal_sig[,2],pch=16,cex=0.4,col="#324946")
points(BH_sig[,1], BH_sig[,2],pch=19,cex=0.6,col="#CD523C")
legend("topright", legend=c("not_significant", "p_val < 0.05", "BH_p_val <  0.05"),col=c("#1D2322", "#324946", "#CD523C"),pch = 19)


#

final_ttest_corr_anno <-data.frame(rownames(final_ttest_corr),final_ttest_corr)
colnames(final_ttest_corr_anno)[1] <- "IlmnID"
final_ttest_corr_anno <- merge(final_ttest_corr_anno, Illumina450Manifest_clean,by="IlmnID")
input_Manhattan <- data.frame(ID=final_ttest_corr_anno$IlmnID,CHR=final_ttest_corr_anno$CHR, MAPINFO=final_ttest_corr_anno$MAPINFO, PVAL=final_ttest_corr_anno$t_test_p_val)
head(input_Manhattan)
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "X"] <- "23"               
levels(input_Manhattan$CHR)[levels(input_Manhattan$CHR) == "Y"] <- "24"              
input_Manhattan$CHR <- as.numeric(as.character(input_Manhattan$CHR))

library(qqman)
manhattan(input_Manhattan, snp="ID",chr="CHR", bp="MAPINFO", p="PVAL",annotatePval = 0.00001,col=rainbow(24),suggestiveline=F,genomewideline=-log10(0.00001) )

library(gplots)
library(viridis)
input_heatmap=as.matrix(final_ttest[1:100,1:8])

group_color = c()
i = 1
for (name in colnames(beta)){
  if (name %in% wt){group_color[i]="#393939"}
  else{group_color[i]="#CBCBCB"}
  i = i+1
}

col_pal=colorRampPalette(c("#29A6AD","#FAB319"))(100)
heatmap.2(input_heatmap,col=viridis(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=group_color,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage",key.xlab='beta-val',key.title=NA,keysize=1,labRow=NA)
legend("topright", legend=levels(targets$Group),col=c('#CBCBCB','#393939'),pch = 19,cex=0.7)

## Single
heatmap.2(input_heatmap,col=viridis(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=group_color,density.info="none",trace="none",scale="none",symm=F,main="Single linkage",key.xlab='beta-val',key.title=NA,keysize=1,labRow=NA)
legend("topright", legend=levels(targets$Group),col=c('#CBCBCB','#393939'),pch = 19,cex=0.7)

heatmap.2(input_heatmap,col=viridis(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=group_color,density.info="none",trace="none",scale="none",symm=F,main="Average linkage",key.xlab='beta-val',key.title=NA,keysize=1,labRow=NA)
legend("topright", legend=levels(targets$Group),col=c('#CBCBCB','#393939'),pch = 19,cex=0.7)
