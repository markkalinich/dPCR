#The purpose of this script is to generate ROC curves and establish which genes are co-correlated from Hong et al.'s 2018 PNAS paper. Questions, comments, concerns 
#can be directed to Mark Kalinich, markkalinich [at] gmail [dot] com. This script was initially written on 20160803 and uploaded to github on 20180202. 

#NOTE: This code requires a package available from Ben Wittner; heat map will not work without it. 

#Step 1: Load Required libraries -----------------------------------------
library("bswHeat")
library(ROCR)
library(corrplot)
library(stargazer)
date <-Sys.Date()
setwd("/Users/markkalinich/")
#Step 2: Import melanoma datasets + make heatmap -----------------
master_df <- read.csv("Dropbox (MIT)/PhD/R/Melanoma/20160920ROC.csv")
master_df$ID <- 1
master_df$ID[grepl("HD", master_df$Patient_ID)] <-0 #iff Brx then label as cancer

set.seed(5739)
rowDF <- subset(master_df, select = c(Patient_ID,FAT1:SignalTransduction))
rowDF$Patient_ID <-paste(master_df$Patient_ID,master_df$DrawDate,sep = "_")
names(rowDF)[names(rowDF) == 'TotalTranscriptsPerWellPerMlBlood'] <- 'Total' 
names(rowDF)[names(rowDF) == 'LineageSpecific'] <- 'Lineage' 
names(rowDF)[names(rowDF) == 'CarcinoEmbryonicAntigen'] <- 'CEA' 
names(rowDF)[names(rowDF) == 'SignalTransduction'] <- 'SigTrans' 
size = dim(rowDF)
namesCol <- colnames(rowDF)[-1] #gene names 
colCat <- data.frame(name = namesCol,
                     nameToPlot = namesCol,
                     category = c(rep('Genes',size[2]-4),rep('Groupings',3)), #c(rep('HD',26),rep('HBV',11),rep('HCV',4),rep('EtOH',6),rep('Other',4),rep('Untreated',16),rep('Treated',32),rep('NED',14),rep('ICC',6),rep('PDAC',6),rep('Breast',5),rep('Lung',8),rep('Prostate',6),rep('Melanoma',8),rep('MET',4)),
                     stringsAsFactors = FALSE)
rownames <- as.character(rowDF$Patient_ID) #gene names
rowCat <- data.frame(name = rownames,
                     category = c(rep("Patient",size[1] )),
                     stringsAsFactors = FALSE)

mydata <- rowDF[,-1]
mtx <- as.matrix(mydata)
rownames(mtx) <- rowDF$Patient_ID
mtxLog <- log10(1+mtx)     

maxval <-0
for (k in 1:ncol(mtxLog)){
  maxval[k]<-max(mtxLog[,k],na.rm = TRUE)
}
mtxLog_norm <-t(t(mtxLog)/maxval)
rownames(mtxLog_norm) <- rowDF$Patient_ID

pdfname<-paste(date,"MelanomaHeatmap")
bswHeat.pdf(pdfname, "numeric", mtxLog, rowCat = rowCat, colCat = colCat,showRampXlab = 'log10(transcript/mL+1)',autoZlimMethod = 'range')
pdfname<-paste(date,"MelanomaHeatmapNORM")
bswHeat.pdf(pdfname, "numeric", mtxLog_norm, rowCat = rowCat, colCat = colCat,showRampXlab = 'Normalized log10(transcript/mL+1)',autoZlimMethod = 'range')

#Step 3: ROC Curve analysis -----------------

op <- par() # save default par settings
pdfname = paste(date,"Melanoma ROC GENE SETS",sep=" ")
pdf(pdfname) # define name of pdf. Output after this will be saved into the pdf until dev.off()
#par(mfrow=c(4,6)) # i want this page to show 3 graphs x 4 graphs 

ROC <- subset(master_df,select=c(Patient_ID,ID,FAT1:SignalTransduction))
ROC[is.na(ROC)] <- 0
p_val_vector =  0
for (k in 3:ncol(ROC)){
  geneDF <- data.frame(cbind(Gene=(log10(ROC[,k]+1)),ID=ROC$ID))
  model<-glm(ID ~ Gene, data = geneDF, family=binomial(link="logit"))
  pred <- prediction(geneDF$Gene,geneDF$ID)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr",cex.lab=7) 
  plot(perf, col=rainbow(10),main = colnames(ROC[k]), xlab="FPR",ylab="TPR",mgp=c(2.3,1,0),cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  auc <- round(unlist(slot(performance(pred,"auc"), "y.values")),digits = 2)  # get AUC for my ROC curve 
  auct <- paste(c("AUC="),auc,sep="")
  pval = round(coef(summary(model))[2,4],digits = 4) #extract p-value from logistic regression for my gene of interest
  pvalt <-paste(c("p="),pval,sep="")
  p_val_vector[k]<-pval 
  # Add text using the following Corner_text function:
  Corner_text <- function(text, location="bottomright"){
    legend(location,legend=text, bty ="n", pch=NA,cex=1.8)}
  Corner_text(text=c(auct,pvalt))} #for loop for all the genes
dev.off() # stop saving stuff to pdf. This will also reset par() to default settings. 

#Step 4: Correlogram -----------
dPCRcor = as.matrix(na.omit(log10(subset(master_df,select = (FAT1:SignalTransduction))+1)))
mcor = cor(dPCRcor, method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
} #determine significance of correlation
p.mat <- cor.mtest(dPCRcor) # matrix of the p-value of the correlation
hmcols<-colorRampPalette(c("blue","white","red"))(256)

pdfname = paste(date,"Melanoma Correlogram GENE SETS",sep=" ")
pdf(pdfname) # define 
corrplot(mcor, method="color", p.mat = p.mat, sig.level = 0.05, col=hmcols,tl.col="black") #plot numbers
dev.off() # stop saving stuff to pdf. This will also reset par() to default settings.

#Step 5: Make a multi-predictor model out of the genes ------
p_val_vector[1:2] <-1.1
sig_genes <- ROC[,which(p_val_vector<.05)]
sig_genes$ID <-ROC$ID 
multi_model<-glm(ID ~ ., data = sig_genes, family=binomial(link="logit")) 
stargazer(multi_model, type="html",title="Multi-gene Logistic Regression Results", align=TRUE)

pdfname = paste(date,"Melanoma Multi-Model GENE SETS",sep=" ")
pdf(pdfname) # define name of pdf. Output after this will be saved into the pdf until dev.off()
#See what model does without cross-validation NOTE: problem extracting pvalue 
NO_CV_Prob <-predict(multi_model,sig_genes[,1:(ncol(sig_genes)-1)],type="response")
pred <- prediction(NO_CV_Prob, ROC$ID)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="No Cross-Val",xlab="FPR",ylab="TPR",mgp=c(2.3,1,0),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
auc <- round(unlist(slot(performance(pred,"auc"), "y.values")),digits = 2)  # get AUC for my ROC curve 
auct <- paste(c("AUC="),auc,sep="")
#pval = round(coef(summary(multi_model))[2,4],digits = 4) #extract p-value from logistic regression for my gene of interest
#p_val_vector[k]<-pval
#pvalt <-paste(c("p="),pval,sep="")
# Add text using the following Corner_text function:
Corner_text <- function(text, location="bottomright"){
  legend(location,legend=text, bty ="n", pch=NA,cex=1.1)}
Corner_text(text=c(auct))
#Perform LOOCV on dataset
LOOCV <- data.frame(matrix(0, ncol = 1, nrow = nrow(sig_genes)))
colnames(LOOCV)[1]<-'Cxr Probability'
for (i in 1:nrow(sig_genes)){
  TrainDF = sig_genes[-i, ]
  TestDF <- sig_genes[i, ]
  GLMTrain <-glm(ID ~ ., data = TrainDF, family=binomial(link="logit")) 
  LOOCV[i,1] <- predict(GLMTrain,TestDF,type="response")}
ROC$Odds <-sapply(LOOCV, function(p){p/(1-p)})

pred <- prediction(LOOCV$`Cxr Probability`, ROC$ID)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10), main="LOOCV",xlab="FPR",ylab="TPR",mgp=c(2.3,1,0),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
auc <- round(unlist(slot(performance(pred,"auc"), "y.values")),digits = 2)  # get AUC for my ROC curve 
auct <- paste(c("AUC="),auc,sep="")
#pval = round(coef(summary(model))[2,4],digits = 4) #extract p-value from logistic regression for my gene of interest
#p_val_vector[k]<-pval
#pvalt <-paste(c("p="),pval,sep="")
# Add text using the following Corner_text function:
Corner_text <- function(text, location="bottomright"){
  legend(location,legend=text, bty ="n", pch=NA,cex=1.1)}
Corner_text(text=c(auct))
dev.off() # stop saving ROC curves into .pdf NOTE:
