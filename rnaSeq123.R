#!/bin/Rscript
#WrittenBy:liang.wan, 2019-12-09
#https://www.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html


rm(list=ls())
args <- commandArgs(T)
argNum <- length(args)
#f(argNum <= 6){
#   cat("\nError Options!\n\n")
#   cat("Usage: Rscript rnaSeq123.R <Human|Mouse> <countsFileFolder> <countsFilePattern> <prefix>")
#   quit()
#}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#load library
library(limma)
library(Glimma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(Mus.musculus)
#record log
sink(paste0(prefix, "rnaSeq123Analysis.log"), append = TRUE)
#
ppi <- 300
#tmp
countsFileFolder <- '/home/liang.wan/workdir/Dec/rnaSeq-test/rawdata'
countsFilePattern <- '^GSM.*'
if(length(files = list.files(path=countsFileFolder, pattern = paste0(countsFilePattern,'.txt$'), full.names=TRUE)) > 0){
    files = list.files(path=countsFileFolder, pattern = paste0(countsFilePattern,'.txt$'), full.names=TRUE)
}else{
    cat("countsFileFolder or countsFilePattern is ERROR! Job failed!")
    quie()
}
#read DGE counts matrix
x <- readDGE(files, columns=c(1,3))
dim(x)
samplenames <- substring(colnames(x), nchar(colnames(x))+13, nchar(colnames(x)))
colnames(x) <- samplenames
#组织样本信息，包括分组，不同lane测序信息
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
#
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
#
x$samples
#组织基因注释
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
head(genes)
genes <- genes[!duplicated(genes$ENTREZID),]
x$genes <- genes
x
#
#数据预处理
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE, prior.count=2)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)
#删除低表达的基因
table(rowSums(x$counts==0)==9)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
#data Fillter
png(paste0(prefix, "_DataFilter.png"), width = 8*ppi, height = 4*ppi, res = ppi)
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()
#归一化基因表达数据，去除系统误差
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
png(paste0(prefix, "_Normalised.png"), width = 4*ppi, height = 4*ppi, res=ppi)
lcpm <- cpm(x, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main=paste0(prefix, " normalised data"),ylab="Log-cpm")
dev.off()
#对样本进行无监督聚类MDS
png(paste0(prefix, "_MDS.png"), width = 4*ppi, height = 4*ppi, res=ppi)
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main=paste0(prefix, " sample groups"))
dev.off()
#差异表达分析
#创建设计矩阵和对比
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design
tmpMatrix <- t(combn(unique(as.character(group)), 2))
tmpVector <- c()
for(row in 1:dim(tmpMatrix)[1]){
    tmpVector <- append(tmpVector, paste0(tmpMatrix[row,1], "-", tmpMatrix[row,2]))
}
contr.matrix <- makeContrasts(contrasts = tmpVector, levels = colnames(design))

#从表达计数数据中删除异方差
png(paste0(prefix, "_Model.png"), width = 4*ppi, height = 4*ppi, res=ppi)
v <- voom(x, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()
#检查DE基因数量
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
#
png(paste0(prefix, "_DEGVenn.png"), width = 8*ppi, height = 8*ppi, res=ppi)
vennDiagram(dt, circle.col=rainbow(7))
dev.off()
#
write.fit(tfit, dt, file=paste0(prefix, "_AllDEGResults.txt"))
#逐个检查单个DE基因，写入文件
tmp <- colnames(summary(dt))
for(i in 1:dim(summary(dt))[2]){
    write.csv(topTreat(tfit, coef=i, n=Inf), file = paste0(tmp[i], "_DEGResult.csv"), row.names = FALSE)
}
#差异表达结果的实用图形表示
for(i in 1:dim(summary(dt))[2]){
    png(paste0(tmp[i], "_MD.png"), width = 8*ppi, height = 8*ppi, res=ppi)
    plotMD(tfit, column=i, status=dt[,i], main=colnames(tfit)[i])
    dev.off()
}
#DETop基因进行聚类
for(i in 1:dim(summary(dt))[2]){
    png(paste0(tmp[i], "_ClusterHeatMap.png"), width = 6*ppi, height = 15*ppi, res=ppi)
    topGenes <- topTreat(tfit, coef=i, n=Inf)$ENTREZID[1:100]
    g <- which(v$genes$ENTREZID %in% topGenes)
    mycol <- colorpanel(1000,"blue","white","red")
    heatmap.2(lcpm[g,], scale="row",
              labRow=v$genes$SYMBOL[g], labCol=group, 
              col=mycol, trace="none", density.info="none", 
              margin=c(8,6), lhei=c(2,10), dendrogram="column")
    dev.off()
}
#






























