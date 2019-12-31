#R包安装 R 3.6以下版本
#source("http://bioconductor.org/biocLite.R")
#biocLite(“affy","annotate")

#bioconductor安装Ｒ包。３.６版本以上
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("limma")



###############
#差异分析
###############
library(limma)
#data<-expr.matrix #数据读入
data<-read.table("uniquenorm.txt",sep = "\t")
data<- as.matrix(data)#将data转换为数值Number格式
#data<- log2(data)  # log的目的是对数据标准化。


mode(data)


group<-c(rep(1,3),rep(2,3)) #分组 ,变动最大
#group2<- c(rep(1,3),rep(2,3),rep(3,3))
#group<-c(gl(2,3))
#group<- c(2,2,1,1,1,2,2,1,1,2)
#group<-c(1,1,1,2,2,2)
#group<-c(rep(c(1,2),5))
design <- model.matrix(~ -1+factor(group)) #设计矩阵
#design <- model.matrix(~ -1+factor(group2)) #设计矩阵
colnames(design)<-c("case","con") #命名.
#colnames(design)<-c("Tumor","Normal") #命名
#colnames(design)<-c("case","con","treat") 
design #查看确认分组顺序是否正确
#配对信息

#case vs control 的差异分析
contrast.matrix <- makeContrasts(case-con,levels=design) # 差异分析的关键，case，con是前面命名的组名。
#contrast.matrix <- makeContrasts(Tumor-Normal,levels=design) # 差异分析的关键，case，con是前面命名的组名。

#建立线性模型，并进行分组差异分析和p值计算、校正
fit <- lmFit(data, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
#results<-decideTests(fit2, method="global", adjust.method="none", p.value=0.05, lfc=0.585)
#summary(results)
dif <- topTable(fit2, coef = 1,  n = nrow(fit2), lfc = 0) 
#dif <- topTable(fit2,coef=1,number=nrow(fit2) ,adjust.method="BH",sort.by="B",resort.by="M")
write.table(dif, "allgene.dif.txt", row.names = TRUE, sep = "\t")#所有基因的差异表达情况输出。

#根据自己需求调整差异基因筛选阈值，也可进入excel筛选
dif1 <- dif[dif[, "P.Value"] < 0.05, ] 
#dif1 <- dif[dif[, "P.Value"] < 0.01, ] 
dif1 <- dif1[abs(dif1$logFC)> 1,]
write.table(dif1, "dif1.txt", row.names = TRUE, sep = "\t")


#plot
x<- dif$logFC
y<- -log10(dif$P.Value)

pdf("Fig3_plot.pdf", height= 8, width=8)
plot(x,y,xlab = "logFC",ylab = "-lgP", main = "GSE108810")
dev.off()


