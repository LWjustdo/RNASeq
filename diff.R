#!/Rscript
#author:lxq
#date:20191205
#description:差异分析、火山图、差异结果注释

rm(list=ls())
args = commandArgs(T)
argNum <- length(args)
if(argNum < 2){
    cat("Usage: Rscript diff.R [type:gene/transcript] [repeatNum:1/2/3] [species:GRCh37/GRCm38]\n")
    quit()
}

type = args[1]
repeatNum = args[2]
species = args[3]


#判断目录Diff存在，否则新建
#if (file_test("-d","DiffExp")){print("DiffExp is exist")}else{dir.create("DiffExp")}


#差异分析
library(DESeq2)
deseq2File<-function(){
    tname <- ifelse(type=="gene","gene_id","transcript_id")
    file <-paste("Quantification/",type,"_count_matrix.csv",sep="")
    database <- as.matrix(read.csv(file, row.names= tname))
    condition <- factor(c(rep("control",repeatNum),rep("treat",repeatNum)))
    coldata <- data.frame(row.names = colnames(database), condition)
    dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
    keep <- rowSums(counts(dds)) >= 1 #过滤表达量均为0的数据
    dds <- dds[keep,]
    dds2 <- DESeq(dds)
    res <- results(dds2)
    resdata <- merge(as.data.frame(res), as.data.frame(counts(dds2, normalized=TRUE)),by="row.names", sort=FALSE)
    resdata <- na.omit(resdata)
    resdata$style <- ifelse(abs(resdata$log2FoldChange) >=1,ifelse(resdata$log2FoldChange >=1,'up','down'),'normal')
    names(resdata)[1]<-paste(type,"_id",sep="")
    outfile <- paste("DiffExp/",type,"_rawresult.csv", sep="")
    write.csv(resdata,file=outfile,row.names=FALSE)
    #return(resdata)
}
deseq2File()

#火山图
library("ggplot2")
library("ggrepel")
volcanoPlot<-function(){
    #data<-deseq2File()
    infile<-paste("DiffExp/",type,"_rawresult.csv", sep="")
    data <- read.csv(infile)
    resdata <- na.omit(data)  #删除掉包含NA的行
    threshold <- as.factor(ifelse(resdata$padj < 1 & abs(resdata$log2FoldChange) >= 1 ,ifelse(resdata$log2FoldChange >= 1 ,'Up','Down'),'Normal'))
    deg_img <- ggplot(resdata,aes(x=resdata$log2FoldChange,y=-log10(resdata$pvalue),colour=threshold)) + xlab("log2(Fold Change)")+ylab("-log10(qvalue)") + geom_point(size = 0.5,alpha=1) + ylim(0,5) + xlim(-15,15) + scale_color_manual(values=c("green","grey", "red"))
    ##添加阈值线
    line_valcao <- deg_img+geom_hline(yintercept = 40,linetype="dotted")+geom_vline(xintercept=c(-1,1),linetype="dotted")
    outfig <- paste("DiffExp/",type,"_volcano.png", sep="")
    png(outfig)
    print(deg_img)
    dev.off()
    #return(data)
}
volcanoPlot()


#利用ensembl数据库注释差异结果
#library('biomaRt')
#library("curl")

biomaRt_gene<-function(){
    #data<-volcanoPlot()
    data <- read.csv("DiffExp/gene_rawresult.csv")
    ensembl_gene_id <- gsub("\\.\\d*", "", data$gene_id)
    data2<-cbind(ensembl_gene_id,data)
    data2 <- data2[data2$style !='normal',]
    Dataset <- ifelse(species=="GRCh37","hsapiens_gene_ensembl","mmusculus_gene_ensembl")
    #mart = useMart("ensembl",dataset=Dataset)#有时报错biomaRt expected a character string of length 1
    mart <-useEnsembl(biomart = "ensembl",dataset = Dataset,mirror = "asia")
    my_ensembl_gene_id<-data2$ensembl_gene_id
    symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name', 'start_position','end_position','strand'),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
    data3<-merge(data2,symbols,by="ensembl_gene_id") #差异结果与注释结果合并
    data3$length <-(data3$end_position-data3$start_position+1) #计算基因长度
    data4<-data3[,c("gene_id","external_gene_name","chromosome_name","start_position","end_position","length","strand","log2FoldChange","pvalue","padj","style")]
    write.csv(data4,file="DiffExp/gene_diff_results.csv",row.names=FALSE) #提取指定列并按顺序生成
}

biomaRt_transcript<-function(){
    #data<-volcanoPlot()
    data <- read.csv("DiffExp/transcript_rawresult.csv")
    ensembl_transcript_id <- gsub("\\.\\d*", "", data$transcript_id)
    data2<-cbind(ensembl_transcript_id,data)
    data2 <- data2[data2$style !='normal',]
    Dataset <- ifelse(species=="GRCh37","hsapiens_gene_ensembl","mmusculus_gene_ensembl")
    mart <-useEnsembl(biomart = "ensembl",dataset = Dataset,mirror = "asia")
    my_ensembl_transcript_id<-data2$ensembl_transcript_id
    mms_symbols<- getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id',"external_gene_name",'transcript_length','transcript_biotype','external_gene_name','transcript_length','transcript_biotype'),filters = 'ensembl_transcript_id', values = my_ensembl_transcript_id, mart = mart)
    data3<-merge(data2,mms_symbols,by="ensembl_transcript_id") 
    data3<-data3[,c("transcript_id",'ensembl_gene_id',"external_gene_name",'transcript_biotype','transcript_length',"log2FoldChange","pvalue","padj","style")]
    write.csv(data3,file="DiffExp/transcript_diff_results.csv",row.names=FALSE)
}
#if(type=="gene"){biomaRt_gene()}else{biomaRt_transcript()}




