#!/Rscript
#author:lxq
#date:20191205
#description:差异结果注释

rm(list=ls())
args = commandArgs(T)
argNum <- length(args)
if(argNum < 2){
    cat("Usage: Rscript diff.R [type:gene/transcript] [species:GRCh37/GRCm38]\n")
    quit()
}

type = args[1]
species = args[2]

#利用ensembl数据库注释差异结果
library('biomaRt')
library("curl")

biomaRt_gene<-function(){
    #data<-volcanoPlot()
    data <- read.csv("DiffExp/gene_rawresult.csv")
    ensembl_gene_id <- gsub("\\.\\d*", "", data$gene_id)
    data2<-cbind(ensembl_gene_id,data)
    data2 <- data2[data2$style !='normal',]
    Dataset <- ifelse(species=="GRCh37","hsapiens_gene_ensembl","mmusculus_gene_ensembl")
    mart = useMart("ensembl",dataset=Dataset)#有时报错biomaRt expected a character string of length 1
    #mart <-useEnsembl(biomart = "ensembl",dataset = Dataset,mirror = "asia")
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


if(type=="gene"){biomaRt_gene()}else{biomaRt_transcript()}
