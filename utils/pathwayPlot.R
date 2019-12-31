#!/bin/env Rscript
rm(list=ls())
#library(optparse)
args=commandArgs(T)
argNum <- length(args)

if(argNum < 3){
    cat("\nArguments is not right!\n\n")
    cat("Usage: Rscript pathwayPlot.R <infile> <prefix> <outdir>\n\n")
    quit()
}

infile <- args[1]
prefix <- args[2]
outdir <- args[3]
ppi <- 300
#
library(ggplot2)
setwd(outdir)
# 设置工作路径到数据存放的文件夹下
# 读数据
pathway <- read.table(infile, header=T, sep="\t")
#
png(paste0(prefix, "_pathwayEnrich.png"), height = 6*ppi, width = 7*ppi, res=ppi)
p <- ggplot(pathway, aes(richFactor, Term))
#pp = ggplot(pathway,aes(Count,Term))
# 改变点的大小
#p + geom_point(aes(size=Count))

# 四维数据的展示
#pbubble = pp + geom_point(aes(size=Count,color=-log10(PValue)))
p <- p + geom_point(aes(size=Count, color=PValue))
#
p <- p + scale_colour_gradient(low="blue", high="red")

# 绘制pathway富集散点图
p <- p + scale_colour_gradient(low="blue", high="red") + labs(color=expression(PValue),
                                                                   size="Gene number",
                                                                   x="Rich factor",
                                                                   y="Pathway name",
                                                                   title="Pathway enrichment")
#pr = pbubble + scale_colour_gradient(low="green",high="red") + labs(color=expression(-log[10](PValue)),size="Gene number",x="Rich factor",y="Pathway name",title="pathway enrichment P<0.05")
# 改变图片的样式（主题）
p <- p + theme_bw()
p
## 保存图片
#png(paste0(prefix, "_pathwayEnrich.png",height = 6, width = 7)   # 保存为pdf格式，支持 pdf，png，svg多重格式
dev.off()
#ggsave("out.png")  # 保存为png格式
#ggsave("out2.png",width=7,height=6)   # 设定图片大小

#pathway.txt
#Term	Count	richFactor	PValue
#hsa04210:Apoptosis	6	0.028571429	0.039430005
#hsa05210:Colorectal cancer	7	0.033333333	0.010926012
#hsa05213:Endometrial cancer	8	0.038095238	9.05E-04
#hsa05230:Central carbon metabolism in cancer	8	0.038095238	0.003089379
#hsa04114:Oocyte meiosis	9	0.042857143	0.019337874
#hsa04611:Platelet activation	9	0.042857143	0.043893444
#hsa04540:Gap junction	9	0.042857143	0.005077648
#hsa05215:Prostate cancer	10	0.047619048	0.001325294
#hsa04725:Cholinergic synapse	10	0.047619048	0.006437549
#hsa04310:Wnt signaling pathway	11	0.052380952	0.009075737
#hsa04514:Cell adhesion molecules (CAMs)	11	0.052380952	0.010985732
#hsa05152:Tuberculosis	13	0.061904762	0.007473111
#hsa04390:Hippo signaling pathway	14	0.066666667	6.44E-04
#hsa04015:Rap1 signaling pathway	14	0.066666667	0.0113959
#hsa05205:Proteoglycans in cancer	15	0.071428571	0.002989335
#hsa04510:Focal adhesion	15	0.071428571	0.003894125
#hsa04810:Regulation of actin cytoskeleton	16	0.076190476	0.001751016
#hsa04144:Endocytosis	17	0.080952381	0.002622779
#hsa04151:PI3K-Akt signaling pathway	23	0.10952381	7.37E-04
#hsa05200:Pathways in cancer	25	0.119047619	7.83E-04
#

#备注，绘图数据的说明：
#（1）绘图数据来自我们公司KEGG富集分析的结果，相应文件是结题报告中存在的，略作调整即可；
#（2）绘图数据每一列的意思：
#1)Pathway      : 通路的名称        
#2)Count            ：差异表达基因中，属于这个通路的基因的数量
#3)All_Unigene        ：所有基因中属于这个通路的基因的数量  
#4)PValue            ：富集分析p值
#5)Qvalue                ：富集分析的Q值
#6)richFactor        ：在我们分析报告中，没有提供这一列，但很容易计算。是 第二列 除以 第三列得到；
#7)Pathway ID        ：通路ID  
#8)Genes                ：通路中基因的ID
#9)KOs                  ：通路中基因的KO号

#（1）在pathway名称中如有重名，这会导致错误。在表中每个pathway只能出现一次；
#（2）文本中，出现了引号会导致错误。例如，  Alzheimer's disease， Huntington's disease这样的名称。
#这两个pathway 名称中的引号需要删除。引号的出现，会导致R无法识别两个引号间的其他符号（退格符，换行符等），导致文件读取错误。
#如果一定要保留引号，则引号的内容再用引号囊括起来，例如：
#"Alzheimer's disease"，"Huntington's disease"，从而避免单个出现的引号对其他字符的影响 。
