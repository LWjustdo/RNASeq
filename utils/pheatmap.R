#1. pheatmap包的安装
#install.packages("pheatmap")
#如果该包已经安装成功，那么就不用这两步了。

#2. 表达谱矩阵文件的导入。
library(pheatmap)  # pheatmap包调入

#提取差异表达基因的表达矩阵
dif<-read.table("allgene.dif.txt",sep = "\t")
expr.matrix<- read.table("uniquenorm.txt",sep = "\t")
dif_heat<-rownames(dif[which(abs(dif$logFC)>0.585 & dif$P.Value<0.05),]) #提取差异基因名称列表
#dif_heat<- read.table("DEG.symbol.txt",row.names = 1) #这个txt文件可以直接是差异基因名称列表.
#dif_heat<- rownames(dif_heat)
#expr.matrix<- read.table("uniquenorm.txt",sep = "\t")
dif_exprs_heatmap<-expr.matrix[which(rownames(expr.matrix) %in% dif_heat),]# 获得差异基因的表达矩阵。
write.table(dif_exprs_heatmap,"dif_heatmap.txt",sep = "\t")
#
#GSM2913531_3B_CTL1.CEL.gz	GSM2913532_3B_CTL2.CEL.gz	GSM2913533_3B_CTL3.CEL.gz	GSM2913534_3B_M1.CEL.gz	GSM2913535_3B_M2.CEL.gz	GSM2913536_3B_M3.CEL.gz
#A1BG	2.32565320059249	1.94436249466558	1.8966063340388	2.11023162204871	1.89695763599491	2.81366763883655
#A1BG-AS1	4.51257193628592	4.13918065905322	4.51257193628592	5.20972290469991	4.77360018276897	5.45323139753985
#A1CF	2.13984357587726	1.50233442472211	1.83779999198207	2.0836808636128	2.30907378456679	2.1267195426419
#A2M	3.01779334140216	3.01429146814953	2.88037548045945	2.59977935701409	2.97962832578069	2.9396863598125
#A2M-AS1	1.22472663421002	1.47017757348227	1.66322733780963	1.16718022918204	1.11187638174351	1.43213552740084
#A2ML1	1.95386828213689	2.1236281982021	2.06942576100061	2.44177267997536	2.39941759072464	2.51410969953112
#A2MP1	3.71634371957077	3.66006381581805	4.12738497687362	3.19429042926572	3.71634371957077	3.7726236233235
#A4GALT	3.88419070560459	4.1109175700953	3.77533913156783	5.19519511790899	4.56785553030293	5.08952682003976

#dif_exprs_heatmap<- read.table("dif_heatmap.txt",header = T, row.names = 1)

#预设置热图注释
annotation_row <- data.frame(Group = c(rep("Case", 3), rep("Ctrl", 3)))
row.names(annotation_row) <- colnames(dif_exprs_heatmap)
annotation_col <- list( Group = c("Ctrl"="orange","Case"="blue") )
#4. 热图构建
pdf("dif_heatmap.pdf", onefile = F, width = 6, height = 6)
pheatmap(dif_exprs_heatmap,
         scale = "row", #是否做Z值转换
         color = colorRampPalette(c("green", "black", "red"))(80),
         #cutree_rows = 2,     #将热图分为高低表达两部分
         cluster_cols=T, #是否横向聚类
         cluster_distance_cols ="correlation", #纵向聚类方式
         clustering_distance_rows = "correlation",# 横向聚类方式 
         annotation = annotation_row,
         show_rownames = F, fontsize_col = 6, #是否显示纵轴基因名称？
         annotation_colors = annotation_col, #分组对应颜色
         border_color = NA) 
dev.off()
