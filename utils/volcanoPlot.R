#!/bin/Rscript
library(ggplot2)
data<- read.table("allgene.dif.txt",header = T, row.names = 1, sep = "\t")
#"logFC"	"AveExpr"	"t"	"P.Value"	"adj.P.Val"	"B"
#"IRF2BPL"	-4.49475248073427	6.78342680319745	-37.4279677343214	8.58519539610212e-11	1.73274998679529e-06	14.3372634481319
#"ZBTB18"	-3.38423495825277	3.63053075195588	-32.5748646789161	2.79892471591668e-10	2.47328861298219e-06	13.5638377769904
#"MGC70870"	-4.13748581745129	2.87164834313146	-30.594992073769	4.76872837061703e-10	2.47328861298219e-06	13.1868441561119
#"PDCD10"	-3.31681156399585	5.32424438398808	-30.0085582815514	5.62057483559919e-10	2.47328861298219e-06	13.0671411066062
#"EPC2"	-3.85154305193362	4.75730039796338	-29.7050797082289	6.12715803642222e-10	2.47328861298219e-06	13.0036549449308
#"RAP2C"	-4.36007946149376	5.14179695122444	-28.53468437609	8.61803742676494e-10	2.89896415640661e-06	12.7484893741725
#"TUBA1C"	-2.5559184826665	9.04009950681825	-27.1975838340864	1.29460318961692e-09	3.73271088229119e-06	12.4355130176857
#"DHX15"	-2.61265061063463	6.52744888765385	-25.026297859527	2.61874248139515e-09	6.0618230410919e-06	11.8725549095367
#"CHIC2"	-3.80374522660045	5.87216541360614	-24.8083483680335	2.81983536348457e-09	6.0618230410919e-06	11.811941173012
#"STRN3"	-2.6744529705801	4.15960471627323	-24.3049142620383	3.35356475149442e-09	6.0618230410919e-06	11.6688496530652
#"TMEM60"	-2.7234394243256	5.31764733126209	-24.1678924533591	3.51770878841514e-09	6.0618230410919e-06	11.62914215453
#"CALM3"	-2.57909800325415	6.012166857344	-24.0985982784826	3.6041161617749e-09	6.0618230410919e-06	11.6089349974241


#plot(logFC,-log10(FDR),pch=16,cex=0.5,xlab="logFC",ylab="-lg(FDR)",cex.axis=0.7,cex.lab=0.9,main="Volcano plot",
#     xlim = c(-3,3),ylim = c(0,20),
#     col="#BCBABE") #col="#BCBABE" 灰色
#pch=16的意思是选择了一种点的样式，不同样式的点对应着不同的编号。
#用xlim=c(-4,4)参数限定X轴只取[-4,4]； 用ylim=c(0,4)参数限定Y轴只取[0,4]
down<- intersect(which(data$P.Value<0.05) , which(data$logFC<=(-1)))
up<- intersect(which(data$P.Val<0.05) , which(data$logFC>=1))

DEGs<- rep("normal",times= nrow(data))

DEGs[down]<-"down"
DEGs[up]<-"up"
DEGs<- factor(DEGs,levels = c("up","down","normal"))

p<- qplot(x=data$logFC,y=-log10(data$P.Value),xlab = "logFC",ylab = "-lg(Pvalue)",size=I(0.7),col=DEGs)
#p<- qplot(x=data$logFC,y=-log10(data$adj.P.Val),xlab = "logFC",ylab = "-lgFDR",size=I(0.7),col=DEGs,xlim = c(-3,3),ylim = c(0,12))

p<- p+ scale_color_manual(values = c("up"="red","down"="green","normal"="black"))
xline=c(-1,1)
p<- p+geom_vline(xintercept = xline, lty=2,size=I(0.2),col="grey11")
yline= -log(0.05,10)
#yline= 1.5
p<- p+geom_hline(yintercept = yline, lty=2,size=I(0.2),col="grey11")

p<- p+theme_bw()+theme(panel.background = element_rect(colour = "black",size = 1,fill = "white"),panel.grid = element_blank())

pdf("vlocana.plot.pdf",height = 5,width = 6 )
print(p)
dev.off()

