#!/home/siyuan.wang/anaconda3/envs/rnaseq/bin/Rscript
# rm(list=ls())

version="1.0.0"

doEnrichment<-function(genes,species=NA,idType=NA){
    # print(head(genes))
    if(length(genes)==0){
        stop('No gene provided. Check filter setting.')
    }
    if(is.na(idType)){
        idType='ENSEMBL,SYMBOL'
        cat('Id type is not specified. Set to ENSEMBL,SYMBOL by default\n')
    }
    idType=unlist(strsplit(idType,split=','))
    isENS=(startsWith(genes,'ENS') & nchar(genes)>14)
    ensID=gsub('\\.\\d*$','',genes[isENS])
    nonEnsID=ensID[!isENS]
    if(is.na(species)){
        cat('Species is not specified. ')
        if('ENSEMBL' %in% idType & length(ensID)!=0){
            # isENS=(startsWith(genes,'ENS') & nchar(genes)>14)
            # ensID=gsub('\\.\\d*$','',genes[isENS])
            if(sum(startsWith(ensID,'ENSG'))>sum(startsWith(ensID,'ENSMUSG'))){
                species='hsa'
                cat('Auto detected as hsa by ENSEMBL ID.\n')
            }else if (sum(startsWith(ensID,'ENSMUSG'))>sum(startsWith(ensID,'ENSG'))) {
                species='mmu'
                cat('Auto detected as mmu by ENSEMBL ID.\n')
            }else{
                species='hsa'
                cat('Set to hsa by default.\n')
            }
        }else{
            species='hsa'
            cat('Set to hsa by default.\n')
        }
    }

    if(species=='hsa'){
        library(org.Hs.eg.db)
        db=org.Hs.eg.db
    }else if (species=='mmu') {
        library(org.Mm.eg.db)
        db=org.Mm.eg.db
    }else {
        stop('species should be one of (hsa|mmu)')
    }

    # if(is.na(idType)){
    #     isENS=(startsWith(genes,'ENS') & nchar(genes)>14)
    #     ensID=gsub('\\.\\d*$','',genes[isENS])
    #     geneSymbol=genes[!isENS]
    #     ensID=bitr(ensID, fromType = "ENSEMBL",
    #                 toType = c("ENTREZID"),
    #                 OrgDb = db)[['ENTREZID']]
    #     geneSymbol=tryCatch(bitr(geneSymbol,
    #                             fromType = "SYMBOL",
    #                             toType = c("ENTREZID"),
    #                             OrgDb = db)[['ENTREZID']],
    #                         error=function(e) {return(NULL)}
    #                     )

    #     geneID=c(ensID,geneSymbol)
    # }else{
    #     geneID=c()
    #     for(it in strsplit(idType,split=',')){
    #         if(it=='ENTREZID'){
    #             geneID=c(geneID,tryCatch(bitr(genes, fromType = idType,
    #                         toType = c("ENSEMBL"),
    #                         OrgDb = db)[['ENTREZID']],
    #                     error=function(e) {return(NULL)}
    #                 ))
    #         }else{
    #             geneID=c(geneID,tryCatch(bitr(genes, fromType = idType,
    #                         toType = c("ENTREZID"),
    #                         OrgDb = db)[['ENTREZID']],
    #                     error=function(e) {return(NULL)}
    #                 ))
    #         }
    #     }
    # }


    geneID=c()
    for(it in idType){
        # print(idType)
        if(it=='ENTREZID'){
            geneID=c(geneID,tryCatch(bitr(nonEnsID, fromType = it,
                                            toType = c("ENSEMBL"),
                                            OrgDb = db)[['ENTREZID']],
                                        error=function(e) {return(NULL)}
                                    )
                    )
        }else if (it=='ENSEMBL') {
            # print(head(genes))
            # print(head(ensID))
            geneID=c(geneID,tryCatch(bitr(ensID, fromType = it,
                                            toType = c("ENTREZID"),
                                            OrgDb = db)[['ENTREZID']],
                                        error=function(e) {print(e);return(NULL)}
                                    )
                    )
        }else{
            geneID=c(geneID,tryCatch(bitr(nonEnsID, fromType = it,
                                            toType = c("ENTREZID"),
                                            OrgDb = db)[['ENTREZID']],
                                        error=function(e) {print(e);return(NULL)}
                                    )
                    )
        }
    }
    # print(geneID)

    if(length(geneID)==0){
        stop('No ENTREZID. Check id type')
    }

    getGoResult(geneID,prefix,db=db)
    getKeggResult(geneID,prefix,species = species)
}

getGoResult<-function(genes,prefix,db=org.Hs.eg.db){
    cat('Start GO analysis\n')
    goResult<-enrichGO(genes, OrgDb = db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05)
    if(is.null(goResult)){
        cat('No GO result\n')
        return()
    }
    goResultDF<-as.data.frame(goResult)
    write.table(goResultDF,sep='\t',file=paste0(prefix,"go_result.tsv"),quote=F,row.names=F)
    goOnto<-levels(goResultDF$ONTOLOGY)
    goResultTop<-NULL
    for(i in 1:length(goOnto)){
        tempGoResult <- goResultDF[goResultDF$ONTOLOGY==goOnto[i],]
        goResultTop <- rbind(goResultTop,tempGoResult[1:min(dim(tempGoResult)[1],10),])
    }
    goResultTop$p.adjust <- -log10(goResultTop$p.adjust)
    goResultTop<-goResultTop[dim(goResultTop)[1]:1,]
    goResultTop$Description<-factor(goResultTop$Description,levels=factor(goResultTop$Description,levels=goResultTop$Description))

    png(paste0(prefix,"go_result_dotPlot.png"),width=800,height=1200)
    print(dotplot(goResult,showCategory= min(c(dim(goResult)[1],30))))
    dev.off()

    png(paste0(prefix,"go_result_parts_barplot.png"),width=1200,height=1200)
    p<-ggplot(goResultTop,aes(x=Description,y=p.adjust,fill = ONTOLOGY)) + # ggtitle(paste0(clName,'')) +
            geom_bar(stat = "identity", position = "dodge", width = 0.8) + coord_flip()+
            theme_light(base_size=12) +
            theme(text=element_text(size=14,  face="bold"),
                strip.background = element_blank(),
                panel.border     = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.position = "right")+
            xlab("Description") +
            ylab("-log10(p.adjust)")
    print(p)
    dev.off()

}


getKeggResult<-function(genes,prefix,species = 'hsa'){
    cat('Start KEGG analysis\n')
    keggResult<-enrichKEGG(gene=genes,organism = species,pvalueCutoff=0.05)
    if(is.null(keggResult)){
        cat('No KEGG result\n')
        return()
    }
    keggResultDF<-as.data.frame(keggResult)
    write.table(keggResultDF,quote=F,sep='\t',file=paste0(prefix,'kegg_result.tsv'),row.names=F)


    png(paste0(prefix,'KEGG.png'),width=800,height=1200)
    print(dotplot(keggResult,showCategory=min(c(dim(keggResult)[1],30))))
    dev.off()


    library(pathview)
    # genePd=as.integer(genes)
    cwd=getwd()
    setwd(outdir)
    pathways=rownames(keggResultDF)

    pwIndex=1
    successedNum=0
    while(TRUE){
        pw=pathways[pwIndex]
        isError=FALSE
        tryCatch(pathview(gene.data = genes,
                            pathway.id = pw,
                            species = species),
                    error=function(e){cat(paste('Fail to get pathway',pw,'\n'));print(e);isError=TRUE})
        if(!isError){
            successedNum=successedNum+1
        }
        pwIndex=pwIndex+1
        if(successedNum>=10 | pwIndex>length(pathways)){
            cat(paste('Successfully get',successedNum,'pathway graph\n'))
            break
        }
    }
    setwd(cwd)
    # print('kegg done')
}

# source('newArgu.R')

# short option(string|NA), long option(string|NA), description(string), default value(string|NA|bool), required(bool)
spec = matrix(c(
    'outdir', 'o', 1, 'character', 'output directory',
    'genelist', 'i', 1, 'character', 'input gene list',
    'prefix', 'n', 1, 'character', 'output file prefix',
    'minabslogfc', 'm', 1, 'double', 'min absolute value log2fc (1)',
    'maxpadj', 'M', 1, 'double', 'max p.adj (0.05)',
    'offlogfc', 'L', 0, 'logical', 'disable filt gene by log2fc (enable)',
    'offpadj', 'P', 0, 'logical', 'disable filt gene by padj (enable)',
    'species', 's', 1, 'character', 'species hsa|mmu (auto detect|hsa)',
    'noheader', 'H', 0, 'logical', 'no header for genelist (disable)',
    'idtype', 't', 1, 'character', 'gene id type, splited by \',\' (ENSEMBL,SYMBOL)',
    'help', 'h', 0, 'logical', 'print this message',
    'version', 'V', 0, 'logical', 'print version message'
), byrow=TRUE, ncol=5)
# optionList=data.frame(
#         outdir=c('-o','--outdir','output directory',NA,TRUE),
#         geneList=c('-i','--infile','input gene list',NA,TRUE),
#         prefix=c('-n','--prefix','output file prefix','',FALSE),
#         minAbsLogFc=c('-f','--abslogfc','min absolute value log2fc','1',FALSE),
#         maxPadj=c('-p','--padj','max p.adj','0.05',FALSE),
#         logFcFilter=c('-l','--disable-logfc','disable filt gene by log2fc',TRUE,FALSE),
#         padjFilter=c('-f','--disable-padj','disable filt gene by padj',TRUE,FALSE),
#         species=c('-s','--species','species hsa|mmu','hsa',FALSE),
#         noHeader=c('-H','--no-header','no header',FALSE,FALSE),
#         idType=c('-t','--id-type','gene id type',NA,FALSE),
#         help=c('-h','--help','print this message',FALSE,FALSE),
#         version=c('-V','--version','print version message',FALSE,FALSE),
#         row.names=c('shortOpt','longOpt','description','defaultValue','required')
#     )
defaultValue=list(
    outdir=NA,
    genelist=NA,
    prefix='',
    minabslogfc=1,
    maxpadj=0.05,
    offlogfc=FALSE,
    offpadj=FALSE,
    species=NA,
    noheader=FALSE,
    idtype=NA,
    help=FALSE,
    version=FALSE
)
isrequired=c('outdir','genelist')

# showUsage(options)
library(getopt)

argv=commandArgs(T)
# argv='-o test2/ -i genes.txt -s mmu -n test -f -F'
# argv=strsplit(args,split=' ')
# argv=c('-o','./test','-i','./gene-deseq2-results.csv','-n','test','--species','mmu')
# print(length(argv))
opts=getopt(spec,opt=argv)

# print(argv)
# print(opts)
# print(isrequired %in% names(opts))
if (!is.null(opts$version)) {
    cat(paste0('enrich.R\nVersion:  ',version,'\n'))
    q()
}else if(length(argv)==0 | !is.null(opts$help)){
    cat(getopt(spec,usage=TRUE,command="enrich.R"))
    q()
}else if (sum(isrequired %in% names(opts))!=length(isrequired)) {
    cat(paste(isrequired[!(isrequired %in% names(opts))],'is required.\n'))
    # cat(getopt(spec,usage=TRUE,command="enrich.R"))
    q()
}

for(opt in names(opts)){
    defaultValue[[opt]]=opts[[opt]]
}
opts=defaultValue
# opts=getOptions(args,optionList)
# print(opts)
outdir=opts$outdir
geneList=opts$genelist
prefix=opts$prefix
if(prefix!=''){
    prefix=paste0(prefix,'_')
}
species=opts$species
idType=opts$idtype
minAbsLogFc=opts$minabslogfc
maxPadj=opts$maxpadj
logFcFilter=!opts$offlogfc
padjFilter=!opts$offpadj
noHeader=opts$noheader

# print(opts)
if(!endsWith(outdir,'/')){
    outdir=paste0(outdir,'/')
}

if(!dir.exists(outdir)){
    dir.create(outdir)
}
if(noHeader){
    geneTable=read.csv(geneList,header=F,row.names=1)
    # head(geneTable)
    # rownames(geneTable)=geneTable$V1
}else{
    geneTable=read.csv(geneList,header=T,row.names=1)
}
geneTableFilter=geneTable

if(logFcFilter & 'log2FoldChange' %in% colnames(geneTableFilter)){
    cat('Filter by logFc\n')
    geneTableFilter=geneTableFilter[!is.na(geneTableFilter$log2FoldChange) & (geneTableFilter$log2FoldChange<=-minAbsLogFc | geneTableFilter$log2FoldChange>=minAbsLogFc),]
}
if(padjFilter & 'padj' %in% colnames(geneTableFilter)){
    cat('Filter by padj\n')
    geneTableFilter=geneTableFilter[!is.na(geneTableFilter$padj) & geneTableFilter$padj<=maxPadj,]
}
genes=rownames(geneTableFilter)
prefix=paste0(outdir,prefix)


library(clusterProfiler)
library(enrichplot)
library(ggplot2)


doEnrichment(genes,species,idType)

