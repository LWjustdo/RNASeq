#!/home/siyuan.wang/anaconda3/envs/rnaseq/bin/Rscript

version="1.0.0"

doCluste=function(mtr,samples,outdir='./',cmethods='average',prefix=''){
    # if(doAll){
    #     countMtrFiltExp=mtr[c(conGrp,expGrp)]
    #     for(method in cmethods){
    #         getcluste(countMtrFiltExp,paste0(outdir,'ALLEXPvsCON_',method),method)
    #     }
    # }else{
        # expCount=1
        # for(eachExp in expGrp){
    # countMtrFiltExp=mtr[samples]
    for(method in cmethods){
        getcluste(mtr,samples,paste0(outdir,prefix,'cluster',ifelse(length(cmethods)==1,'',paste0('_',method))),method)
        # expCount=expCount+1
    }
        # }
    # }
}

getcluste=function(countMtrFilt,samples,prefix,method,expIndex=NA){
    library(pheatmap)
    library(dynamicTreeCut)
    countMtrFiltExp=countMtrFilt[samples]
    # expNum=dim(countMtrFiltExp)[2]-1
    countMtrFiltLog=log1p(countMtrFiltExp)
    # if(is.na(expIndex)){
    #     annotDF=data.frame(Group=c('Con',paste0('Exp',1:expNum), row.names = colnames(countMtrFiltLog)))
    # }else{
    #     annotDF=data.frame(Group=c('Con',paste0('Exp',expIndex), row.names = colnames(countMtrFiltLog)))
    # }
    cat(paste0('Cluste by method:  ',method,'\n'))
    dist.mtr=dist(countMtrFiltLog)
    my.tree=hclust(dist.mtr, method = method)
    cluster=tryCatch(
        cutreeDynamic(my.tree, distM=as.matrix(dist.mtr),verbose=0),
        error=function(e){cat(paste('error in method',method,e,'\n'));return(NULL)})
    if(is.null(cluster)) return(NULL)
    # print(class(cluster))
    cluster.factor=as.factor(cluster)
    # cat(paste0(method,'\n'))
    # cat(levels(cluster.factor))
    # cat('\n')
    annot_row=data.frame(Cluster=cluster.factor,row.names=rownames(countMtrFiltLog))
    annot_color=rainbow(length(levels(cluster.factor)))
    # cat(annot_color)
    # cat('\n')
    names(annot_color)=levels(cluster.factor)
    # groupColor=c('lightblue',rainbow(expNum))
    # names(groupColor)=c('Con',paste0('Exp',1:expNum))
    annot_color=list(Cluster=annot_color)#,
                        # Group=groupColor)

    countMtrFilt$Cluster=cluster
    write.csv(countMtrFilt,file=paste0(prefix,'.csv'))

    pheatmap(countMtrFiltLog,
            show_rownames=F,
            cluster_cols=F,
            cluster_rows=my.tree,
            filename=paste0(prefix,'.png'),
            angle_col=45,
            border=F,
            annotation_row=annot_row,
            annotation_colors=annot_color,
            # annotation_col=annotDF,
            color=colorRampPalette(c('navy','white','firebrick3'))(500),
            height=7,width=6)
}

spec = matrix(c(
    'outdir', 'o', 1, 'character', 'output directory',
    'infile', 'i', 1, 'character', 'input file',
    'samples','s',1,'character','sample colname, splited by \',\'',
    'prefix','n',1,'character','output file prefix',
    'minabslogfc', 'm', 1, 'double', 'min absolute value log2fc (1)',
    # 'congrp','c',1,'character','control colname',
    # 'expgrp','e',1,'character','experiment colname. Splited by \',\'',
    'cmethods','M',1,'character','cluste method. Splited by \',\' (average)',
    # 'allinone','A',0,'logical','cluste all exp group in one graph (disable)',
    'help', 'h', 0, 'logical', 'print this message',
    'version', 'V', 0, 'logical', 'print version message'
), byrow=TRUE, ncol=5)

defaultValue=list(
    outdir=NA,
    infile=NA,
    minabslogfc=1,
    samples=NA,
    prefix='',
    # congrp=NA,
    # expgrp=NA,
    cmethods='average',
    # allinone=FALSE,
    help=NA,
    version=NA
)
isrequired=c('outdir','infile','samples')#'congrp','expgrp')

library(getopt)
argv=commandArgs(T)

opts=getopt(spec,opt=argv)
if (!is.null(opts$version)) {
    cat(paste0('docluster\nVersion:  ',version,'\n'))
    q()
}else if(length(argv)==0 | !is.null(opts$help)){
    cat(getopt(spec,usage=TRUE,command="docluster"))
    q()
}else if (sum(isrequired %in% names(opts))!=length(isrequired)) {
    cat(paste0(isrequired[!(isrequired %in% names(opts))],' is required.\n'))
    # cat(getopt(spec,usage=TRUE,command="docluster"))
    q()
}

for(opt in names(opts)){
    defaultValue[[opt]]=opts[[opt]]
}
opts=defaultValue

outdir=opts$outdir
if(!endsWith(outdir,'/')){
    outdir=paste0(outdir,'/')
}
if(!dir.exists(outdir)){
    dir.create(outdir)
}
infile=opts$infile
minAbsLogFc=opts$minabslogfc
samples=opts$samples
samples=unlist(strsplit(samples,','))
prefix=opts$prefix
if(prefix!=''){
    prefix=paste0(prefix,'_')
}
# conGrp=opts$congrp
# expGrp=opts$expgrp
# expGrp=unlist(strsplit(expGrp,','))
cmethods=opts$cmethods
cmethods=unlist(strsplit(cmethods,','))
# allinone=opts$allinone
countMtr=read.csv(infile,header=T,row.names=1)
# row.names(countMtr)=countMtr$Row.names
dim(countMtr)
# if(sum(c(conGrp,expGrp) %in% colnames(countMtr))!=length(c(conGrp,expGrp))){
#     cat(paste0('Cannot find column named ',c(conGrp,expGrp)[c(conGrp,expGrp) %in% colnames(countMtr)],'\n'))
#     q()
# }
if(sum(samples %in% colnames(countMtr))!=length(samples)){
    cat(paste0('Cannot find column named ',samples[!(samples %in% colnames(countMtr))],'\n'))
    q()
}

countMtrFilt=countMtr[!is.na(countMtr$log2FoldChange) & abs(countMtr$log2FoldChange)>=minAbsLogFc,]
dim(countMtrFilt)

doCluste(countMtrFilt,samples,outdir=outdir,cmethods=cmethods,prefix=prefix)
