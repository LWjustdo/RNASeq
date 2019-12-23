#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/12/13 
@Description :
'''
import  os,sys,time
from profile import Profile,timefly
from argparse import ArgumentParser
from multiprocessing import Pool

var_path = Profile()

@timefly
def GeneBodyCoverage(prefix,species):
    path = 'BamQC/'+prefix
    os.makedirs(path,exist_ok=True)
    os.system("{RSeQC_path}/geneBody_coverage.py -r {ref_path}/bed12/{s}.bed12 -i Mapping/{p}_Aligned.sortedByCoord.out.bam -o {pa}/{p} -f png "
              ">>BamQC/rseqc.log 2>&1".format(s=species,pa=path,p=prefix,**var_path))

@timefly
def InnerDistance(prefix,species):
    path = 'BamQC/' + prefix
    os.makedirs(path, exist_ok=True)
    os.system("{RSeQC_path}/inner_distance.py -r {ref_path}/bed12/{s}.bed12 -i Mapping/{p}_Aligned.sortedByCoord.out.bam -o {pa}/{p} "
              ">>BamQC/rseqc.log 2>&1".format(s=species,pa=path,p=prefix,**var_path))

@timefly
def JunctionSaturation(prefix, species):
    path = 'BamQC/' + prefix
    os.makedirs(path, exist_ok=True)
    os.system("{RSeQC_path}/junction_saturation.py -r {ref_path}/bed12/{s}.bed12 -i Mapping/{p}_Aligned.sortedByCoord.out.bam -o {pa}/{p} "
              ">>BamQC/rseqc.log 2>&1".format(s=species,pa=path, p=prefix, **var_path))

@timefly
def ReadDistribution(prefix, species):
    path = 'BamQC/' + prefix
    os.makedirs(path, exist_ok=True)
    os.system("{RSeQC_path}/read_distribution.py -r {ref_path}/bed12/{s}.bed12 -i Mapping/{p}_Aligned.sortedByCoord.out.bam > {pa}/{p}.featureDistribution.xls "
              ">>BamQC/rseqc.log 2>&1".format(s=species,pa=path, p=prefix, **var_path))

@timefly
def RPKMsaturation(prefix, species):
    path = 'BamQC/' + prefix
    os.makedirs(path, exist_ok=True)
    os.system("{RSeQC_path}/RPKM_saturation.py -r {ref_path}/bed12/{s}.bed12 -i Mapping/{p}_Aligned.sortedByCoord.out.bam -o {pa}/{p} -d '1++,1--,2+-,2-+' "
              ">>BamQC/rseqc.log 2>&1".format(s=species,pa=path, p=prefix, **var_path))

@timefly
def Bam2wg(prefix,threads):
    os.system("source /home/liang.wan/anaconda3/bin/activate /home/liang.wan/anaconda3/envs/scasat && "
              "bamCoverage -p {t} -b Mapping/{p}_Aligned.sortedByCoord.out.bam -o BamQC/{p}.bw >>BamQC/correlate.log 2>&1".format(p=prefix,t=threads, **var_path))

@timefly
def PlotCorrelation(allSamps,threads):
    path = 'BamQC'
    os.makedirs(path, exist_ok=True)
    bws =' '.join( ['BamQC/'+prefix+'.bw' for prefix in allSamps])
    num = len(allSamps)
    pool = Pool(num)
    for pre in allSamps:
        pool.apply_async(Bam2wg,args=(pre,threads,))
    pool.close()
    pool.join()
    os.system(" source /home/liang.wan/anaconda3/bin/activate /home/liang.wan/anaconda3/envs/scasat && "
              "multiBigwigSummary bins -p {t} -b {b} -o BamQC/bw.npz >>BamQC/correlate.log 2>&1 && "
              "plotCorrelation -in BamQC/bw.npz -c pearson -p heatmap -o BamQC/sample_correlation.png >>BamQC/correlate.log 2>&1 && "
              "source /home/liang.wan/anaconda3/bin/deactivate".format(b=bws,t=threads))


def BamQC(prefix,species):
    pool =Pool(5)
    todolist = [GeneBodyCoverage,InnerDistance,JunctionSaturation,ReadDistribution,RPKMsaturation]
    for i in todolist:
        pool.apply_async(i,args=(prefix,species,))
    pool.close()
    pool.join()

def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--model', dest='model', type=str, default='', help='输入运行功能[BamQC,PlotCorrelation]，必须！')
    parser.add_argument('--prefix', dest='prefix', type=str, default='', help='[BamQC] 输入样本名称，需要！')
    parser.add_argument('--control_samples', dest='con_samp', type=str, default='', help='[PlotCorrelation] 对照组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--experimental_samples', dest='exp_samp', type=str, default='', help='[PlotCorrelation] 室验组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--species', dest='species', type=str, default='GRCh37', help='[BamQC] 输入物种名称[GRCh37,GRCm38]，默认GRCh37！')
    parser.add_argument('--threads', dest='threads', type=int, default=16, help='输入配置cpu线程数，默认16 !')
    return parser

if __name__  == '__main__':
    if len(sys.argv) < 2:
        print('\nusage:  python {} -h \n'.format(sys.argv[0]))
        sys.exit(1)
    parser = Parser_opt()
    args = parser.parse_args()
    print('\n' + '#' * 60)
    print('python ' + ' '.join(sys.argv))
    print('{0} 开始分析......'.format(time.ctime()))
    if args.model == 'BamQC':
        BamQC(args.prefix,args. species)
    else:
        conSamps = args.con_samp.split(',')
        expSamps = args.exp_samp.split(',')
        allSamps = conSamps + expSamps
        PlotCorrelation(allSamps, args.threads)
