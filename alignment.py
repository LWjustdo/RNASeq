#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/1 
@Description :比对，排序，去重，碱基质量控制
'''

import os
import time
import sys
from profile import Profile,timefly
from argparse import ArgumentParser

var_path = Profile()

@timefly
def Star(prefix,species,threads):
    os.makedirs('Mapping',exist_ok=True)
    os.system("{STAR} --runThreadN {th} --genomeDir {ref_path}/star_index/{s} --readFilesCommand zcat "
              "--readFilesIn DataQC/{T}_R1.clean.fastq.gz DataQC/{T}_R2.clean.fastq.gz "
              "--twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 "
              "--alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 "
              "--alignSJstitchMismatchNmax 5 -1 5 5 --outSAMstrandField intronMotif "
              "--outFileNamePrefix Mapping/{T}_  --outSAMtype BAM SortedByCoordinate "
              "--outBAMsortingThreadN {th}  --chimOutJunctionFormat 1 " #--quantMode TranscriptomeSAM GeneCounts
              ">Mapping/{T}.star.log 2>&1".format(T=prefix, th=threads, s=species,**var_path))
    os.system("{samtools} index Mapping/{T}_Aligned.sortedByCoord.out.bam".format(T=prefix, **var_path))
    current_path = os.getcwd()
    os.system("ln -s {cp}/Mapping/{T}_Aligned.sortedByCoord.out.bam Mapping/{T}.sort.bam".format(cp=current_path,T=prefix))

@timefly
def Hisat2(prefix,species,threads):
    os.makedirs('Mapping', exist_ok=True)
    os.system("{hisat2} --dta -t -p {th} -x {ref_path}/hisat2_index/{s}/genome -1 DataQC/{T}_R1.clean.fastq.gz "
        " -2 DataQC/{T}_R2.clean.fastq.gz -S Mapping/{T}.sam  >Mapping/{T}.hisat2.log 2>&1".format(T=prefix, th=threads,s=species, **var_path))
    os.system("{samtools} view -@ {th} -bS Mapping/{T}.sam -o Mapping/{T}.bam >>Mapping/{T}.hisat2.log 2>&1".format(T=prefix,th=threads,**var_path))
    os.system("{samtools} sort -@ {th} Mapping/{T}.bam Mapping/{T}.sort >>Mapping/{T}.hisat2.log 2>&1".format(T=prefix,th=threads,**var_path))
    os.system("rm Mapping/{T}.sam Mapping/{T}.bam ".format(T=prefix,  **var_path))

def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--prefix', dest='prefix', type=str, default='', help='输入样本名称，必须！')
    parser.add_argument('--species', dest='species', type=str, default='GRCh37', help='输入物种名称[GRCh37,GRCm38]，默认GRCh37！')
    parser.add_argument('--alignTool', dest='alignTool', type=str, default='star', help='输入软件名称[star,hisat2]，默认star！')
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
    print('{0} 开始对样本{1}进行分析......'.format(time.ctime(), args.prefix))
    if args.alignTool == 'star':
        Star(args.prefix,args.species,args.threads)
    else:
        Hisat2(args.prefix,args.species,args.threads)



