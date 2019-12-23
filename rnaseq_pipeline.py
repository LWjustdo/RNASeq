#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/12/10 
@Description :
'''

from argparse import ArgumentParser
from multiprocessing import Pool,Process
import sys,time,os
from qc import QC
from alignment import Star,Hisat2
from quantify import Stringtie_single,Matrix
from diff import Diff,DiffTransFig
from fusion import starFusion
from splice import rMATS,rmats2sashimiplot
from clusterEnrich import Cluster,Enrich
from bamqc import BamQC,PlotCorrelation



def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--control_samples', dest='con_samp', type=str, default='', help='对照组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--experimental_samples',dest='exp_samp',type=str,default='',help='室验组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--species', dest='species', type=str, default='GRCh37', help='物种名称[GRCh37,GRCm38]，默认GRCh37！')
    parser.add_argument('--threads', dest='threads', type=int, default=16, help='配置cpu线程数，默认16 !')
    parser.add_argument('--alignTool', dest='alignTool', type=str, default='star', help='输入软件名称[star,hisat2]，默认star！')

    return parser

def basicAnalyze(sample,species,threads,alignTool):
    QC(sample,threads)
    if alignTool == 'star':
        Star(sample,species,threads)
    else:
        Hisat2(sample,species,threads)
    Stringtie_single(sample,species,threads)
    starFusion(sample, species, threads)
    BamQC(sample, species)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('\nUsage:   python {} -h \n'.format(sys.argv[0]))
        sys.exit(1)
    parser = Parser_opt()
    args = parser.parse_args()
    conSamps = args.con_samp.split(',')
    expSamps = args.exp_samp.split(',')
    allSamps = conSamps + expSamps
    print('{0} 开始进行分析......'.format(time.ctime()))
    print('python ' + ' '.join(sys.argv))
    if  len(conSamps) != len(expSamps):
        print('两组样本数量不一致，请统一！')
        sys.exit(1)
    proc_list = []
    for samp in allSamps:
        proc = Process(target=basicAnalyze, args=(samp, args.species, args.threads, args.alignTool,))
        proc_list.append(proc)
    for p in proc_list:
        p.start()
    for p in proc_list:
        p.join()
    PlotCorrelation(allSamps,args.threads)
    Matrix(allSamps)
    try:
        Diff(len(conSamps),args.species)
        DiffTransFig()
    except:
        print('Diff error !')
    finally:
        Cluster(allSamps)
        Enrich()
    rMATS(conSamps,expSamps,args.species)
    rmats2sashimiplot(conSamps,expSamps)


    print('分析结束')
    t1 = time.time()
    print('\n' + '#' * 60)

