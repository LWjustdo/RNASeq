#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/12/13 
@Description :
'''

import os,sys,time
from profile import timefly
from argparse import ArgumentParser

@timefly
def Cluster(sampList):
    os.makedirs('Cluster', exist_ok=True)
    path = os.path.dirname(__file__)
    rscript = path + '/docluster.R'
    allsamp = ','.join(sampList)
    os.system("Rscript "+rscript+" -o Cluster -i DiffExp/gene_rawresult.csv -s {} -n gene >Cluster/cluster.gene.log 2>&1".format(allsamp))
    os.system("Rscript " + rscript + " -o Cluster -i DiffExp/transcript_rawresult.csv  -s {} -n transcript >Cluster/cluster.trans.log 2>&1".format(allsamp))

@timefly
def Enrich():
    os.makedirs('Enrich',exist_ok=True)
    path = os.path.dirname(__file__)
    rscript = path + '/doenrich.R'
    os.system("Rscript " + rscript +" -o Enrich -i DiffExp/gene_rawresult.csv -P -n enrich >Enrich/enrich.log 2>&1")


def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--model', dest='model', type=str, default='', help='输入运行功能[Cluster,Enrich]，必须！')
    parser.add_argument('--control_samples', dest='con_samp', type=str, default='', help='[Cluster] 对照组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--experimental_samples',dest='exp_samp',type=str,default='',help='[Cluster] 室验组样本名称，如果含重复用逗号分隔，必需！')
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
    if args.model == 'Cluster':
        conSamps = args.con_samp.split(',')
        expSamps = args.exp_samp.split(',')
        allSamps = conSamps + expSamps
        Cluster(allSamps)
    else:
        Enrich()
