#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/11/29 
@Description :
'''
import os
import time
import sys
from profile import Profile,timefly
from argparse import ArgumentParser


var_path = Profile()

@timefly
def stringtie1(prefix,species,threads):
    os.makedirs('Assembly',exist_ok=True)
    os.system("{stringtie} -p {th} -G {ref_path}/gtf/{s}/{s}.gtf  -o Assembly/{T}.gtf  Mapping/{T}.sort.bam "
              ">Assembly/{T}.stringtie.log 2>&1".format(T=prefix, th=threads,s=species,**var_path))

@timefly
def stringtie_merge(species,threads):
    if not os.path.exists('Assembly'):
        print("请确认是否已对每个样本单独运行stringtie！")
    os.system("ls Assembly/*gtf >Assembly/gtflist.txt")
    os.system("{stringtie} --merge -p {th} -G {ref_path}/gtf/{s}/{s}.gtf  -o Assembly/stringtie_merged.gtf Assembly/gtflist.txt "
              ">>Assembly/merge.stringtie.log 2>&1".format(th=threads,s=species,**var_path))

@timefly
def stringtie2(prefix,threads):
    os.system("{stringtie} -p {th} -G Assembly/stringtie_merged.gtf -e -B -o ballgown/{T}/{T}.gtf  Mapping/{T}.sort.bam "
              ">>Assembly/{T}.stringtie.log 2>&1".format(T=prefix, th=threads, **var_path))

@timefly
def Stringtie_single(prefix,species,threads):
    os.makedirs('Quantification',exist_ok=True)
    os.system("{stringtie} -p {th} -G {ref_path}/gtf/{s}/{s}.gtf -e -B -o Quantification/{T}/{T}.gtf  Mapping/{T}.sort.bam "
              ">>Quantification/stringtie.log 2>&1".format(T=prefix, th=threads,s=species, **var_path))

@timefly
def Matrix(sampsList):
    sample_list_file = open('Quantification/sample_list.txt','w',encoding='utf-8')
    for samp in sampsList:
        sample_list_file.write("{0} Quantification/{0}/{0}.gtf\n".format(samp))
    sample_list_file.close()
    os.system("source /home/xueqiang.liu/anaconda3/bin/activate python27 && "
              "python {prepDE}  -i Quantification/sample_list.txt -g Quantification/gene_count_matrix.csv -t Quantification/transcript_count_matrix.csv && "
              "source /home/xueqiang.liu/anaconda3/bin/deactivate".format( **var_path))


def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--model', dest='model', type=str, default='', help='输入运行功能[stringtie1,stringtie_merge,stringtie2,stringtie_single,matrix]，必须！')
    parser.add_argument('--prefix', dest='prefix', type=str, default='', help='输入样本名称，(stringtie1,stringtie2,stringtie_single)需要！')
    parser.add_argument('--control_samples', dest='con_samp', type=str, default='', help='[matrix] 对照组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--experimental_samples',dest='exp_samp',type=str,default='',help='[matrix] 室验组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--species', dest='species', type=str, default='GRCh37', help='输入物种名称[GRCh37,GRCm38]，(stringtie1,stringtie_merge,stringtie_single)需要，默认GRCh37！')
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
    if args.model == 'stringtie1':
        stringtie1(args.prefix,args.species,args.threads)
    elif args.model == 'stringtie_merge':
        stringtie_merge(args.species, args.threads)
    elif args.model == 'stringtie2':
        stringtie2(args.prefix, args.threads)
    elif args.model == 'stringtie_single':
        Stringtie_single(args.prefix, args.species, args.threads)
    else:
        conSamps = args.con_samp.split(',')
        expSamps = args.exp_samp.split(',')
        allSamps = conSamps + expSamps
        Matrix(allSamps)
