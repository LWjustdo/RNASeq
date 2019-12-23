#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/12/12 
@Description :差异分析/火山图/差异注释
'''

import os,sys,time
from profile import timefly
from collections import defaultdict
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from argparse import ArgumentParser


@timefly
def Diff(repeatNum,species):
    os.makedirs('DiffExp',exist_ok=True)
    path = os.path.dirname(__file__)
    rscript = path +'/diff.R'
    if repeatNum == 1:
        os.system("Rscript "+rscript+" gene "+str(repeatNum)+" "+species+" >>DiffExp/deseq2.log 2>&1")
        os.system("Rscript " + rscript + " transcript " + str(repeatNum) + " " + species+" >>DiffExp/deseq2.log 2>&1")
    else:
        os.system("source /home/xueqiang.liu/anaconda3/bin/activate python37 && "+
                 "Rscript "+rscript+" gene "+str(repeatNum)+" "+species +" >>DiffExp/deseq2.log 2>&1 && "+
                 "Rscript " + rscript + " transcript " + str(repeatNum) + " " + species +
                  " >>DiffExp/deseq2.log 2>&1 && source /home/xueqiang.liu/anaconda3/bin/deactivate")
    rscript2 = path + '/bimoRt.R'
    try:
        os.system("Rscript "+rscript2+" gene "+species+" >>DiffExp/bimoRt.log 2>&1")
    except:
        print('ERROR: please rerun Rscript path/to/bimoRt.R gene [species]\n')
    try:
        os.system("Rscript " + rscript2 + " transcript " + species+" >>DiffExp/bimoRt.log 2>&1")
    except:
        print('ERROR: please rerun Rscript path/to/bimoRt.R transcript  [species]\n')



def DiffTransFig():
    file = 'DiffExp/transcript_diff_results.csv'
    trans_type_dict = defaultdict(int)
    lncRNA_length_dict =defaultdict(int)
    with open(file) as f:
        for line in f:
            if not line.startswith('transcript_id'):
                lin = line.strip().split(',')
                type = lin[3].strip('"')
                trans_type_dict[type] += 1 #转录本类型及对应数量
                if type =='lncRNA':
                    lncRNA_length_dict[int(lin[4])] += 1#lncRNA长度及对应数量

    #转录本类型饼图
    total = sum([trans_type_dict[k] for k in trans_type_dict.keys()])
    typelist=['protein_coding','lncRNA','miRNA','retained_intron']
    num = [trans_type_dict[k] for k in typelist]
    other = total - sum(num)
    typelist.append('other')
    num.append(other)
    explode =(0.1,0,0,0,0)
    fig1,ax1 = plt.subplots()
    ax1.pie(num,explode=explode,labels=typelist,autopct='%1.1f%%',startangle=90)
    ax1.axis('equal')
    # plt.show()
    plt.savefig('DiffExp/trans_type.png')

    #lncRNA长度分布
    lncx =sorted([k for k in lncRNA_length_dict.keys()])
    lncy = [lncRNA_length_dict[k] for k in lncx ]
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2.set_xlabel('lncRNA length')
    ax2.set_ylabel('counts')
    plt.sca(ax2)
    plt.plot(lncx,lncy)
    # plt.show()

    # lncRNA累计长度分布
    cum_length_dict = defaultdict(int)
    max_length = max([x for x in lncRNA_length_dict])
    all_num = sum(lncy)

    for n in range(max_length):  # 从0到最大长度，补充测序深度出现次数为0的数据
        if not n in lncRNA_length_dict:
            lncRNA_length_dict[n] = 0

    for n in range(max_length):  # 累积测序深度-对应read数量
        for m in range(max_length):
            if m > n:
                cum_length_dict[n] += lncRNA_length_dict[m]
    lncx2 = sorted([k for k in cum_length_dict])
    lncy2 = [cum_length_dict[k] for k in lncx2]
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.set_xlabel('lncRNA cum length')
    ax3.set_ylabel('counts')
    plt.sca(ax3)
    plt.plot(lncx2, lncy2)
    # plt.show()
    fig2.savefig( 'DiffExp/lncRNA_length.png')
    fig3.savefig( 'DiffExp/lncRNA_cum_length.png')


def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--model', dest='model', type=str, default='', help='输入运行功能[Diff,DiffTransFig]，必须！')
    parser.add_argument('--control_samples', dest='con_samp', type=str, default='', help='[Diff] 对照组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--experimental_samples', dest='exp_samp', type=str, default='', help='[Diff] 室验组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--species', dest='species', type=str, default='GRCh37', help='[Diff] 输入物种名称[GRCh37,GRCm38]，默认GRCh37！')
    return parser


if __name__  == '__main__':
    if len(sys.argv) < 2:
        print('\nusage:  python {} -h \n'.format(sys.argv[0]))
        sys.exit(1)
    parser = Parser_opt()
    args = parser.parse_args()
    print('\n' + '#' * 60)
    print('python ' + ' '.join(sys.argv))
    print('{0} 开始进行分析......'.format(time.ctime()))
    conSamps = args.con_samp.split(',')
    if args.model == 'Diff':
        Diff(len(conSamps), args.species)
    else:
        DiffTransFig()

