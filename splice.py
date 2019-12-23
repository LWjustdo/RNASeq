#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/12/4 
@Description :
无重复样本使用rMATS3版本，重复样本使用rMATS4版本
'''
import os
import time
import sys
from profile import Profile,timefly
from argparse import ArgumentParser
import pandas as pd
from multiprocessing import Process,Pool
import shutil

var_path = Profile()

@timefly
def asprofile(prefix,species):
    if not os.path.exists('Splice'):
        os.mkdir('Splice')
    os.system("{ASprofile_path}/extract-as  Assembly/{T}.gtf {ref_path}/fasta/{s}/{s}.fa.hdrs > Splice/{T}.as.txt "
              ">Splice/{T}.asprofile.log 2>&1".format(T=prefix,s=species, **var_path))
    os.system("perl {ASprofile_path}/summarize_as.pl Assembly/{T}.gtf Splice/{T}.as.txt -p Splice/{T} "
              ">>Splice/{T}.asprofile.log 2>&1".format(T=prefix, **var_path))
    os.system("{ASprofile_path}/extract-as-fpkm Assembly/{T}.gtf {ref_path}/fasta/{s}/{s}.fa.hdrs Splice/{T}.as.nr > Splice/{T}.as.fpkm.txt "
              ">>Splice/{T}.asprofile.log 2>&1".format(T=prefix,s=species, **var_path))

@timefly
def rMATS(conSamps,expSamps,species):
    os.makedirs('Splice',exist_ok=True)
    if len(conSamps) == 1:
        prefix1,prefix2 = conSamps+expSamps
        os.system("source /home/xueqiang.liu/anaconda3/bin/activate python27 && "
              "python {rMATS3} -b1 Mapping/{p1}.sort.bam -b2 Mapping/{p2}.sort.bam -c 0.05 -gtf {ref_path}/gtf/{s}/{s}.gtf -o Splice -t paired -len 150 >>Splice/splice.log 2>&1 && "
              "source /home/xueqiang.liu/anaconda3/bin/deactivate".format(p1=prefix1,p2=prefix2,s=species,  **var_path))
        try:
            shutil.rmtree('Splice/SAMPLE_1')
            shutil.rmtree('Splice/SAMPLE_2')
        except:
            print('Splice/SAMPLE_1文件不存在')
    else:
        conBam = open('Mapping/conBam.txt', 'w', encoding='utf-8')
        info = ','.join(["Mapping/"+consamp+".sort.bam" for consamp in conSamps])
        conBam.write(info)
        conBam.close()
        expBam = open('Mapping/expBam.txt', 'w', encoding='utf-8')
        info2 = ','.join(["Mapping/" + expsamp + ".sort.bam" for expsamp in expSamps])
        expBam.write(info2)
        expBam.close()
        os.system("source /home/xueqiang.liu/anaconda3/bin/activate python27 && "
                  "python {rMATS4} --b1 Mapping/conBam.txt --b2 Mapping/expBam.txt --gtf {ref_path}/gtf/{s}/{s}.gtf --od Splice -t paired --readLength 150 >>Splice/splice.log 2>&1 && "
                  "source /home/xueqiang.liu/anaconda3/bin/deactivate".format( s=species,**var_path))

def rmatsPlot(conBams,expBams,asevent):
    os.system("source /home/xueqiang.liu/anaconda3/bin/activate python27 &&  "
              "python {rmatsPlot} --b1 {cb} --b2 {eb} --l1 control --l2 experimental --exon_s 1 --intron_s 5 -t {ase} "
              "-e Splice/{ase}.MATS.JC.filtered.txt -o Splice/{ase} >>Splice/plot.log 2>&1 && "
              "source /home/xueqiang.liu/anaconda3/bin/deactivate".format(cb=conBams,eb=expBams,ase=asevent,** var_path))

@timefly
def rmats2sashimiplot(conSamps, expSamps):
    ASevents=['A3SS','A5SS','MXE','RI','SE']
    for ase in ASevents:
        asfile = 'Splice/MATS_output/' +ase + '.MATS.JunctionCountOnly.txt' if len(conSamps) == 1 else 'Splice/'+ase+'.MATS.JC.txt' #样本是否重复文件路径、名称均不同
        df = pd.read_csv(asfile,sep='\t')
        df2 = df[ abs(df.IncLevelDifference) >= 0.5]  if len(conSamps) == 1 else df[(df.FDR <= 0.05) & (abs(df.IncLevelDifference) >= 0.5)] #提取P值小于0.05，差异大于50%的数据,无重复不考虑P值
        df2['absIncLD'] = abs(df2['IncLevelDifference'])
        df3 = df2.sort_values('absIncLD',ascending=False)
        df3 = df3.drop('absIncLD',axis=1)
        df3.to_csv('Splice/{}.MATS.JC.filtered.txt'.format(ase),sep='\t',index=False)
    conBams = ','.join(["Mapping/" + consamp + ".sort.bam" for consamp in conSamps])
    expBams = ','.join(["Mapping/" + expsamp + ".sort.bam" for expsamp in expSamps])
    pool=Pool(5)
    for asevent in ASevents:
        pool.apply_async(rmatsPlot,args=(conBams,expBams,asevent,))
    pool.close()
    pool.join()



@timefly
def cuffcompare(prefix,species):
    os.system("{cuffcompare} -r {ref_path}/gtf/{s}/{s}.gtf -s {ref_path}/fasta/{s}/{s}.fa -o Splice/{T} Assembly/{T}.gtf "
              ">Splice/{T}.cuffcompare.log 2>&1" .format(T=prefix,s=species, **var_path))


def Parser_opt():
    parser = ArgumentParser()
    parser.add_argument('--control_samples', dest='con_samp', type=str, default='', help='对照组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--experimental_samples', dest='exp_samp', type=str, default='', help='室验组样本名称，如果含重复用逗号分隔，必需！')
    parser.add_argument('--model', dest='model', type=str, default='', help='输入运行模块[rMATS,rmats2sashimiplot]，必须！')
    parser.add_argument('--species', dest='species', type=str, default='GRCh37', help='输入物种名称[GRCh37,GRCm38]，默认GRCh37！')
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
    expSamps = args.exp_samp.split(',')
    if args.model == 'rMATS':
        rMATS(conSamps, expSamps, args.species)
    else:
        rmats2sashimiplot(conSamps, expSamps)