#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/6/24 
@Description :fastq文件质控
'''

import os
import time
import sys
from profile import Profile,timefly

var_path = Profile()


@timefly
def QC(prefix, threads):
    path = 'DataQC'
    os.makedirs(path,exist_ok=True)
    fastq = [ i for i in os.listdir() if prefix in i]
    fastq1 = ''.join([i for i in fastq if '1.f' in i])
    fastq2 = ''.join([i for i in fastq if '2.f' in i])
    os.system("{fastp} -w {th} --in1 {f1} --out1 {p}/{T}_R1.clean.fastq.gz "
              "--in2 {f2} --out2 {p}/{T}_R2.clean.fastq.gz "
              "--low_complexity_filter --correction --length_required=70 "
              "--html {p}/{T}.QCReport.html --json {p}/{T}.json --report_title {p}/{T}.QCReport "
              " >{p}/{T}.fastp.log 2>&1".format(T=prefix, p=path,f1=fastq1, f2=fastq2, th=threads,**var_path))
    os.system("python {summary4fastp} {p}/{T}.json > {p}/{T}.QCsummary.xls ".format(T=prefix,p=path,**var_path))


if __name__  == '__main__':
    if len(sys.argv) < 2:
        print('\nusage:  python {} [prefix] [threads]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    threads = sys.argv[2]
    QC(prefix, threads)

