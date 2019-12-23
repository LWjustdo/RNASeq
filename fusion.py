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

var_path = Profile()

@timefly
def starFusion(prefix,species,threads):
    if not os.path.exists('Fusion'):
        os.mkdir('Fusion')
    os.system("{STAR-Fusion} --genome_lib_dir {ref_path}/starFusion_index/{s}/ctat_genome_lib_build_dir "
              "--CPU {th} -J Mapping/{T}_Chimeric.out.junction --output_dir Fusion/{T} "
              ">>Fusion/fusion.log 2>&1" .format(T=prefix, s=species,th=threads,**var_path))



if __name__  == '__main__':
    if len(sys.argv) < 3:
        print('\nusage:  python {} [prefix] [species:GRCh37,GRCm38] [threads]\n'.format(sys.argv[0]))
        sys.exit(1)
    prefix = sys.argv[1]
    species = sys.argv[2]
    threads = sys.argv[3]
    starFusion(prefix,species, threads)
