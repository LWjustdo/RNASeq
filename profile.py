#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@Author      : xueqiang.liu
@contact     : xqliu@dragondropinn.com
@Date        : 2019/7/1 
@Description : RNAseq   Profile.txt文件
'''
import os
import functools
import time

#脚本所在路径
path = os.path.dirname(__file__)
file = path +'/Profile.txt'

def Profile():
    dic = {}
    with open(file) as f:
        for line in f:
            lin = line.strip().split(':')
            if len(lin) == 2:
                dic[lin[0]] = lin[1]
    return dic

def timefly(func):
    @functools.wraps(func)
    def wrapper(*args,**kw):
        s=time.time()
        res = func(*args,**kw)
        e=time.time()
        print('{} runtime: {}'.format(func.__name__,TransTime(e-s)))
        return res
    return wrapper

def TransTime(seconds):
    h = seconds//3600
    m = seconds%3600//60
    s = seconds%60
    return '{}h {}min {:.0f}s'.format(h,m,s)


if __name__ == '__main__':
    var = Profile()
    print(var)
