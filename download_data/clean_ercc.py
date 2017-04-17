#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np

pd.read_table(sys.stdin, 
            names = ['number','id','group','mix1','mix2','fold','log2fold'],
            skiprows=1) \
    .drop('number', axis=1) \
    .assign(label = lambda d: np.where(d.log2fold==0, 'DE','notDE'))\
    .to_csv(sys.stdout,index=False, sep='\t')
