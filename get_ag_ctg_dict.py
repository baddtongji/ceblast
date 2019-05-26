# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:00:38 2015

@author: Administrator
"""

"""read ag_ctg from file"""
import pandas as pd
from django.utils.datastructures import SortedDict
df = pd.read_table(r'/home/djangoapps/work/mysite/epiblast_program/ag_ctg_dict.txt')
ag_ctg_dict = SortedDict()
rev_ag_ctg_dict = SortedDict()
for idx, ctg_no in enumerate(df['category']):
    ctg_name = df['name'][idx]
    ag_ctg_dict[ctg_name] = ctg_no
    rev_ag_ctg_dict[ctg_no] = ctg_name
#print ag_ctg_dict, len(ag_ctg_dict)
