# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:46:52 2015

@author: Administrator
"""

"""
This program is to get glycosylated site
"""

def get_glyco_site(site, next_site1, next_site2):
    if site == 'N' and next_site1 != 'P' and next_site2 in ['S', 'P']:
        return 1
    return 0
    