# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:02:14 2015

@author: Administrator
"""
import os
from Bio.PDB.PDBParser import PDBParser
from get_fp import Complex

def get_fp_files(flu_file, fp_out_path, epi_resi_list=None, horizontal=5, vertical=4, physchem_type_num=3):
    p = PDBParser(PERMISSIVE=1)
#    print 'pdb file', flu_file
    epi_resi_list = map(int, epi_resi_list)
    struct = p.get_structure(os.path.basename(flu_file), flu_file)
    c = Complex(struct)
    if not os.path.exists(fp_out_path):
#        print 'no',fp_out_path
        #physchem_type_num = 4表示要糖基化;3表示不要糖基化
        c.get_fp(spin_image_radius_step = horizontal, spin_image_height_step = vertical, physchem_type_num = physchem_type_num)
        #epi_resi_list is a int list
        if not epi_resi_list is None:
            c.write_fp_to_file_with_given_res(fp_out_path, epi_resi_list)
        else:
            c.write_fp_to_file(fp_out_path)
#    else:
#        print 'have', fp_out_path
