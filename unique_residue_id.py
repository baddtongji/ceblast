# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:59:22 2015

@author: Administrator
"""
#remove resdiues' repeated id, e.g. 100C, 100D, 100E etc
#to get unique residues' id like 100, 101, 102 etc
"""
将pdb文件中不确定的residue id变成唯一id号，以使得用户在输入id后与我们的程序输入相一致
"""

from basic_pdb_related_func import preprocessFile, getPdbStruct, writePdbStruct2File

def unique_residue_id(pdb_path, out_path):
    pdb_path = preprocessFile(pdb_path)
    struct = getPdbStruct(pdb_path)
    
    for chn in struct.get_chains():
        pre_resi_id = -1
        for res in chn.get_residues():
            cur_resi_id = res.id[1]
            #first resi id
            if pre_resi_id == -1:
               pre_resi_id = cur_resi_id
               continue
            #rest resi id
            if cur_resi_id <= pre_resi_id:
                res.id = (res.id[0], pre_resi_id+1, ' ')
                pre_resi_id += 1
            else:
                pre_resi_id = cur_resi_id
            
    writePdbStruct2File(struct, out_path)

if __name__ == '__main__':
    pdb_path = r'/pywork/aaa.pdb'
    out_path = r'/pywork/aaa_uniq_id.pdb'
    unique_residue_id(pdb_path, out_path)