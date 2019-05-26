# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:57:25 2015

@author: Administrator
"""
"""
因为只能有一条抗原链，对于多条要进行合并
"""
from basic_pdb_related_func import preprocessFile, getPdbStruct, writePdbStruct2File

def covertMoreChn2SingleChn(pdb_path, chn_lst, merge_chn_name, out_path):
    """
    chn_lst: chain list to be merged
    """
    pdb_path = preprocessFile(pdb_path)
    struct = getPdbStruct(pdb_path)

    #change residue id
    cur_resi_id = 0
    for chn in struct.get_chains():
        if chn.get_id() in chn_lst and chn.get_id() != merge_chn_name:
            for res in chn.get_residues():
                cur_resi_id += 1
                res.id = (res.id[0], cur_resi_id, res.id[2])
            chn.id = merge_chn_name
    #save to file
    writePdbStruct2File(struct, out_path)

if __name__ == '__main__':
    pdb_path = r'/pywork/aaa_uniq_id.pdb'
    chn_lst = ['A', 'B']
    merge_chn_name = 'X'
    out_path = r'/pywork/aaa_std.pdb'
    covertMoreChn2SingleChn(pdb_path, chn_lst, merge_chn_name, out_path)
    