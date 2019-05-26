# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:05:20 2015

@author: Administrator
"""
"""
根据复合物结构提取抗原表位，只允许指定一条抗原链，并且所有residue的Id要唯一
"""
import os, itertools
from dist import atom_dist
from Bio.PDB import Select
from basic_pdb_related_func import preprocessFile, getPdbStruct, writePdbStruct2File

class ChainSelect(Select):
    #select signle chain
    def __init__(self, chn_lst):
        self.chn_lst = chn_lst
    def accept_chain(self, chain):
        if chain.get_id() in self.chn_lst:
            return 1
        else:
            return 0

class ResiSelect(Select):
    def __init__(self, resi_lst):
        self.resi_lst = resi_lst
    def accept_residue(self, residue):
#        chain_id = residue.get_parent().get_id()
        res_id = str(residue.get_id()[1])
#        res_uni_id = chain_id + '_' + res_id
        res_uni_id = res_id
        if res_uni_id in self.resi_lst:
            return 1
        else:
            return 0

def creatAbAgPdb(pdb_path, ag_chn, ab_chn_lst):
    #读取一个ag_chain和一个ab_chain的pdb信息，存成pdb文件并返回存储路径
    pdb_path = preprocessFile(pdb_path)
    pdb_doc = os.path.dirname(pdb_path)
    pdb_basename = os.path.basename(pdb_path)
    
    struct = getPdbStruct(pdb_path)   
    
    ag_chn_lst = [ag_chn]
    if ab_chn_lst is None:
        ab_chn_lst = [chn.get_id() for chn in struct.get_chains() if chn.get_id()!=ag_chn]
    else:
        ab_chn_lst = ab_chn_lst
    ag_path = os.path.join(pdb_doc, pdb_basename.split('.')[0] + '_' + '_'.join(ag_chn_lst) + '.temp')
    ab_path = os.path.join(pdb_doc, pdb_basename.split('.')[0] + '_' + '_'.join(ab_chn_lst) + '.temp')
    
    getChainFromFile(pdb_path, ag_chn_lst, ag_path)
    getChainFromFile(pdb_path, ab_chn_lst, ab_path)
    return ag_path, ab_path

def getChainFromFile(pdb_path, chn_lst, out_path):
    struct = getPdbStruct(pdb_path)
    if type(chn_lst) != list:
        chn_lst = list(chn_lst)
    writePdbStruct2File(struct, out_path, ChainSelect(chn_lst))

def getAgEpitopeResidue(ag_path, ab_path, cutoff=4):
    ag_epitope_resi_lst = []
    ag_struct = getPdbStruct(ag_path)
    ab_struct = getPdbStruct(ab_path)
    agab_resi_product = itertools.product(ag_struct.get_residues(), ab_struct.get_residues())
    for ag_res, ab_res in agab_resi_product:
        atom_pair_product = itertools.product(ag_res.get_list(), ab_res.get_list())
        atom_pair_record = 0
        for atom1, atom2 in atom_pair_product:
            if atom_dist(atom1, atom2) <= cutoff:
                atom_pair_record = 1
                break
        if atom_pair_record == 1:
#            chain_id = ag_res.get_parent().get_id()
            res_id = str(ag_res.get_id()[1])
#            res_uni_id = chain_id + '_' + res_id
            res_uni_id = res_id
            ag_epitope_resi_lst.append(res_uni_id)
    return ag_epitope_resi_lst

def createAgEpitopeFile(ag_path, ag_epi_resi_lst, ag_epitope_path):
    #将一条chain筛选的resi从pdb文件中找出来并重新写到pdb文件
    struct = getPdbStruct(ag_path)
    writePdbStruct2File(struct, ag_epitope_path, ResiSelect(ag_epi_resi_lst))
    
def extractEpitopeFromComplex(pdb_path, ag_chn, epitope_path, ab_chn_lst=None, cutoff=4):
    ag_path, ab_path = creatAbAgPdb(pdb_path, ag_chn, ab_chn_lst=ab_chn_lst)
    ag_epi_resi_lst = getAgEpitopeResidue(ag_path, ab_path, cutoff=cutoff)
    createAgEpitopeFile(ag_path, ag_epi_resi_lst, epitope_path)
    os.remove(ag_path)
    os.remove(ab_path)
    
if __name__ == '__main__':
    pdb_path = r'F:\1A2Y.pdb'
    ag_chn = 'C'
    epitope_path = r'F:\1A2Y_epi.pdb'
    extractEpitopeFromComplex(pdb_path, ag_chn, epitope_path)
    print 'Done!'
