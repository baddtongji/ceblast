# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 12:51:14 2015

@author: Administrator
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 09:09:41 2015

@author: Administrator
"""

#average 5 pdb structures
import os
from Bio.PDB import PDBParser, PDBIO, Residue
import numpy as np

def get_atom_info(a):
    name = a.get_name()
    coord = a.get_coord()
    bfac = a.get_bfactor()
    occ = a.get_occupancy()
    altloc = a.get_altloc()
    fullname = a.get_fullname()
    ser_num = a.get_serial_number()
    ele = a.element
    return name, coord, bfac, occ, altloc, fullname, ser_num, ele

def id_plus_one(struct, level="atom"):
    """
    struct: atom and residue is already marked as '-1'
    """
    #一旦发现-1这个标志，就将后面的id依次+1
    if level == 'atom':
        change = 0
        atom_list = [atom for atom in struct.get_atoms()]
        for idx, atom in enumerate(atom_list):
            if atom.serial_number == '-1':
                change = 1
                pre_atom = atom_list[idx-1]
                atom.set_serial_number(str(int(pre_atom.serial_number) + 1))
                continue
            if change:
                atom.set_serial_number(str(int(atom.serial_number) + 1))
    elif level == 'residue':
        change = 0
        chn = struct[0].child_list[0]
        for idx, res in enumerate(chn.child_list):
            if res.id[1] == -1:
                change = 1
                pre_res = chn.child_list[idx-1]
                res.id = (res.id[0], pre_res.id[1]+1, res.id[2])
                continue
            if change:
                res.id = (res.id[0], res.id[1]+1, res.id[2])
        #update child_dict
        chn.child_dict = {}
        for res in chn.child_list:
            chn.child_dict[res.id] = res
    return struct
                
def del_gap(pre_resi_id, struct):
    """
    pre_resi_id: an integer
    """
    chn = [chn for chn in struct.get_chains()][0]
    pre_res = None
    post_res = None
    #get pre resi and post resi and their CA atom
    for res in struct.get_residues():
        if res.id[1] == pre_resi_id:
            pre_res = res
        if res.id[1] == pre_resi_id + 1:
            post_res = res
    if pre_res is None or post_res is None:
        return
    else:
        ave_coord = []
        ave_bfac = []
        pre_res_ca = [atom for atom in pre_res.get_list() if atom.id == 'CA'][0]
        post_res_ca = [atom for atom in post_res.get_list() if atom.id == 'CA'][0]
        ave_coord.append(get_atom_info(pre_res_ca)[1])
        ave_coord.append(get_atom_info(post_res_ca)[1])
        ave_bfac.append(get_atom_info(pre_res_ca)[2])
        ave_bfac.append(get_atom_info(post_res_ca)[2])
        ave_coord = np.mean(ave_coord, axis=0)
        ave_bfac = np.mean(ave_bfac)
        #get new CA atom coord & bfac & No.
        new_res_ca = pre_res_ca.copy()
        new_res_ca.set_coord(ave_coord)
        new_res_ca.set_bfactor(ave_bfac)
        new_res_ca.set_serial_number('-1')#mark
        #get paranet residue name & id
        new_id = (pre_res.id[0], -1, pre_res.id[2])#mark
        new_res = Residue.Residue(new_id, 'UNK', pre_res.segid)
        new_res.add(new_res_ca)
        pos = chn.child_list.index(post_res)
        chn.insert(pos, new_res)
        #change all the atoms and residues id
        struct = id_plus_one(struct, level="atom")
        struct = id_plus_one(struct, level="residue")
        return struct

def get_no_gap_pdb(std_res_list, res_list, pdb_path, out_path):
    std_res_list = map(int, std_res_list)
    new_res_list = []
    for res in res_list:
        if res.isdigit():
            new_res_list.append(int(res))
        else:
            new_res_list.append(res)
    res_list = new_res_list
    zip_list = zip(std_res_list, res_list)
    zip_list.sort()
    std_res_list = [str(i) for i, j in zip_list]
    res_list = [str(j) for i, j in zip_list]
    
    p = PDBParser(PERMISSIVE=1)
    struct = p.get_structure(os.path.basename(pdb_path), pdb_path)
    
    while 'g' in res_list:
        for idx, res_id in enumerate(res_list):
            std_id = std_res_list[idx]
            pre_std_id = std_res_list[idx-1]
            pre_res_id = res_list[idx-1]
            id_diff = int(pre_std_id) - int(pre_res_id)
            if res_id == 'g':
                new_resi_id = int(std_id) - id_diff
                pre_resi_id = new_resi_id - 1
                struct = del_gap(pre_resi_id, struct)
                break
        #update res_list
        update = 0
        if_del_gap = 1
        for idx, Id in enumerate(res_list):
            if if_del_gap and Id == 'g':
                update = 1
                if_del_gap = 0
                res_list[idx] = str(new_resi_id)
                continue
            if update and Id != 'g':
                res_list[idx] = str(int(res_list[idx]) + 1)
#        print std_res_list, res_list
                  
    #write struct into file
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path)
#    print 'done!'
    return res_list
                 
if __name__ == '__main__':
    std_res_list = []
    res_list = []
    resi_path = r'/try/AAA16811.txt'
    f = open(resi_path)
    for line in f:
        r1, r2 = line.strip().split('\t')
        std_res_list.append(r1)
        res_list.append(r2)
    print std_res_list
    print res_list
    pdb_path = r'/try/AAA16811.pdb'
    out_path = r'/try/AAA16811_no_gap.pdb'
    get_no_gap_pdb(std_res_list, res_list, pdb_path, out_path)
    
