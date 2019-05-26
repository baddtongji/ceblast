import os
from Bio.PDB import PDBParser, PDBIO, Select

abbrev_mapping = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B", \
"CYS": "C", "GLN": "Q", "GLU": "E", "GLX": "Z", "GLY": "G", "HIS": "H", \
"ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", \
"SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "UNK": "X"}


def preprocessFile(pdbFile):
    f = open(pdbFile)
    txt = ''
    for line in f:
        if line.startswith('ATOM') or line.startswith('TER'):
            txt += line
    f.close()
    g = open(pdbFile, 'w')
    g.write(txt)
    return pdbFile

def getPdbStruct(pdbFile):
    p = PDBParser(PERMISSIVE=1)
    struct = p.get_structure(os.path.basename (pdbFile), pdbFile)
    return struct

def getPdbResiId(pdb_path):
    pdb_path = preprocessFile(pdb_path)
    struct = getPdbStruct(pdb_path)
    resi_id_list = []
    if len(struct[0].child_dict) == 1:#only one chain
        for resi in struct.get_residues():
            resi_id_list.append({'id':resi.get_id()[1],'name':resi.get_resname()})
    return resi_id_list

def writePdbStruct2File(struct, out_path, select_class=None):
    io = PDBIO()
    io.set_structure(struct)
    if select_class is None:
        io.save(out_path)
    else:
        io.save(out_path, select_class)

def pdb2seq(pdb_path):
    """pdb_path: it should better store only one chain. 
    And make sure that the residues' ids are continuous."""
    global abbrev_mapping
    
    struct = getPdbStruct(pdb_path)
    chains = [chn for chn in struct.get_chains()]
    chain_num = len(chains)
    
    if chain_num == 1:
        chain = chains[0]
        seq = ''
        for res in chain.get_residues():
            if res.resname in abbrev_mapping:
                abbr_aa = abbrev_mapping[res.resname]
            else:
                abbr_aa = 'X'
            seq += abbr_aa
        return seq
    elif chain_num > 1:
        print 'Moer than one chain, only the first one is selected!'
    else:
        print 'No chain found!'
        
def pdb2fasta(pdb_path, fas_path):
    """pdb_path: it should better store only one chain. 
    And make sure that the residues' ids are continuous."""
    global abbrev_mapping
    
    struct = getPdbStruct(pdb_path)
    Id = struct.id
    chains = [chn for chn in struct.get_chains()]
    chain_num = len(chains)
    res_list = []   #record residue's ID, name and aa_type
    
    if chain_num == 1:
        chain = chains[0]
        seq = ''
        for res in chain.get_residues():
            if res.resname in abbrev_mapping:
                abbr_aa = abbrev_mapping[res.resname]
            else:
                abbr_aa = 'X'
            res_list.append((res.id[1], abbr_aa))
            seq += abbr_aa
        g = open(fas_path, 'w')
        g.write('>'+Id+'\n')
        g.write(seq+'\n')
        g.close()
        return res_list
    elif chain_num > 1:
        print 'Moer than one chain, only the first one is selected!'
    else:
        print 'No chain found!'


def get_resi_name(pdb_path, id_list):
    global abbrev_mapping
    
    id_name_list = []
    struct = getPdbStruct(pdb_path)
    for res in struct.get_residues():
        if res.id[1] in id_list:
            id_name = str(res.id[1]) + abbrev_mapping[res.resname]
            id_name_list.append(id_name)
    return id_name_list

if __name__ == '__main__':
    pdb_path = r'F:\aaa.pdb'
    pdb_path = preprocessFile(pdb_path)
    struct = getPdbStruct(pdb_path)
    out_path = r'F:\bbb.pdb'
    
    class ResiSelect(Select):
        def __init__(self, Id_lst):
            self.Id_lst = Id_lst
        def accept_residue(self, residue):
            if residue.get_id()[1] in self.Id_lst:
                return 1
            else:
                return 0
                
    writePdbStruct2File(struct, out_path, ResiSelect([82, 83, 100]))
