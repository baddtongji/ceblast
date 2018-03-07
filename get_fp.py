#-*-coding:gbk-*-
"""
Finger print generation(110 bit version)
"""

import os, sys, glob, warnings, math
import numpy as np

from collections import defaultdict

from vec_geom import angle, length

from dist import coord_dist, res_ca_dist

from get_glycosylation_site import get_glyco_site

class Residue(object):
    """
    The residue class used for fingerprint generation
    """
    abbrev_mapping = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B", "CYS": "C", "GLN": "Q", "GLU": "E", "GLX": "Z", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "UNK": "X"}
    #下面是三种物化性质，可见是针对20中氨基酸残基打分的
    hydro_dict = {'A':0.61,'C':1.07,'D':0.46,'E':0.47,'F':2.02,'G':0.07,'H':0.61,'I':2.22,'K':1.15,'L':1.53,'M':1.18,'N':0.06,'P':1.95,'Q':0.0,'R':0.6,'S':0.05,'T':0.05,'V':1.32,'W':2.65,'Y':1.88, 'X': 1.00, 'B':0.26}
    charged_dict={'A':-0.01,'C':0.12,'D':0.15,'E':0.07,'F':0.03,'G':0.0,'H':0.08,'I':-0.01,'K':0.0,'L':-0.01,'M':0.04,'N':0.06,'P':0.0,'Q':0.05,'R':0.04,'S':0.11,'T':0.04,'V':0.01,'W':0.0,'Y':0.03, 'X': 0.04, 'B':0.105}
    h_bond_dict={'A':0,'C':0,'D':1,'E':1,'F':0,'G':0,'H':1,'I':0,'K':2,'L':0,'M':0,'N':2,'P':0,'Q':2,'R':4,'S':1,'T':1,'V':0,'W':1,'Y':1, 'X': 0.85, 'B':1.5}

    def __init__(self, res, comp):
        
        self.c = comp
        self.fp = None 
        self.ca = None
        self.body = res
        self.resnum = res.get_id ()
        
        for a in res: #iterate over the atoms in the residue
            if a.get_name () == 'CA':
                self.ca = a
                break
        if self.ca is None:
            warnings.warn ("resnum %s has no CA atom" %str(res.get_id ()))
            
    def get_id (self):
        """ residue id"""
        return self.body.get_id () [1] #get residue id number e.g.(' ', 129, ' ')[1] = '129'

    def get_list (self):
        """
        return the atom list
        """
        return self.body.get_list () #get atom list

    def get_resname_abbrev (self):
        """
        get the abbreviation code of residue name
        """
        return Residue.abbrev_mapping[self.body.get_resname()]

    def set_fp_length (self, length):
        """
        set the length of the fingerprint
        """
        self.fp = [0] * length
        
    def is_valid (self):
        """
        if the residue is valid (contain CA atom or not)
        """
        return self.ca is not None
        
    def turn_on_bit(self,bit_num,count):
        """
        for the bit position, `bit_num`, set the value to `count`
        """
        self.fp[bit_num] = count#第几位的fingerprint小格子里有多少residue

    def get_surrounding_res(self):
        """
        get surrouding residues
        """
        for res in self.c.get_residues ():
            #if res.ca == self.ca:continue
            yield res
    
    def get_next_2_res(self):
        """
        get next two residues for determining if the site is glycosylated
        """
        all_res_list = self.c.get_all_residues()
        res_num = self.get_id()
        ind = all_res_list.index(self)
        #check the index
        if ind + 2 < len(all_res_list):
            next_res1 = all_res_list[ind+1]
            next_res2 = all_res_list[ind+2]
            #check if the next two residues exist
            if next_res1.get_id() == res_num + 1 and next_res2.get_id() == res_num + 2:
                return next_res1, next_res2
    
    def glyco_site(self):
        next_2_res = self.get_next_2_res()
        if next_2_res is not None:
            res1, res2 = next_2_res
            res_resname_abbr = self.get_resname_abbrev()
            next_res1_resname_abbr = res1.get_resname_abbrev()
            next_res2_resname_abbr = res2.get_resname_abbrev()
            return get_glyco_site(res_resname_abbr, next_res1_resname_abbr, next_res2_resname_abbr)
        #无论是非糖基化位点还是末尾的氨基酸，都认为不是糖基化位点
        return 0
        
    def get_struct_fp(self, spin_image_radius_step, spin_image_height_step, 
                      spin_image_radius_ind_min, spin_image_radius_ind_max, 
                      spin_image_height_ind_min, spin_image_height_ind_max,
                      spin_image_radius_seg_cnt, spin_image_height_seg_cnt,
                      offset = 0):
        """
        get the structural based fingerprints, given a set of spinimage parameters and the fingerprint bit starting position#就是offset了吧
        """
        def get_region_number(res_ca):#计算residue是属于哪个小格子里的
            """
            according to the given residue CA atom, calculate which spin image region this residue should belong to
            """
            my_point = np.array(self.ca.get_coord ())#self.ca是当前residue的alpha碳原子
            v1 = my_point - np.array([0,0,0])#the original point
            v2 = np.array(res_ca.get_coord ()) - my_point#res_ca是参数给的alpha碳原子
            
            ang = angle(v1, v2)

            v2_len = length(v2)
            
            x, y = v2_len * math.sin(ang), v2_len * math.cos(ang)#投影！！！
            
            if x < 0:
                raise ValueError("x cannot be negative")

            x_ind  = math.floor( x / spin_image_radius_step)#返回浮点型的整数部分e.g.24.0
            y_ind  = math.floor( y / spin_image_height_step)
            
            if x_ind >= spin_image_radius_ind_min and x_ind <= spin_image_radius_ind_max and\
               y_ind >= spin_image_height_ind_min and y_ind <= spin_image_height_ind_max:
                #within the range
                return int((y_ind - spin_image_height_ind_min) * spin_image_radius_seg_cnt + x_ind)
            else:
                return None
            
        #type one fp
#        print 'min_max', spin_image_radius_ind_min, spin_image_radius_ind_max, spin_image_height_ind_min, spin_image_height_ind_max
        d_ = defaultdict(int)#告诉你它的默认value一定是int类型的
        for res in self.get_surrounding_res():
            num = get_region_number(res.ca)
            if num is not None:
                d_[num] += 1
        
        for bit_num,count in d_.iteritems():#将d_储存到bit里
            self.turn_on_bit(offset + bit_num, count)#1~80

        return self.fp
    
    def get_surrounding_fp(self, dist_step, dist_ind_min, dist_ind_max, include_glyco = 0, offset = 0):#算剩下的30(or 40)bit
        """
        get the surrouding-residue based fingerprints
        include_glyco: 决定是否输出糖基化描述符
        """
#        print self, dist_ind_min, dist_ind_max, dist_step
        #type two fp
        d_ = defaultdict(list)
        for other in self.get_surrounding_res():
            dist = res_ca_dist(other,self)#算两个alpha C的距离
            dist_ind = math.floor( dist / dist_step )
    
            if dist_ind >= dist_ind_min and dist_ind < dist_ind_max + 1:#说明这个res在目标res邻近范围内
                #within range
                d_[dist_ind].append(other)
                
        slice_count = dist_ind_max - dist_ind_min #how many slices for the sphere;=9-0=9
        
        for i in xrange(slice_count + 1):
            # for layer i
            h_bond, charged , hydro = 0 , 0 , 0
            for res in d_[i]:
                # for res in layer i
                code = res.get_resname_abbrev()
                hydro += Residue.hydro_dict[code]#求和
                charged += Residue.charged_dict[code]
                h_bond += Residue.h_bond_dict[code]
            
            #fp for layer i,in the 3 aspects
            self.turn_on_bit(offset + i , hydro)#hydro是x+0~x+9
            self.turn_on_bit(offset + (slice_count + 1) + i , charged)#x+10~x+19
            self.turn_on_bit(offset + (slice_count + 1)*2 + i , h_bond)#x+20~x+29
        
        if include_glyco:
            for i in xrange(slice_count + 1):
                # for layer i
                glyco = 0
                for res in d_[i]:
                    # for res in layer i
                    if res.glyco_site():
                        glyco += 1
                        
                self.turn_on_bit(offset + (slice_count + 1)*3 + i , glyco)#x+30~x+39

    def __repr__(self):
        #return "ca atom index:%d" % (self.ca.index)
        return "%s" % self.get_id()

class Complex(object):
    def __init__(self, pdb_fp, epitope = []):
        """
        pdb_fp: the pdb structure
        epitope: the epitope residue number list
        """
        self.st = pdb_fp

        all_epitope = (len(epitope) == 0)

        self.residues = []
        self.all_residues = []
            
        for res in self.st.get_residues():
            res = Residue(res, self)
            self.all_residues.append(res)
            # if res.is_valid () and (all_epitope or (len(epitope) != 0 and res.resnum in epitope)):
            if res.is_valid () and (all_epitope or (len(epitope) != 0 and res.get_id () in epitope)): #if the res is valid and filter those residue in the epitope
                self.residues.append(res)

    def get_residues (self):
        return self.residues
    
    def get_all_residues(self):
        return self.all_residues
    
    def cont_residues(self):
        if epitope == []:
            return 0
        return 1
        
    def get_fp(self, 
               spin_image_radius_range = (0, 20),
               spin_image_height_range =  (-30, 10),
               sphere_radius_range = (0, 20),
               spin_image_radius_step = 2,
               spin_image_height_step = 5,
               sphere_radius_step = 2,
               physchem_type_num = 3):
        """
        get the 110-bit fingerprint
        physchem_type_num: 用户提供的理化性质个数，默认为3
        """
        #radius part
        spin_image_radius_min , spin_image_radius_max = spin_image_radius_range
        spin_image_radius_seg_cnt = ( spin_image_radius_max - spin_image_radius_min ) / spin_image_radius_step
        spin_image_radius_ind_min , spin_image_radius_ind_max = int(spin_image_radius_min / spin_image_radius_step), int(spin_image_radius_max / spin_image_radius_step - 1) #index min and max must be integer. Maybe some warning should be put here.

        #height part
        spin_image_height_min , spin_image_height_max = spin_image_height_range
        spin_image_height_seg_cnt = ( spin_image_height_max - spin_image_height_min ) / spin_image_height_step
        spin_image_height_ind_min , spin_image_height_ind_max = int(spin_image_height_min / spin_image_height_step), int(spin_image_height_max / spin_image_height_step - 1) #index min and max must be integer. Maybe some warning should be put here.

        cylinder_slice_count = (spin_image_height_ind_max - spin_image_height_ind_min + 1) * (spin_image_radius_ind_max - spin_image_radius_ind_min + 1)
        
        #sphere part
        dist_min, dist_max = sphere_radius_range
        dist_step = sphere_radius_step
        
        dist_ind_min, dist_ind_max =  int(dist_min / dist_step), int(dist_max / dist_step - 1)#index min and max must be integer. Maybe some warning should be put here.
        
        sphere_slice_count = dist_ind_max - dist_ind_min + 1
        
        fp_size = cylinder_slice_count + sphere_slice_count * physchem_type_num
        self.res_list = []
            
        for i, res in enumerate(self.residues):
            
            res.set_fp_length (fp_size)#initialize the fingerprint
            
            #print "residue %d" %i
            if physchem_type_num == 3:
                include_glyco = 0
            elif physchem_type_num == 4:
                include_glyco = 1
            
            res.get_surrounding_fp(dist_step, dist_ind_min, dist_ind_max, include_glyco = include_glyco, offset = cylinder_slice_count)
            
            res.get_struct_fp(spin_image_radius_step, spin_image_height_step,
                              spin_image_radius_ind_min, spin_image_radius_ind_max,
                              spin_image_height_ind_min, spin_image_height_ind_max,
                              spin_image_radius_seg_cnt, spin_image_height_seg_cnt)
            
#            print res.fp,len(res.fp)
            self.res_list.append(res)
        
        return self.res_list
    
    def fp2str (self):
        """
        print the fingerprint string representation
        """
        return "\n".join("%s\t%s" %(res.get_id()," ".join(map(str,res.fp))) for res in self.res_list)

    def write_fp_to_file(self,path):
        """
        write the fp string to file given path
        """
        if not os.path.exists( os.path.dirname(path) ):
            os.mkdir(os.path.dirname(path))
        with open(path,'w') as f:
            f.write(self.fp2str())
    
    def fp2str_with_given_res(self, given_res_list):
        """
        print the fingerprint string representation if it belongs to given res list[int list]
        """
#        print 'all_fp_len:', len(self.res_list)
#        print 'part_fp_len:', len([res for res in self.res_list if res.get_id() in given_res_list])
        return "\n".join("%s\t%s" %(res.get_id()," ".join(map(str,res.fp))) for res in self.res_list if res.get_id() in given_res_list)

    def write_fp_to_file_with_given_res(self, path, given_res_list):
        """
        write the fp string to file given path
        """
        if not os.path.exists( os.path.dirname(path) ):
            os.mkdir(os.path.dirname(path))
        with open(path,'w') as f:
            f.write(self.fp2str_with_given_res(given_res_list))
        
if __name__ == "__main__":
    import sys
    
    from Bio.PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)
    
    try:
        path = 'F:\\20141013_agab\\epi_para_tope_whole\\ab_ag_l3\\1A2Y\\1A2Y_C.pdb' #sys.argv[1]
        horizontal = 5  #int(sys.argv[2])
        vertical = 4    #int(sys.argv[3])
        
        epitope = [13,18,19,21,22,23,24,25,26,27,102,103,116,117,118,119,120,121,122,124,125,126,127,128,129]#[] if len(sys.argv) == 4 else map(int, sys.argv[-1].split (','))
    except:
        print "Expecting format like: \n\tpython get_fp.py path/to/pdb/file horizontal-step vertical-step [a list of epitope residue numbers]"
        print 'For example: \n\tpython get_fp.py test/data/sample1.pdb  5 4 211,213,214,224,225,226,227,228,229'
        sys.exit (-1)
        
    struct = p.get_structure(os.path.basename (path), path)
    c = Complex(struct, epitope)
#    print type(c.get_residues()[2].get_id())
    c.get_fp (spin_image_radius_step = horizontal, spin_image_height_step = vertical, physchem_type_num = 3)
    
    print c.fp2str ()
    print str(c.res_list[0])
    print type(str(c.res_list[0]))
