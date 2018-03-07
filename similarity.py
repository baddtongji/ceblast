"""
define the `similarity` of two chains
"""

import re
from collections import OrderedDict,defaultdict
from UserDict import UserDict
import math
from numpy import corrcoef,array,vstack,zeros
import numpy as np

class residue_fp(object):
    """
    The fingerprint for residue
    """
    def __init__(self,fp_str,comp,residue_id):
        """
        Parameter:
        fp_str: the string representation of the fingerprint
        comp: the complex associated with the residue
        residue_id: the residue id
        """
        self.fp_str = fp_str
        self.complex = comp
        self.residue_id = residue_id

        self.res = filter (lambda r: r.get_id() == self.residue_id, self.complex.residues) [0] #search for the residue in the complex that has id = self.residue_id
        
    def get_edit_dist_to(self,residue_fp1):
        """
        Measure the disntace of self.fp_str to the other finger print using the correlation coefficient

        Paramter: 
        residue_fp1: the other finger print

        Return:
        The distance, float
        """
        return corrcoef(vstack((self.fp_str,residue_fp1.fp_str)))[0,1]

    def within_range_of(self, o_res, range_dist):
        """
        determine if the residue is within distance `range_dist` of the `o_res`

        Parameter:
        o_res: the object residue
        range_dist: the distance, int

        Return: 
        boolean
        """
        def distance(xyz1,xyz2):
            x1,y1,z1 = tuple(xyz1)
            x2,y2,z2 = tuple(xyz2)
            return math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

        outer_c = 0
        
        o_res_total = len(o_res.get_list())
        res_total = len(self.res.get_list())
        
        # o_res_total = len(o_res.atom)
        # res_total = len(self.res.atom)
        for a1 in self.res.get_list ():
            inner_c = 0
            for a2 in o_res.get_list ():
                if distance(a1.get_coord(), a2.get_coord()) < range_dist:
                    inner_c += 1
            if inner_c >= o_res_total / 2.:
                outer_c += 1
        if outer_c >= res_total / 2.:
            return True
        return False
    
    def get_string(self):
        """
        Get the fingerprint string representation
        """
        return self.fp_str
    
    def __repr__(self):
        return str(self.fp_str)
    __str__ = __repr__

class residue_fp_list(UserDict):
    """
    The fingerprints of a list of residues
    """
    def __init__(self, fp_string = '', comp = None):
        """
        Parameter:
        fp_string: the total/big string of the fingerprints, in the format: res_id1\t1,0,1,1\nres_id2\t1,1,0,0
        comp: the associating complex
        """
        UserDict.__init__(self)

        self.comp = comp
        self.data = OrderedDict()#the place to store the fingerprints
        
        from StringIO import StringIO
        for line in StringIO (fp_string).readlines():
            s_line = re.split(r"[ ,\t]",line) 
            residue_id,fp = int(s_line[0]), map(float, s_line[1:])
            
            self.data[residue_id] = residue_fp(fp, self.comp, residue_id)

class dist_mat(UserDict):
    """
    The distance matrix that measures pairwise distance between residue fingerprints from complex A and residue fingerprints from complex B 
    """
    
    def __init__(self,fp1,fp2,resi_match=0):
        from scipy import isnan
        """
        Parameter:
        fp1: the first fingerprint list 
        fp2: the second fingerprint list
        resi_match: boolean value (default 0), indicating whether the residues in fp1 and fp2 are matched/corresponding.
        """
        self.data = defaultdict(dict)

        self.fp1 = fp1.fp
        self.fp2 = fp2.fp
        
        if resi_match:
            resi_pair = [(fp1.fp.keys()[idx], fp2.fp.keys()[idx]) for idx in range(len(fp1.fp))]
            for res1, fp1 in self.fp1.items():
                for res2, fp2 in self.fp2.items():
                    if (res1, res2) in resi_pair:#filter self.data
                        self.data[res1][res2] = fp1.get_edit_dist_to(fp2)
                        #self.data[res2][res1] = self.data[res1][res2]
                        if isnan(fp1.get_edit_dist_to(fp2)):
                            print res1, fp1 
                            print res2, fp2
                            raw_input()
        else:
            for res1, fp1 in self.fp1.items():
                for res2, fp2 in self.fp2.items():
                    self.data[res1][res2] = fp1.get_edit_dist_to(fp2)
                    if isnan(fp1.get_edit_dist_to(fp2)):
                        print res1, fp1 
                        print res2, fp2
                        raw_input()

        self.clustered_fp1_res = set()
        self.clustered_fp2_res = set()
        self.clusters = []
        
    def find_closest_tuple(self):
        """
        find the closest residue pairs in terms of their fingerprint distance meanwhile ignoring those already in the clusters or considered not suitable
        
        Return: 
        the residue pair with the **maximum** fingerprint similarity with their similarity score
        """
        from itertools import chain
        
        all_pairs = list(chain.from_iterable(map(lambda (row, cols): map (lambda (col, val): (row, col, val), cols.items ()), self.data.items ())))
        #from defaultdict(<type 'dict'>, {res1: {res2: 3, res3: 3}, res2: {res4: 6}}) to [(res1, res2, 3), (res1, res3, 3), (res2, res4, 6)]
        
        tuples = filter(lambda tpl: 
                        tpl not in self.not_suitable_tuple and 
                        tpl [0] not in self.clustered_fp1_res and 
                        tpl [1] not in self.clustered_fp2_res,
                        all_pairs)  #filter out those in the not_suitable list
                
        if len (tuples) == 0:
            return None
        else:
            return max (tuples, key = lambda (_,__,num): num)


    def find_clusters(self, cutoff):
        def helper():
            """
            Perform one round of residue cluster discovery and add the discoverd cluster to the class variable
            
            steps: 
            1. find the closest pair as the cluster base
            2. use the base to find within-range pairs in a greedy-manner, always consider those pairs whose are more similar with each other
            """
            
            cluster = set()
            
            self.not_suitable_tuple = set()
            
            self.extending_tuple = self.find_closest_tuple() #we can the closest residue as the expansion point
            
            if self.extending_tuple is None: #if no more, it's done
                return

            cluster.add(self.extending_tuple); ext_res1 , ext_res2 , dist = self.extending_tuple #add the first tuple to the cluster and unpack it
            
            self.clustered_fp1_res.add(ext_res1); self.clustered_fp2_res.add(ext_res2)#add residues to corresponding sides
            
            center_res1, center_res2 = self.fp1[ext_res1].res, self.fp2[ext_res2].res #get the expansion residues
            
            while True:
                t = self.find_closest_tuple() #get the next closest(edit distance) tuple
                if t is None:#cannot find any appropriate tuple, quit
                    break
                
                res1 , res2 , dist = t
                if self.fp1[res1].within_range_of(center_res1, cutoff) and self.fp2[res2].within_range_of(center_res2, cutoff):
                    #check if it is within the range of extending tuple,if so, add it to the cluster
                    cluster.add(t)

                    self.clustered_fp1_res.add(res1)
                    self.clustered_fp2_res.add(res2)
                    #print "cluster size: %d,total residue number: f1 = %d, f2 = %d" %(len(cluster),len(self.fp1),len(self.fp2))
                else:
                    #print "out of range"
                    self.not_suitable_tuple.add(t)
                                
            self.clusters.append(cluster)
            self.extending_tuple = None
#            print 'cluster...', cluster
            return cluster

        while helper():
            pass
#        if self.fp1 != self.fp2:
#            print 'final cluster...', self.clusters
        return self.clusters
        
#res_sim_mat = {"AA" : 1,"AC" : 0,"AD" : 0,"AE" : 0,"AF" : 0,"AG" : 0,"AH" : 0,"AI" : 0,"AK" : 0,"AL" : 0,"AM" : 0,"AN" : -2,"AP" : -1,"AQ" : -1,"AR" : -1,"AS" : 1,"AT" : 0,"AV" : 0,"AW" : -3,"AY" : -2,"CC" : 9,"CD" : -3,"CE" : -4,"CF" : -2,"CG" : -3,"CH" : -3,"CI" : -1,"CK" : -3,"CL" : -1,"CM" : -1,"CN" : -3,"CP" : -3,"CQ" : -3,"CR" : -3,"CS" : -1,"CT" : -1,"CV" : -1,"CW" : -2,"CY" : -2,"DD" : 6,"DE" : 2,"DF" : -3,"DG" : -1,"DH" : -1,"DI" : -3,"DK" : -1,"DL" : -4,"DM" : -3,"DN" : 1,"DP" : -1,"DQ" : 0,"DR" : -2,"DS" : 0,"DT" : -1,"DV" : -3,"DW" : -4,"DY" : -3,"EE" : 5,"EF" : -3,"EG" : -2,"EH" : 0,"EI" : -3,"EK" : 1,"EL" : -3,"EM" : -2,"EN" : 0,"EP" : -1,"EQ" : 2,"ER" : 0,"ES" : 0,"ET" : -1,"EV" : -2,"EW" : -3,"EY" : -2,"FF" : 6,"FG" : -3,"FH" : -1,"FI" : 0,"FK" : -3,"FL" : 0,"FM" : 0,"FN" : -3,"FP" : -4,"FQ" : -3,"FR" : -3,"FS" : -2,"FT" : -2,"FV" : -1,"FW" : 1,"FY" : 3,"GG" : 6,"GH" : -2,"GI" : -4,"GK" : -2,"GL" : -4,"GM" : -3,"GN" : 0,"GP" : -2,"GQ" : -2,"GR" : -2,"GS" : 0,"GT" : -2,"GV" : -3,"GW" : -2,"GY" : -3,"HH" : 8,"HI" : -3,"HK" : -1,"HL" : -3,"HM" : -2,"HN" : 1,"HP" : -2,"HQ" : 0,"HR" : 0,"HS" : -1,"HT" : -2,"HV" : -3,"HW" : -2,"HY" : 2,"II" : 4,"IK" : -3,"IL" : 2,"IM" : 1,"IN" : -3,"IP" : -3,"IQ" : -3,"IR" : -3,"IS" : -2,"IT" : -1,"IV" : 3,"IW" : -3,"IY" : -1,"KK" : 5,"KL" : -2,"KM" : -1,"KN" : 0,"KP" : -1,"KQ" : 1,"KR" : 2,"KS" : 0,"KT" : -1,"KV" : -2,"KW" : -3,"KY" : -2,"LL" : 4,"LM" : 2,"LN" : -3,"LP" : -3,"LQ" : -2,"LR" : -2,"LS" : -2,"LT" : -1,"LV" : 1,"LW" : -2,"LY" : -1,"MM" : 5,"MN" : -2,"MP" : -2,"MQ" : 0,"MR" : -1,"MS" : -1,"MT" : -1,"MV" : 1,"MW" : -1,"MY" : -1,"NN" : 6,"NP" : -2,"NQ" : 0,"NR" : 0,"NS" : 1,"NT" : 0,"NV" : -3,"NW" : -4,"NY" : -2,"PP" : 7,"PQ" : -1,"PR" : -2,"PS" : -1,"PT" : -1,"PV" : -2,"PW" : -4,"PY" : -3,"QQ" : 5,"QR" : 1,"QS" : 0,"QT" : -1,"QV" : -2,"QW" : -2,"QY" : -1,"RR" : 5,"RS" : -1,"RT" : -1,"RV" : -3,"RW" : -3,"RY" : -2,"SS" : 4,"ST" : 1,"SV" : -2,"SW" : -3,"SY" : -2,"TT" : 5,"TV" : 0,"TW" : -2,"TY" : -2,"VV" : 4,"VW" : -3,"VY" : -1,"WW" : 11,"WY" : 2,"YY" : 7}
aa_20 = "FLIMVSPTAYHQNKDCWRGE"
res_sim_mat = {'GW': -2.0, 'GV': 0.0, 'GT': 1.0, 'GS': 0.0, 'GR': -2.0, 'GQ': -2.0, 'GP': -2.0, 'GY': -3.0, 'GG': 6.0, 'GF': -3.0, 'GE': -2.0, 'GD': -1.0, 'GC': -3.0, 'GA': 0.0, 'GN': -2.0, 'GM': -3.0, 'GL': -4.0, 'GK': -2.0, 'GI': -4.0, 'GH': -2.0, 'ME': -2.0, 'MD': -3.0, 'MG': -3.0, 'MF': 0.0, 'MA': -1.0, 'MC': -1.0, 'MM': 5.0, 'ML': 2.0, 'MN': -2.0, 'MI': 1.0, 'MH': -2.0, 'MK': -1.0, 'MT': -1.0, 'MW': -1.0, 'MV': -2.0, 'MQ': 0.0, 'MP': -2.0, 'MS': -1.0, 'MR': -1.0, 'MY': -1.0, 'FP': -4.0, 'FQ': -3.0, 'FR': -3.0, 'FS': -2.0, 'FT': -2.0, 'FV': -1.0, 'FW': 1.0, 'FY': 3.0, 'FA': -2.0, 'FC': -2.0, 'FD': -3.0, 'FE': -3.0, 'FF': 6.0, 'FG': -3.0, 'FH': -1.0, 'FI': 0.0, 'FK': -3.0, 'FL': 0.0, 'FM': 0.0, 'FN': -3.0, 'SY': -2.0, 'SS': 4.0, 'SR': -1.0, 'SQ': 0.0, 'SP': -1.0, 'SW': -3.0, 'SV': -2.0, 'ST': 1.0, 'SK': 0.0, 'SI': -2.0, 'SH': -1.0, 'SN': 1.0, 'SM': -1.0, 'SL': -2.0, 'SC': -1.0, 'SA': 1.0, 'SG': 0.0, 'SF': -2.0, 'SE': 0.0, 'SD': 0.0, 'YI': -1.0, 'YH': 2.0, 'YK': -2.0, 'YM': -1.0, 'YL': -1.0, 'YN': -2.0, 'YA': -2.0, 'YC': -2.0, 'YE': -2.0, 'YD': -3.0, 'YG': -3.0, 'YF': 3.0, 'YY': 7.0, 'YQ': -1.0, 'YP': -3.0, 'YS': -2.0, 'YR': -2.0, 'YT': -2.0, 'YW': 2.0, 'YV': -1.0, 'LF': 0.0, 'LG': -4.0, 'LD': -4.0, 'LE': -3.0, 'LC': -1.0, 'LA': -1.0, 'LN': -3.0, 'LL': 4.0, 'LM': 2.0, 'LK': -2.0, 'LH': -3.0, 'LI': 2.0, 'LV': 3.0, 'LW': -2.0, 'LT': -2.0, 'LR': -2.0, 'LS': -2.0, 'LP': -3.0, 'LQ': -2.0, 'LY': -1.0, 'RT': -1.0, 'RV': -3.0, 'RW': -3.0, 'RP': -2.0, 'RQ': 1.0, 'RR': 5.0, 'RS': -1.0, 'RY': -2.0, 'RD': -2.0, 'RE': 0.0, 'RF': -3.0, 'RG': -2.0, 'RA': -1.0, 'RC': -3.0, 'RL': -2.0, 'RM': -1.0, 'RN': 0.0, 'RH': 0.0, 'RI': -3.0, 'RK': 2.0, 'VH': -3.0, 'IP': -3.0, 'EM': -2.0, 'EL': -3.0, 'EN': 0.0, 'EI': -3.0, 'EH': 0.0, 'EK': 1.0, 'EE': 5.0, 'ED': 2.0, 'EG': -2.0, 'EF': -3.0, 'EA': -1.0, 'EC': -4.0, 'IT': -2.0, 'EY': -2.0, 'IW': -3.0, 'ET': 0.0, 'EW': -3.0, 'EV': -3.0, 'EQ': 2.0, 'EP': -1.0, 'ES': 0.0, 'ER': 0.0, 'II': 4.0, 'IH': -3.0, 'VR': -3.0, 'IM': 1.0, 'IN': -3.0, 'KC': -3.0, 'KA': -1.0, 'KG': -2.0, 'KF': -3.0, 'KE': 1.0, 'KD': -1.0, 'KK': 5.0, 'KI': -3.0, 'KH': -1.0, 'KN': 0.0, 'KM': -1.0, 'KL': -2.0, 'KS': 0.0, 'KR': 2.0, 'KQ': 1.0, 'KP': -1.0, 'KW': -3.0, 'KV': -3.0, 'KT': 0.0, 'KY': -2.0, 'DN': 1.0, 'DL': -4.0, 'DM': -3.0, 'DK': -1.0, 'DH': -1.0, 'DI': -3.0, 'DF': -3.0, 'DG': -1.0, 'DD': 6.0, 'DE': 2.0, 'DC': -3.0, 'DA': -2.0, 'DY': -3.0, 'DV': -3.0, 'DW': -4.0, 'DT': 1.0, 'DR': -2.0, 'DS': 0.0, 'DP': -1.0, 'DQ': 0.0, 'QQ': 5.0, 'QP': -1.0, 'QS': 0.0, 'QR': 1.0, 'QT': 0.0, 'QW': -2.0, 'QV': -2.0, 'QY': -1.0, 'QA': -1.0, 'QC': -3.0, 'QE': 2.0, 'QD': 0.0, 'QG': -2.0, 'QF': -3.0, 'QI': -3.0, 'QH': 0.0, 'QK': 1.0, 'QM': 0.0, 'QL': -2.0, 'QN': 0.0, 'WG': -2.0, 'WF': 1.0, 'WE': -3.0, 'WD': -4.0, 'WC': -2.0, 'WA': -3.0, 'WN': -4.0, 'WM': -1.0, 'WL': -2.0, 'WK': -3.0, 'WI': -3.0, 'WH': -2.0, 'WW': 11.0, 'WV': -3.0, 'WT': -3.0, 'WS': -3.0, 'WR': -3.0, 'WQ': -2.0, 'WP': -4.0, 'WY': 2.0, 'PR': -2.0, 'PS': -1.0, 'PP': 7.0, 'PQ': -1.0, 'PV': -2.0, 'PW': -4.0, 'PT': 1.0, 'PY': -3.0, 'PC': -3.0, 'PA': -1.0, 'PF': -4.0, 'PG': -2.0, 'PD': -1.0, 'PE': -1.0, 'PK': -1.0, 'PH': -2.0, 'PI': -3.0, 'PN': -1.0, 'PL': -3.0, 'PM': -2.0, 'CK': -3.0, 'CI': -1.0, 'CH': -3.0, 'CN': -3.0, 'CM': -1.0, 'CL': -1.0, 'CC': 9.0, 'CA': 0.0, 'CG': -3.0, 'CF': -2.0, 'CE': -4.0, 'CD': -3.0, 'CY': -2.0, 'CS': -1.0, 'CR': -3.0, 'CQ': -3.0, 'CP': -3.0, 'CW': -2.0, 'CV': -1.0, 'CT': -1.0, 'IY': -1.0, 'VA': 0.0, 'VC': -1.0, 'VD': -3.0, 'VE': -2.0, 'VF': -1.0, 'VG': -3.0, 'IQ': -3.0, 'VI': 3.0, 'IS': -2.0, 'IR': -3.0, 'VL': 1.0, 'VM': 1.0, 'VN': -3.0, 'IV': 1.0, 'VP': -2.0, 'VQ': -2.0, 'IK': -3.0, 'VS': -2.0, 'VT': -2.0, 'IL': 2.0, 'VV': 4.0, 'VW': -3.0, 'IA': -1.0, 'VY': -1.0, 'IC': -1.0, 'IE': -3.0, 'ID': -3.0, 'IG': -4.0, 'IF': 0.0, 'HY': 2.0, 'HR': 0.0, 'HS': -1.0, 'HP': -2.0, 'HQ': 0.0, 'HV': -2.0, 'HW': -2.0, 'HT': 0.0, 'HK': -1.0, 'HH': 8.0, 'HI': -3.0, 'HN': 1.0, 'HL': -3.0, 'HM': -2.0, 'HC': -3.0, 'HA': -2.0, 'HF': -1.0, 'HG': -2.0, 'HD': 1.0, 'HE': 0.0, 'NH': -1.0, 'NI': -3.0, 'NK': 0.0, 'NL': -3.0, 'NM': -2.0, 'NN': 6.0, 'NA': -2.0, 'NC': -3.0, 'ND': 1.0, 'NE': 0.0, 'NF': -3.0, 'NG': 0.0, 'NY': -2.0, 'NP': -2.0, 'NQ': 0.0, 'NR': 0.0, 'NS': 1.0, 'NT': 0.0, 'NV': -3.0, 'NW': -4.0, 'TY': -2.0, 'TV': -2.0, 'TW': -3.0, 'TT': 4.0, 'TR': -1.0, 'TS': 1.0, 'TP': 1.0, 'TQ': 0.0, 'TN': 0.0, 'TL': -2.0, 'TM': -1.0, 'TK': 0.0, 'TH': 0.0, 'TI': -2.0, 'TF': -2.0, 'TG': 1.0, 'TD': 1.0, 'TE': 0.0, 'TC': -1.0, 'TA': -1.0, 'AA': 4.0, 'AC': 0.0, 'AE': -1.0, 'AD': -2.0, 'AG': 0.0, 'AF': -2.0, 'AI': -1.0, 'AH': -2.0, 'AK': -1.0, 'AM': -1.0, 'AL': -1.0, 'AN': -1.0, 'AQ': -1.0, 'AP': -1.0, 'AS': 1.0, 'AR': -1.0, 'AT': -1.0, 'AW': -3.0, 'AV': -2.0, 'AY': -2.0}
#deal with missing key
res_sim_mat['VK'] = -1.16
#get B?, ?B value
for i in aa_20:
    res_sim_mat['B'+i] = []
    res_sim_mat[i+'B'] = []
    for j in ['D', 'N']:
        res_sim_mat['B'+i].append(res_sim_mat[j+i])
        res_sim_mat[i+'B'].append(res_sim_mat[i+j])
    res_sim_mat['B'+i] = round(np.mean(res_sim_mat['B'+i]), 2)
    res_sim_mat[i+'B'] = round(np.mean(res_sim_mat[i+'B']), 2)
#get BB value
res_sim_mat['BB'] = []
for i in aa_20:
    res_sim_mat['BB'].append(res_sim_mat['B'+i])
    res_sim_mat['BB'].append(res_sim_mat[i+'B'])
res_sim_mat['BB'] = round(np.mean(res_sim_mat['BB']), 2)
#get X?, ?X value
for i in aa_20:
    res_sim_mat['X'+i] = []
    res_sim_mat[i+'X'] = []
    for j in aa_20:
        res_sim_mat['X'+i].append(res_sim_mat[j+i])
        res_sim_mat[i+'X'].append(res_sim_mat[i+j])
    res_sim_mat['X'+i] = round(np.mean(res_sim_mat['X'+i]), 2)
    res_sim_mat[i+'X'] = round(np.mean(res_sim_mat[i+'X']), 2)

#get XB, BX value
res_sim_mat['XB'] = []
for i in aa_20:
    res_sim_mat['XB'].append(res_sim_mat[i+'B'])
res_sim_mat['XB'] = round(np.mean(res_sim_mat['XB']), 2)#get XB value

res_sim_mat['BX'] = []
for i in aa_20:
    res_sim_mat['BX'].append(res_sim_mat['B'+i])
res_sim_mat['BX'] = round(np.mean(res_sim_mat['BX']), 2)#get XB value

#get XX value
res_sim_mat['XX'] = round(np.mean(res_sim_mat.values()), 2)
 
class FPWithComplex(object):
    """
    Fingerprint with complex
    """
    def __init__(self, complex, fp_string):
        """
        Parameter:
        complex: the complex structure
        fp_string: the complex's residues' fingerprint string
        """
        self.complex = complex #reads the structure
        self.fp = residue_fp_list(fp_string, complex) #instantiates the fp list   
    
def similarity_between(c1,c2, cutoff = 20, resi_match=0):
    """
    Parameter: 
    c1: complex fingerprint 1
    c2: complex fingerprint 2

    Return:
    Three part similiarity scores
    """
    #print "generating distance matrix"
    pair = dist_mat(c1,c2, resi_match)
    
    #print "finding clusters"
    clusters = pair.find_clusters(cutoff)
    
    #val1, val2 and val3 starts
    val1,val2,val3 = 0, 0, 0

    #print "calculating value1"
    for c in clusters:
        for t in c:
            val1 += t[2] * 10#the edit distance

            #print "calculating value3"
    for c in clusters:
        val2 += len(c)#pair count

    #print "calculating value2"
    for c in clusters:
        for t in c:
            res1,res2,dist = t
            res1_code = pair.fp1[res1].res.get_resname_abbrev()
            res2_code = pair.fp2[res2].res.get_resname_abbrev()
            try:
                val3 += res_sim_mat[res1_code + res2_code]
            except KeyError:
                val3 += res_sim_mat[res2_code + res1_code]

    return val1 , val2 , val3
