from __future__ import division

import os

from similarity import FPWithComplex, similarity_between

def help ():
    print   """
    Usage:

    Without any specification: python get_sim_score.py --query-pdb test/data/sample1.pdb  --against-pdb test/data/sample2.pdb 
    
    With specification: python get_sim_score.py --query-pdb test/data/sample1.pdb  --query-epitope 211,213,214,224,225,226,227,228,229 --against-pdb test/data/sample2.pdb --against-epitope 216,217,218,219,220,221

    With spinimage configuration: python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb  --spin-image-radius-step=2 --spin-image-height-step=5 --sphere-radius-step=2

    You can specify the fingerprint path and skip the fp calcualtion step by using the option --query-fp and --against_fp, for example: 
    
    python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb --query-fp=fp/h2-v2/sample1.pdb.fp --against-fp=fp/h2-v2/sample2.pdb.fp
    """
def get_similarity_score(query_pdb, against_pdb, query_epitope=[], against_epitope=[],\
                        spin_image_radius_range = (0, 20), spin_image_height_range = (-30, 10), sphere_radius_range = (0, 20),\
                        spin_image_height_step = 5, spin_image_radius_step = 2, sphere_radius_step = 2, cutoff = 20.0, resi_match = 0, glycosylate = 0,\
                        query_fp_path=None, against_fp_path=None):
                            
    query_pdb_path = query_pdb
    against_pdb_path = against_pdb
    if query_epitope != []:
        query_epitope = map(int, query_epitope.split(','))
    if against_epitope != []:
        against_epitope = map(int, against_epitope.split(','))

    #calculate the finger print
    from get_fp import Complex
    from Bio.PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)

    if query_pdb_path is None:
        query_struct = None
    else:
        query_struct = p.get_structure(os.path.basename (query_pdb_path), query_pdb_path)
    
    if against_pdb_path is None:
        against_struct = None
    else:
        against_struct = p.get_structure(os.path.basename (against_pdb_path), against_pdb_path)

    query_complex = Complex (query_struct, query_epitope)
    against_complex = Complex (against_struct, against_epitope)
    
    if glycosylate:
        physchem_type_num = 4
    else:
        physchem_type_num = 3
    
    #if query_fp is not given
    if query_fp_path is None:
        query_complex.get_fp(spin_image_radius_step = spin_image_radius_step, spin_image_height_step = spin_image_height_step, sphere_radius_step = sphere_radius_step, physchem_type_num = physchem_type_num)        
        query_fp_string = query_complex.fp2str ()
    else:
        #if fp is given, read it from file
        with open (query_fp_path, 'r') as f1:
            query_fp_string = f1.read()
    
    #if against_fp is not given
    if against_fp_path is None:
        against_complex.get_fp(spin_image_radius_step = spin_image_radius_step, spin_image_height_step = spin_image_height_step, sphere_radius_step = sphere_radius_step, physchem_type_num = physchem_type_num)
        against_fp_string = against_complex.fp2str ()
    else:
        #if fp is given, read it from file
        with open(against_fp_path, 'r') as f2:
            against_fp_string = f2.read()
        
    query = FPWithComplex (query_complex, query_fp_string)
    against = FPWithComplex (against_complex, against_fp_string)
    
    score1, score2, score3 = similarity_between (query, against, cutoff = cutoff, resi_match = resi_match)
    z1, z2, z3 = similarity_between (query, query, cutoff = cutoff, resi_match = resi_match) #the normalization constant
    #print score1, score2, score3
    #print score1/z1, score2/z2, score3/z3
    return score1, score2, score3, sum ([score1, score2, score3]) / sum ([z1, z2, z3])

if __name__ == "__main__":
    import sys, getopt

    try:
        #parse the cmd argument
        optlist, args = getopt.getopt (sys.argv[1:], "", \
        ['query-pdb=', 'against-pdb=', 'query-fp=', 'against-fp=', \
        'query-epitope=', 'against-epitope=', \
        'spin-image-radius-step=', 'spin-image-height-step=', 'sphere-radius-step=', \
        'glycosylate=', 'resi_match=', 'cutoff='])
    except:
        help ()
        sys.exit (-1)

    #and make them into the right data type
    spin_image_radius_range = (0, 20)
    spin_image_height_range = (-30, 10)
    sphere_radius_range = (0, 20)
    
    spin_image_height_step = 5
    spin_image_radius_step = 2
    sphere_radius_step = 2
    
    glycosylate = 0
    resi_match = 0

    cutoff = 20.0
    
    query_epitope = []
    against_epitope = []

    query_fp_path, against_fp_path = None, None

    while len(optlist) > 0:
        opt, val = optlist.pop ()
        if opt == '--spin-image-radius-step':
            spin_image_radius_step = float (val)
        elif opt == '--spin-image-height-step':
            spin_image_height_step = float (val)
        elif opt == '--sphere-radius-step':
            sphere_radius_step = float (val)
        elif opt == '--cutoff':
            cutoff = float (val)
        elif opt == '--query-pdb':
            query_pdb_path = val
        elif opt == '--against-pdb':
            against_pdb_path = val
        elif opt == '--query-fp':
            query_fp_path = val
        elif opt == '--against-fp':
            against_fp_path = val
        elif opt == '--query-epitope':
            query_epitope = map (int, val.split (','))
        elif opt == '--against-epitope':
            against_epitope = map (int, val.split (','))
        #specific parameters for flu
        elif opt == '--glycosylate':
            glycosylate = int (val) #boolean
        elif opt == '--resi_match':
            resi_match = int (val) #boolean
        else:
            raise Exception ("Invalid option")

    #calculate the fingerprint
    from get_fp import Complex
    from Bio.PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)

    query_struct = p.get_structure(os.path.basename (query_pdb_path), query_pdb_path)
    against_struct = p.get_structure(os.path.basename (against_pdb_path), against_pdb_path)

    query_complex = Complex (query_struct, query_epitope)
    against_complex = Complex (against_struct, against_epitope)

    if glycosylate:
        physchem_type_num = 4
    else:
        physchem_type_num = 3
    
    #if query_fp is not given
    if query_fp_path is None:
        query_complex.get_fp(spin_image_radius_step = spin_image_radius_step, spin_image_height_step = spin_image_height_step, sphere_radius_step = sphere_radius_step, physchem_type_num = physchem_type_num)        
        query_fp_string = query_complex.fp2str ()
    else:
        #if fp is given, read it from file
        with open (query_fp_path, 'r') as f1:
            query_fp_string = f1.read()
    
    #if against_fp is not given
    if against_fp_path is None:
        against_complex.get_fp(spin_image_radius_step = spin_image_radius_step, spin_image_height_step = spin_image_height_step, sphere_radius_step = sphere_radius_step, physchem_type_num = physchem_type_num)
        against_fp_string = against_complex.fp2str ()
    else:
        #if fp is given, read it from file
        with open(against_fp_path, 'r') as f2:
            against_fp_string = f2.read()
    
    query = FPWithComplex (query_complex, query_fp_string)
    against = FPWithComplex (against_complex, against_fp_string)

    score1, score2, score3 = similarity_between (query, against, cutoff = cutoff, resi_match = resi_match)
    z1, z2, z3 = similarity_between (query, query, cutoff = cutoff, resi_match = resi_match) #the normalization constant
    
    #print score1, score2, score3
    #print score1/z1, score2/z2, score3/z3
    print score1, score2, score3, sum ([score1, score2, score3]) / sum ([z1, z2, z3])
