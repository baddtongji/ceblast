# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 09:29:13 2015

@author: Administrator
"""
"""每条抗原链与标准序列做比对，找出抗原上的表位位点"""
import os
from django.conf import settings
from Bio.Blast.Applications import _Ncbiblast2SeqCommandline
from Bio.Blast import NCBIXML
_static_path = settings.STATIC_DOC

def read_as_file(path):
    sites = []
    f = open(path)
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        site = cols[1]
        sites.append(int(site))
    return list(set(sites))

def blast_with_std_seq(fa_path, out_path, e_val_thr=0.9, subtype='h3n2'):
    os.chdir(r'/usr/share/ncbi-blast+/bin')
    if subtype == 'h3n2':
        std_path = _static_path + r'epiblast_data/std_seq/h3n2_std_seq.fas'
    elif subtype == 'h1n1':
        std_path = _static_path + r'epiblast_data/std_seq/h1n1_std_seq.fas'
    blast2seq_cline = _Ncbiblast2SeqCommandline(cmd='blastp', query=fa_path, subject=std_path, outfmt=5, out=out_path)
    stdout, stderr = blast2seq_cline()
    result_handle = open(out_path)
    blast_records = NCBIXML.parse(result_handle)
    max_score = -1000.
    result = {}
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <= e_val_thr:
                    if float(hsp.score) >= max_score:
                        max_score = float(hsp.score)
                        result['e_value'] = hsp.expect
                        result['score'] = hsp.score
#                        result['identities'] = hsp.identities
#                        result['positives'] = hsp.positives
                        result['gaps'] = hsp.gaps
                        result['q_start'] =  int(hsp.query_start)
                        result['s_start'] = int(hsp.sbjct_start)
#                        result['match'] = hsp.match
                        result['q_seq'] = hsp.query
                        result['s_seq'] = hsp.sbjct
    return result

def get_hsp_epi_res(hsp_res, query_resi_lst, site_path):
    """
    hsp_res: a hsp result (a dictionary) from 'blast_with_std_seq'.
    """
    if os.path.exists(site_path):
        sites = read_as_file(site_path)
    else:
        #user defined as_type has no site yet
        sites = [] 
    res = hsp_res
    #user-defined as_type: set query_first_resi_id and highlight_flags empty
    if sites == []:
        return [], [0]*len(res['q_seq'])
        
    query_first_resi_id = query_resi_lst[0][0]
    query_resi_ids = []
    highlight_flags = [0]*len(res['q_seq'])
    #if query & subject both have no gap
    if int(res['gaps']) == 0:
        for idx, aa in enumerate(res['q_seq']):
            if res['s_start'] + idx in sites:
                query_id = query_first_resi_id + res['q_start'] + idx - 1
                query_resi_ids.append(str(query_id))
                highlight_flags[idx] = 1
    else:
        #if query has gaps
        if '-' in res['q_seq']:
            q_gap_num = 0
            for idx, aa in enumerate(res['q_seq']):
                if aa == '-':
                    q_gap_num += 1
                    if res['s_start'] + idx in sites:
                        query_resi_ids.append('-')
                        highlight_flags[idx] = 1
                else:
                    if res['s_start'] + idx in sites:
                        query_id = query_first_resi_id + res['q_start'] + idx -1 - q_gap_num
                        query_resi_ids.append(str(query_id))
                        highlight_flags[idx] = 1
        #if subject has gaps
        elif '-' in res['s_seq']:
            s_gap_num = 0
            for idx, aa in enumerate(res['q_seq']):
                if res['s_seq'][idx] == '-':
                    s_gap_num += 1
                    continue
                if res['s_start'] + idx - s_gap_num in sites:
                    query_id = query_first_resi_id + res['q_start'] + idx - 1
                    query_resi_ids.append(str(query_id))
                    highlight_flags[idx] = 1
    return query_resi_ids, highlight_flags

