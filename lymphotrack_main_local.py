import argparse
import time
import pickle
import re
import numpy as np
import pandas as pd
import openpyxl
import os
import clone_operations
import process_sample
import Kevinsmerge
import process_sample
# import clone_operations
# import subclone
import sys

import subclone
import time

parser = argparse.ArgumentParser()
parser.add_argument('--sample_sheet', type=str, required=True)
parser.add_argument('--AnalyzePositiveControl', type=int, required=True)
parser.add_argument('--Match_Clone_with_FR', type=int, required=True)
parser.add_argument('--Input_Nanograms', type=int, required=True)
parser.add_argument('--skiprows', type=int, required=True)
parser.add_argument('--readrows', type=int, required=True)
parser.add_argument('--clone_length_fr1', type=int, required=True)
parser.add_argument('--clone_length_fr2', type=int, required=True)
parser.add_argument('--clone_length_fr3', type=int, required=True)
parser.add_argument('--independent_kmer', type=int, required=True)
parser.add_argument('--must_be_identical', type=float, required=True)
parser.add_argument('--aggregate', type=int, required=True)
parser.add_argument('--save_path',  type=str, required=True)
args = parser.parse_args()

print('argument parsing')
sample_sheet = pd.read_excel(args.sample_sheet)
if args.skiprows > 0:
    sample_sheet= sample_sheet.iloc[args.skiprows:]
if args.readrows > 0:    
    sample_sheet= sample_sheet.iloc[0:args.readrows]
    
sample_id = []
sample_replicate_number = []
sample_framework_number = []
clone_sequence_original_location = []
clone_length = []
clone_sequence = []
clone_match_to_total = []
clone_total = []
subclone_obs = {}
for index, row in sample_sheet.iterrows():
    print(index)
    print(row)
    for AnalyzePositiveControl in np.arange(0, 1+args.AnalyzePositiveControl):
        print('about to process sample sheet')
        
        all_subclone_tables = []
        fr_clones = []
        sample_path, clone_dictionary = process_sample.process_samplesheet(row,AnalyzePositiveControl)
        K_values = [args.clone_length_fr1, args.clone_length_fr2, args.clone_length_fr3]
        k_spliced_dictionary_fr1 = clone_operations.create_clone_dictionary_with_k_splice(clone_dictionary, args.clone_length_fr1)
        k_spliced_dictionary_fr2 = clone_operations.create_clone_dictionary_with_k_splice(clone_dictionary, args.clone_length_fr2)
        k_spliced_dictionary_fr3 = clone_operations.create_clone_dictionary_with_k_splice(clone_dictionary, args.clone_length_fr3)
        print('made dictionary')
        independent_k = args.independent_kmer
        framework_count = 0
        #if AnalyzePositiveControl:
        #    row['Sample ID']=row['Sample ID']+'PositiveControl'
        for framework_path in sample_path:
               
            framework=re.findall('FR[1-3]',framework_path)[0]    
            framework_count += 1
            
            if framework_count == 1:
                k_spliced_dictionary = k_spliced_dictionary_fr1
                if args.Match_Clone_with_FR:
                    print('match clone with fr')
                    keys=[mm for mm in k_spliced_dictionary_fr1.keys() if re.findall('FR1',mm) !=[]]
                else:
                    keys = k_spliced_dictionary_fr1.keys()   
            elif framework_count == 2:
                k_spliced_dictionary = k_spliced_dictionary_fr2
                if args.Match_Clone_with_FR:
                    keys=[mm for mm in k_spliced_dictionary_fr2.keys() if re.findall('FR2',mm) !=[]]
                else:
                    keys = k_spliced_dictionary_fr2.keys()             
            else:
                k_spliced_dictionary = k_spliced_dictionary_fr3
                if args.Match_Clone_with_FR:
                    keys=[mm for mm in k_spliced_dictionary_fr3.keys() if re.findall('FR3',mm) !=[]]
                else:
                    keys = k_spliced_dictionary_fr3.keys()        
            if framework_count >0:
            
                for key in keys:
                    print(key)
                    #breakpoint()
                    subclone_obj_a1 = subclone.Subclone(framework_path, framework_count, k_spliced_dictionary[key], independent_k, clone_dictionary[key], args.Input_Nanograms,args.must_be_identical, False, row, key, args.save_path,args.aggregate)
                    subclone_obs[row['Sample ID']+'_'+'Rep'+str(row['Replicate Number'])+'_'+framework+'.'+key] = subclone_obj_a1
                    
                    # Record information for final excel sheet
                    sample_id.append(row['Sample ID'])
                    sample_replicate_number.append(row['Replicate Number'])
                    sample_framework_number.append(framework_count)
                    clone_sequence_original_location.append(key)
                    clone_length.append(len(clone_dictionary[key]))
                    clone_sequence.append(clone_dictionary[key])
                    #clone_match_to_total.append(subclone_obj_a1.fr_total)
                    clone_total.append(sum(subclone_obj_a1.cluster['count']))
                    all_subclone_tables.append(subclone_obj_a1.fr_subclone_seqs)
                    fr_clones.append(clone_dictionary[key])
        
                

pickle.dump( subclone_obs, open( "subcloneob.pkl", "wb" ) )
if args.AnalyzePositiveControl != 0:
    process_sample.calculate_quantity(subclone_obs)                


