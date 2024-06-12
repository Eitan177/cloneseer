from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import re

def combine(positivecontrol,subclone_obj,positivecontrolcellequivalents=50):
   
    if len(positivecontrol.fr_subclone_seqs) ==0:
        positivecontrol.fr_subclone_seqs = pd.DataFrame(columns= ['sum of counts'])
    if len(positivecontrol.fr_subclone_seqs_unfiltered) ==0:
        positivecontrol.fr_subclone_seqs_unfiltered = pd.DataFrame(columns= ['sum of counts'])
    if len(positivecontrol.cluster) ==0:  
        positivecontrol.cluster = pd.DataFrame(columns= ['count'])  
    if len(subclone_obj.fr_subclone_seqs) ==0:
        subclone_obj.fr_subclone_seqs = pd.DataFrame(columns= ['sum of counts'])
    if len(subclone_obj.fr_subclone_seqs_unfiltered) ==0:
        subclone_obj.fr_subclone_seqs_unfiltered = pd.DataFrame(columns= ['sum of counts'])
    if len(subclone_obj.cluster) ==0:  
        subclone_obj.cluster = pd.DataFrame(columns= ['count'])  


    positivecontrol_readtotal=positivecontrol.fr_subclone_seqs['sum of counts'].sum()    
    positivecontrol_readtotal_unfiltered=positivecontrol.fr_subclone_seqs_unfiltered['sum of counts'].sum()
    read_total_bcell=subclone_obj.cluster['count'].sum()
    

    subclone_obj.fr_subclone_seqs['cell equivalents'] = positivecontrolcellequivalents*(subclone_obj.fr_subclone_seqs['sum of counts']/positivecontrol_readtotal)
    subclone_obj.fr_subclone_seqs_unfiltered['cell equivalents'] = positivecontrolcellequivalents*(subclone_obj.fr_subclone_seqs_unfiltered['sum of counts']/positivecontrol_readtotal)
    subclone_obj.fr_subclone_seqs['fraction of bcells']=subclone_obj.fr_subclone_seqs['sum of counts']/read_total_bcell
    subclone_obj.fr_subclone_seqs_unfiltered['fraction of bcells']=subclone_obj.fr_subclone_seqs_unfiltered['sum of counts']/read_total_bcell
    positivecontrol.fr_subclone_seqs['fraction of bcells']= positivecontrol.fr_subclone_seqs['sum of counts']/read_total_bcell
    positivecontrol.fr_subclone_seqs_unfiltered['fraction of bcells']= positivecontrol.fr_subclone_seqs_unfiltered['sum of counts']/read_total_bcell
    cell_equivalents_all_dna = (subclone_obj.input_nanograms*1000/6.4)+positivecontrolcellequivalents
    read_total_subclone=subclone_obj.fr_subclone_seqs['sum of counts'].sum()
    read_total_subclone_unfiltered=subclone_obj.fr_subclone_seqs_unfiltered['sum of counts'].sum()
    
    fraction_subclone=read_total_subclone/read_total_bcell
    fraction_subclone_unfiltered=read_total_subclone_unfiltered/read_total_bcell
    
    total_cell_equivalents=subclone_obj.fr_subclone_seqs['cell equivalents'].sum()
    total_cell_equivalents_unfiltered=subclone_obj.fr_subclone_seqs_unfiltered['cell equivalents'].sum()
    fraction_subclone_in_all_dna=total_cell_equivalents/cell_equivalents_all_dna
    fraction_subclone_in_all_dna_unfiltered=total_cell_equivalents_unfiltered/cell_equivalents_all_dna

    stats_combined = pd.DataFrame({'Read Total':[read_total_bcell],'Subclone Read Total':[read_total_subclone],'Positive Control Read Total':positivecontrol_readtotal, 'Cell Equivalents':[total_cell_equivalents],'Fraction of B Cells':[fraction_subclone],'Positive Control Cell Equivalents':[positivecontrolcellequivalents],
                                   'Cell Equivalents All DNA':[cell_equivalents_all_dna],'Fraction of Clone In All DNA':[fraction_subclone_in_all_dna]})
    stats_combined_unfiltered = pd.DataFrame({'Read Total':[read_total_bcell],'Subclone Read Total Unfiltered':[read_total_subclone_unfiltered],'Positive Control Read Total Unfiltered':positivecontrol_readtotal_unfiltered, 'Cell Equivalents Unfiltered':[total_cell_equivalents_unfiltered],'Fraction of B Cells Unfiltered':[fraction_subclone_unfiltered],'Positive Control Cell Equivalents':[positivecontrolcellequivalents],
                                      'Cell Equivalents All DNA':[cell_equivalents_all_dna],'Fraction of Clone In All DNA Unfiltered':[fraction_subclone_in_all_dna_unfiltered]})

    combine_excel_name=subclone_obj.save_path + "/" + subclone_obj.sample_id + "_Rep" + str(subclone_obj.sample_replicate_number) + "_FR" + str(subclone_obj.fr) + "_" + subclone_obj.clone_id + "_Subclonetable_w_PositiveControl.xlsx"
    
    with pd.ExcelWriter(combine_excel_name) as writer: 
        subclone_obj.fr_subclone_seqs.to_excel(writer, sheet_name='subclone_information_vdj_filtered')
        positivecontrol.fr_subclone_seqs.to_excel(writer, sheet_name='positivecontrol_information_vdj_filtered')
        stats_combined.to_excel(writer, sheet_name='stats_combined_vdj_filtered')
        subclone_obj.fr_subclone_seqs_unfiltered.to_excel(writer, sheet_name='subclone_information_not_filtered')
        positivecontrol.fr_subclone_seqs_unfiltered.to_excel(writer, sheet_name='positivecontrol_information_not_filtered')
        stats_combined_unfiltered.to_excel(writer, sheet_name='stats_combined_not_filtered')        

def calculate_quantity(subclone_obs):
    
    sorted_subclone_names=pd.DataFrame({'subclones':subclone_obs.keys()}).sort_values(by='subclones')
    sorted_subclone_names['subclones_remove_endkey']=sorted_subclone_names['subclones'].apply(lambda x: re.sub('[.].+$','',x))
    
    for jj in np.where(sorted_subclone_names.duplicated(['subclones_remove_endkey']))[0]:
        combine_needed = sorted_subclone_names['subclones_remove_endkey'].iloc[jj]
        objs_for_combine = sorted_subclone_names['subclones'][sorted_subclone_names['subclones_remove_endkey']==combine_needed]
        positive_control_keyname=objs_for_combine[[re.findall('PositiveControl',y)!=[] for y in objs_for_combine]].iloc[0]
        sample_keyname=objs_for_combine[[re.findall('PositiveControl',y)==[] for y in objs_for_combine]].iloc[0]
        
        combine(subclone_obs[positive_control_keyname],subclone_obs[sample_keyname])
        




def find_sample_id_path(row):
    """
    Gets the full sample ID path for each framework based upon the sample ID path and the sample ID.

    Row contains Sample ID Path, Sample ID, Replicate Number, Clone File Path, Clone File, and Descriptors
    Returns the full sample path for each framework of the Sample ID. Otherwise Raise Exception
    """
    
    all_files_fr1 = [f for f in listdir(row['Sample ID Path'] + "/IGH_FR1_output") if isfile(join(row['Sample ID Path'] + "/IGH_FR1_output", f))]
    for item in all_files_fr1:
        
        if re.match(str(row['Sample ID']) + ".+"+"_L001_001_combined.fastq_unique_reads.tsv", item):
            #breakpoint()
            sampleID = item#str(row['Sample ID']) + "_L001_001_combined.fastq_unique_reads.tsv"
            return row['Sample ID Path'] + "/IGH_FR1_output/" + sampleID, row['Sample ID Path'] + "/IGH_FR2_output/" + sampleID, row['Sample ID Path'] + "/IGH_FR3_output/" + sampleID
        elif re.match(str(row['Sample ID']) + ".+"+"_combined.fastq_unique_reads.tsv", item): 
            sampleID = item
            return row['Sample ID Path'] + "/IGH_FR1_output/" + sampleID, row['Sample ID Path'] + "/IGH_FR2_output/" + sampleID, row['Sample ID Path'] + "/IGH_FR3_output/" + sampleID
    raise Exception("Sample ID not present in the Sample Path")


def all_clones(clone_path):
    """
    Gets all the clones in the clone path. Opens the file and then creates a dictionary with keys as SampleID_Framework# and values as the Clone Sequence
    """
    #breakpoint()
    try:
        clone_file = pd.read_csv(clone_path)
    except: #UnicodeDecodeError:
        pass
    try:
        clone_file = pd.read_excel(clone_path)
    except: #UnicodeDecodeError:
        pass #raise Exception(" You're file is not a csv, tsv, or excel file. Unable to process")
    try:
        clone_file = pd.read_excel(re.sub('[.]xlsx$','',clone_path))
    except: #UnicodeDecodeError:
        raise Exception(" You're file is not a csv, tsv, or excel file. Unable to process")

    clone_matches = {}
    for index, row in clone_file.iterrows():
        clone_matches[row['Sample ID'] + "_FR" + str(row['Framework']) + "_Clone" + str(row['Clone Number'])] = row['Clone Sequence']
    return clone_matches


def process_samplesheet(row,PositiveControl_analyze=0):
    """
    Process the samplesheet amd return out the path for each one of the frameworks along with a clone dictionary that lists each clone id with its sequence.
    """
    sample_path_fr1, sample_path_fr2, sample_path_fr3 = find_sample_id_path(row)
    clone_path = row['Clone File Path'] + "/" + row['Clone File'] + ".xlsx"
    # Hard coded Lymphotrack positive control sequence, FR1
    if PositiveControl_analyze:
        clone_dictionary = {'PositiveControl':'CTTCTGGAGGCACCTTCAGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATTACCGCGGACGAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGATAGGCGCGGGGAATGGCCTCCCTCGGATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCAC'}
    else:
        clone_dictionary = all_clones(clone_path)
    return [sample_path_fr1, sample_path_fr2, sample_path_fr3], clone_dictionary
