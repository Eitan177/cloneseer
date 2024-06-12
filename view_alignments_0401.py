import time
from random import random
import os
import numpy as np
from bokeh.plotting import figure, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering
# from scipy.cluster.hierarchy import linkage
from fastcluster import linkage
from scipy.spatial.distance import pdist, squareform
import scipy
import streamlit as st
#import streamlit_ext as ste
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import pdb
from KevinsFunction import *
from view_alignments_0326 import *
from multiple_sequence_alignment import create_multiple_sequence_alignment, run_igblast


def make_fasta_file_from_reference_seq(reference_seq, reference_file):
    fasta_file = open(reference_file, "w")
    fasta_lines = []
    fasta_lines.append(">reference_seq" + "\n")
    fasta_lines.append(reference_seq + "\n")
    fasta_file.writelines(fasta_lines)
    fasta_file.close()


def plot_dendrogram(uniq, model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram

    dendrogram(linkage_matrix, **kwargs)


@st.experimental_memo
def getlinkage(uniqdf, distancemetric):  # , optimal_ordering):
    return (linkage(uniqdf, method=distancemetric))  # ,optimal_ordering=optimal_ordering))


@st.cache()
def convert_df(fr_df):
    return fr_df.to_csv(index=False).encode('utf-8')


def view_alignment(aln, use_consensus, fr, sample_id, sample_replicate_number, clone_id, save_path, reads_in_file,input_file_gui=True):
    """Bokeh sequence alignment view"""
    start = time.time()
    consensus_sequence = use_consensus  # Mother clone sequence
    reference_file = save_path + "/" + sample_id + "_reference_seq.fasta"
    make_fasta_file_from_reference_seq(consensus_sequence,reference_file)
    igblast_parsed_output_for_reference=run_igblast(reference_file)
    forfileseqs = [rec.annotations['seq'] for rec in aln]  # Sequences with inserts
    forfileinsert = [rec.annotations['insert'] for rec in aln]  # Insert positions
    seqs = [rec.seq for rec in aln]  # Sequences without the inserts
    counts = [int(np.ceil(int(rec.description))) for rec in aln]  # The count
    alnscores = [int(np.float(rec.name)) for rec in aln]  # The alignment score
    alnscores = np.asarray(alnscores, dtype='int64')
    colors = get_colors_consensus(seqs, consensus_sequence)

    x_axis = np.arange(1, len(seqs[0]) + 1)
    y_axis = np.arange(0, len(seqs), 1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x_axis, y_axis)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    gg = np.array(colors).reshape((-1, np.max(gx)))
    ugg = gg
    countgg = np.array(counts)

    k = np.unique(gg)

    k = k[k != 'white'][::-1]

    v = np.arange(0, len(k)) + 1

    all_nnn = np.empty(shape=gg.shape)

    unique_nnn = np.empty(shape=k.shape)

    unique_nnn = np.empty(shape=ugg.shape)

    for key, val in zip(k, v):
        unique_nnn[ugg == key] = val

    unique_nnn[ugg == 'white'] = 0

    # col_score_kms = st.columns(11)

    # setting distance_threshold=0 ensures we compute the full tree.
    dendo, axis = plt.subplots(1, 1, figsize=(3, 7))

    uniq, counts = unique_nnn, countgg
    if uniq.shape[0] > 1:
        uniqdf = pd.DataFrame(uniq)

        uniqdf = np.array(uniqdf)
        countsforlink = counts.copy()
        if np.argmax(alnscores) == np.argmax(countsforlink) and np.max(countsforlink) > 250:
            countsforlink[np.argmax(countsforlink)] = 1

        dfwcountnomain = np.repeat(uniqdf, countsforlink, axis=0)
        mismatchmutatedspotsuse = np.sum(dfwcountnomain != 0, axis=0) / dfwcountnomain.shape[0] > 0.01
        indelmutatedspotsuse = np.sum(dfwcountnomain == 5, axis=0) / dfwcountnomain.shape[0] > 0.001
        if len(np.where(~indelmutatedspotsuse)) <= 0:
            indexfirstletternotindel = np.where(~indelmutatedspotsuse)[0][0]
            if indexfirstletternotindel > 0:
                indelmutatedspotsuse[0:indexfirstletternotindel] = False
        mutatedspotsuse = mismatchmutatedspotsuse + indelmutatedspotsuse
        """
        with st.spinner('performing distance for ' + str(uniqdf.shape[0]) + ' sequences with ' + str(
                uniqdf.shape[1]) + ' positions and using' + str(np.sum(mutatedspotsuse)) + ' of those positions'):
        """
        Z2 = getlinkage(uniqdf, 'ward')
        # st.write('done with distance computation')

        countsforlink = countsforlink[scipy.cluster.hierarchy.leaves_list(Z2)]
        counts = counts[scipy.cluster.hierarchy.leaves_list(Z2)]

        seqs = np.array(seqs)[scipy.cluster.hierarchy.leaves_list(Z2)]
        forfileseqs = np.array(forfileseqs)[scipy.cluster.hierarchy.leaves_list(Z2)]
        forfileinsert = np.array(forfileinsert)[scipy.cluster.hierarchy.leaves_list(Z2)]
        strc = [str(hhh) if hhh > 10 else '' for hhh in counts]
        ugg = ugg[scipy.cluster.hierarchy.leaves_list(Z2)]
        alnscores = alnscores[scipy.cluster.hierarchy.leaves_list(Z2)]
        
    if len(forfileinsert) == 1:
        table_of_sequences = pd.DataFrame(
            {'sequence': seqs, 'sequenceNotformatted': forfileseqs, 'numberobserved': counts, 'inserts': forfileinsert,
             'alnscores': alnscores, 'mother_clone': consensus_sequence})

    else:
        table_of_sequences = pd.DataFrame(
            {'sequence': seqs, 'sequenceNotformatted': forfileseqs, 'numberobserved': counts,
             'inserts': forfileinsert.tolist(), 'alnscores': alnscores, 'mother_clone': consensus_sequence})
    ## 04/24 we do not need streamlit interactions for the linux version
    ## ste.download_button(label='proper order Download', data=convert_df(table_of_sequences), file_name='unique_data.csv',
    ##                    mime='text/csv')  # ,key=random())
    with st.expander('alignment view'):

        # 03/31 Kevins clustering and feature addition
        
        table_of_sequences_with_clusters_and_features = KevinsFunction(table_of_sequences)
         
        if not input_file_gui:
            table_of_sequences_with_clusters_and_features['closest match'] = np.arange(0, table_of_sequences_with_clusters_and_features.shape[0])
            table_of_sequences_with_clusters_and_features['outliers'] = 1
            table_of_sequences_with_clusters_and_features.sort_values(['alnscores'], inplace=True)

            unfiltered_sequences_ordered=make_alignmentplotwithcluster(table_of_sequences_with_clusters_and_features, consensus_sequence, sample_id, sample_replicate_number, clone_id, save_path=save_path, fr=fr, reorder=False, input_file_gui=input_file_gui,full_plot=True)
            make_alignmentplotwithcluster(table_of_sequences_with_clusters_and_features, consensus_sequence, sample_id, sample_replicate_number, clone_id, save_path=save_path, fr=fr, reorder=False, input_file_gui=input_file_gui,full_plot=False)
        table_of_sequences_with_clusters_and_features_without_unassigned = table_of_sequences_with_clusters_and_features[table_of_sequences_with_clusters_and_features['closest match'] != -1]

        table_of_sequences_with_clusters_and_features_just_unassigned = table_of_sequences_with_clusters_and_features[table_of_sequences_with_clusters_and_features['closest match'] == -1]
        
        table_of_outliers = table_of_sequences_with_clusters_and_features_without_unassigned.groupby(
            ['closest match']).apply(lambda x: np.array((x['sequence'][x['outliers'] > 0].iloc[0],
                                                         sum(x['numberobserved']),
                                                         x['index_mismatch'][x['outliers'] > 0].iloc[0],
                                                         x['homology_against_mother_clone'][x['outliers'] > 0].iloc[0],
                                                         x['inserts'][x['outliers'] > 0].iloc[0],
                                                         x['deletion_in_trinucleotide'][x['outliers'] > 0].iloc[0],
                                                         x['insertion_in_trinucleotide'][x['outliers'] > 0].iloc[0],
                                                         x['sequenceNotformatted'][x['outliers'] > 0].iloc[0],
                                                         len(x['sequence'][x['outliers'] > 0].iloc[0].replace("-", "")),
                                                         len(x['sequenceNotformatted'][x['outliers'] > 0].iloc[0].replace("-", "")))))#,
                                                         #x['v-genes'][x['outliers'] > 0].iloc[0], 
                                                         #x['v-genes score'][x['outliers'] > 0].iloc[0],
                                                         #x['v-genes percent'][x['outliers'] > 0].iloc[0],
                                                         #x['v-genes endmatch'][x['outliers'] > 0].iloc[0],
                                                         #x['v-genes fullmatch query'][x['outliers'] > 0].iloc[0],
                                                         #x['v-genes fullmatch ig element'][x['outliers'] > 0].iloc[0]        )))

        table_w_sequence_and_count = pd.DataFrame({'seq': np.vstack(table_of_outliers).T[0],
                                                   'sum of counts': [np.int(ii) for ii in np.vstack(table_of_outliers).T[1]],
                                                   'discrepant_positions': [str([np.int(jj) for jj in ii]) for ii in np.vstack(table_of_outliers).T[2]],
                                                   'mother_clone_distance': np.vstack(table_of_outliers).T[3],
                                                   'insertpositions': np.vstack(table_of_outliers).T[4],
                                                   'delete_trinucleotide': np.vstack(table_of_outliers).T[5],
                                                   'insert_trinucleotide': np.vstack(table_of_outliers).T[6],
                                                   'notformatted': np.vstack(table_of_outliers).T[7],
                                                   'length': np.vstack(table_of_outliers).T[8],
                                                   'unformattedlength': np.vstack(table_of_outliers).T[9],
                                                   'closest match':table_of_outliers.index})#,
                                                   #'V-gene': np.vstack(table_of_outliers).T[8],
                                                   #'V-gene score': np.vstack(table_of_outliers).T[9],
                                                   #'V-gene percent': np.vstack(table_of_outliers).T[10],
                                                   #'V-gene endmatch': np.vstack(table_of_outliers).T[11],
                                                   #'V-gene fullmatch query': np.vstack(table_of_outliers).T[12],
                                                   #'V-gene fullmatch ig element': np.vstack(table_of_outliers).T[13],})
        if len(table_of_sequences_with_clusters_and_features_just_unassigned) > 0:
            table_w_sequence_and_count_unassigned = pd.DataFrame(
                {'align length sequence ': table_of_sequences_with_clusters_and_features_just_unassigned['sequence'],
                 'sum of counts': table_of_sequences_with_clusters_and_features_just_unassigned['numberobserved'],
                 'discrepant_positions': table_of_sequences_with_clusters_and_features_just_unassigned['index_mismatch'].apply(lambda x: str(x)),
                 'mother_clone_distance': table_of_sequences_with_clusters_and_features_just_unassigned['homology_against_mother_clone'],
                 'insertpositions': table_of_sequences_with_clusters_and_features_just_unassigned['inserts'],
                    'delete_trinucleotide': table_of_sequences_with_clusters_and_features_just_unassigned['deletion_in_trinucleotide'],
                    'insert_trinucleotide': table_of_sequences_with_clusters_and_features_just_unassigned['insertion_in_trinucleotide'],
                 'notformatted': table_of_sequences_with_clusters_and_features_just_unassigned['sequenceNotformatted'],
                 'length': table_of_sequences_with_clusters_and_features_just_unassigned['sequence'].apply(lambda x: len(x.replace("-", ""))),
                 'unformattedlength': table_of_sequences_with_clusters_and_features_just_unassigned['sequenceNotformatted'].apply(lambda x: len(x.replace("-", ""))),
                 'closest match': np.arange(-1,table_of_sequences_with_clusters_and_features_just_unassigned['sequence']-1)* -1})
                 #'V-gene': table_of_sequences_with_clusters_and_features_just_unassigned['v-genes'],
                 #'V-gene score': table_of_sequences_with_clusters_and_features_just_unassigned['v-genes score'],
                 #'V-gene percent': table_of_sequences_with_clusters_and_features_just_unassigned['v-genes percent'],
                 #'V-gene endmatch': table_of_sequences_with_clusters_and_features_just_unassigned['v-genes endmatch'],
                 #'V-gene fullmatch query': table_of_sequences_with_clusters_and_features_just_unassigned['v-genes fullmatch query'],
                 #'V-gene fullmatch ig element': table_of_sequences_with_clusters_and_features_just_unassigned['v-genes fullmatch ig element'],})
            table_w_sequence_and_count = pd.concat([table_w_sequence_and_count, table_w_sequence_and_count_unassigned])

        # make_alignmentplotwithcluster(table_of_sequences_with_clusters_and_features, consensus_sequence, reorder=False)

        if input_file_gui:
            unfiltered_sequences_ordered=make_alignmentplotwithcluster(table_of_sequences_with_clusters_and_features, consensus_sequence, sample_id, sample_replicate_number, clone_id, save_path, fr=fr, reorder=True,full_plot=True)
            make_alignmentplotwithcluster(table_of_sequences_with_clusters_and_features, consensus_sequence, sample_id, sample_replicate_number, clone_id, save_path, fr=fr, reorder=True,full_plot=False)

        """
        # st.write('subclones')
        # st.write(table_w_sequence_and_count)
        ste.download_button(label='dendrogram ordered sequences',
                            data=convert_df(table_of_sequences_with_clusters_and_features),
                            file_name='dendrogramordered.csv', mime='text/csv')
        ste.download_button(label='subclones', data=convert_df(table_w_sequence_and_count),
                            file_name='subclonetable_' + str(fr) + '.csv', mime='text/csv')
        """

        if fr is None:
            table_of_sequences_with_clusters_and_features.to_excel(save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_" + clone_id + "_OrderedDendrogram.xlsx")
            table_w_sequence_and_count.to_excel(save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_" + clone_id + "_Subclonetable.xlsx")
            table_sequence_count_sorted_name = save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_" + clone_id + "_Subclonetable_Sorted.xlsx"
            table_w_sequence_and_count_sorted = table_w_sequence_and_count.sort_values(by=['seq'])
            #dfv, dfj, dfd = create_multiple_sequence_alignment(table_w_sequence_and_count_sorted, save_path, sample_id + "_Rep" + str(sample_replicate_number) + "_" + clone_id)
            table_w_sequence_count_sorted_w_vdj= create_multiple_sequence_alignment(table_w_sequence_and_count_sorted, save_path, sample_id + "_Rep" + str(sample_replicate_number) + "_" + clone_id)
        else:
            table_of_sequences_with_clusters_and_features.to_excel(save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + clone_id + "_OrderedDendrogram.xlsx")
            table_w_sequence_and_count.to_excel(save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + clone_id + "_Subclonetable.xlsx")
            table_w_sequence_and_count_sorted = table_w_sequence_and_count.sort_values(by=['seq'])
            table_sequence_count_sorted_name = save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + clone_id + "_Subclonetable_Sorted.xlsx"
            #dfv, dfj, dfd = create_multiple_sequence_alignment(table_w_sequence_and_count_sorted, save_path, sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + clone_id)
            
            table_w_sequence_count_sorted_w_vdj = create_multiple_sequence_alignment(table_w_sequence_and_count_sorted, save_path, sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + clone_id)
        
        unfiltered_table_w_sequence_count_sorted_w_vdj = pd.merge(unfiltered_sequences_ordered,table_w_sequence_count_sorted_w_vdj)
        unfiltered_table_w_sequence_count_sorted_w_vdj['read_total']=reads_in_file
        unfiltered_table_w_sequence_count_sorted_w_vdj.to_excel(table_sequence_count_sorted_name, index=False)
        igblast_parsed_output_for_reference['v_gene_family']=[re.sub('[*].+$','',y) for y in igblast_parsed_output_for_reference['top_v_gene_match']]
        igblast_parsed_output_for_reference['v_gene_locus']=[re.sub('[-].+$','',y) for y in igblast_parsed_output_for_reference['top_v_gene_match']]
        
        table_w_sequence_count_sorted_w_vdj['v_gene_locus']=[re.sub('[-].+$','',str(y)) for y in table_w_sequence_count_sorted_w_vdj['top_v_gene_match']]
        
        #if fr == 3:
        
        table_w_sequence_count_sorted_w_filtered_vdj=pd.merge(table_w_sequence_count_sorted_w_vdj,igblast_parsed_output_for_reference[['top_v_gene_match','top_j_gene_match','v_gene_family','v_gene_locus','v_end','j_start']],on=['top_j_gene_match','v_gene_locus'],how='left',indicator=True)
        
        table_w_sequence_count_sorted_w_filtered_vdj = table_w_sequence_count_sorted_w_filtered_vdj[table_w_sequence_count_sorted_w_filtered_vdj[['mother_clone_distance','_merge']].apply(lambda x: x[0]<3 or x[1]=='both',axis=1)].drop(columns=["_merge"])
            
        if fr == 3:
            
            table_w_sequence_count_sorted_w_filtered_vdj= table_w_sequence_count_sorted_w_filtered_vdj[table_w_sequence_count_sorted_w_filtered_vdj[['mother_clone_distance','insert_trinucleotide','delete_trinucleotide']].apply(lambda x: x[0]<10 and ((x[1] % 3 == 0 and x[1]<=15) and (x[2] % 3 == 0 and x[2]<=15)),axis=1)]
        else:
            table_w_sequence_count_sorted_w_filtered_vdj= table_w_sequence_count_sorted_w_filtered_vdj[table_w_sequence_count_sorted_w_filtered_vdj[['top_v_gene_match_x','top_v_gene_match_y','v_end_x','v_end_y','j_start_x','j_start_y','insert_trinucleotide','delete_trinucleotide','mother_clone_distance']].apply(lambda x: (re.sub(';.+','',str(x[0]))==re.sub(';.+','',str(x[1])) and str(x[2])==str(x[3]) and str(x[4])==str(x[5])) or (str(x[0])==str(x[1]) and ((x[6] % 3 == 0 and x[6]>2) or (x[7] % 3 == 0 and x[7]>2)) or (x[8]<3)),axis=1)]
        
        
        filtered_table_of_sequences_with_clusters_and_features=pd.merge(table_of_sequences_with_clusters_and_features,table_w_sequence_count_sorted_w_filtered_vdj['closest match'].drop_duplicates())
        clone_id = clone_id + '_filtered'
        
        if filtered_table_of_sequences_with_clusters_and_features.shape[0] > 0:
            
            sequences_ordered=make_alignmentplotwithcluster(filtered_table_of_sequences_with_clusters_and_features, consensus_sequence, sample_id, sample_replicate_number,clone_id, save_path=save_path, fr=fr, reorder=True, input_file_gui=input_file_gui,full_plot=True)
            make_alignmentplotwithcluster(filtered_table_of_sequences_with_clusters_and_features, consensus_sequence, sample_id, sample_replicate_number,clone_id, save_path=save_path, fr=fr, reorder=True, input_file_gui=input_file_gui,full_plot=False)
        else:
            sequences_ordered=filtered_table_of_sequences_with_clusters_and_features
        table_w_sequence_count_sorted_w_filtered_vdj['query id']=['Seq_'+ str(x) for x in np.arange(0,table_w_sequence_count_sorted_w_filtered_vdj.shape[0])]
            
        table_w_sequence_count_sorted_w_filtered_vdj=pd.merge(sequences_ordered,table_w_sequence_count_sorted_w_filtered_vdj)
        table_w_sequence_count_sorted_w_filtered_vdj['read_total']=reads_in_file
        
        table_w_sequence_count_sorted_w_filtered_vdj.to_excel(re.sub('[.]xlsx','filtered.xlsx',table_sequence_count_sorted_name),index=False)
        
     
    return [sum(counts), table_w_sequence_count_sorted_w_vdj, table_w_sequence_count_sorted_w_filtered_vdj] #table_w_sequence_and_count]


def get_colors_consensus(seqs, consensus_sequence):
    """make colors for bases in sequence"""
    # [i for i in range(len(s1)) if s1[i] != s2[i]]
    text = [[s[i], i] for s in list(seqs) for i in range(len(s))]

    def clrs(letter, pos):
        if consensus_sequence[pos] == letter and letter != '-':
            color = 'white'
        elif consensus_sequence[pos] == letter and letter == '-':
            color = 'gray'
        elif (letter != '-' and consensus_sequence[pos] == '-') or (letter == '-' and consensus_sequence[pos] != '-'):
            color = 'black'
        else:
            if letter == 'A':
                color = 'red'
            elif letter == 'T':
                color = 'blue'  # 'lightcoral'
            elif letter == 'G':
                color = 'brown'  # 'crimson'
            elif letter == 'C':
                color = 'turquoise'  # firebrick'
            elif letter == 'N':
                color = '#FFEFDB'
            else:
                color = 'lightblue'
        return color

    colors = [clrs(s, i) for s, i in text]
    return colors


def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', '-': 'white', 'N': 'brown'}
    colors = [clrs[i] for i in text]
    return colors
