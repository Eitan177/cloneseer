import numpy as np
import pandas as pd
from Bio import pairwise2
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import subclone
import streamlit as st
#import streamlit_ext as ste

from view_alignments_0401 import convert_df


def framework_merge(framework_1, framework_2):
    """
    Merges two files between the same framework. Each index in framework number 1 is getting assigned correct merges from framework number 2
    :param framework_1: Arbitrary framework 1
    :param framework_2: Arbitrary framework 2
    :return: correct_merge: A list containing lists of matched merge sequences from framework 2 to each individual framework 1 sequence.
    """

    correct_merge = []
    framework_1_outliers = framework_1['Outliers'].tolist()
    framework_2_outliers = framework_2['Outliers'].tolist()

    for outlier1_position in range(len(framework_1_outliers)):
        position_merge = []
        for outlier2_position in range(len(framework_2_outliers)):
            shortest_length = min(len(framework_1_outliers[outlier1_position]),
                                  len(framework_2_outliers[outlier2_position]))
            alignments = pairwise2.align.globalms(framework_1_outliers[outlier1_position].replace("-", ""),
                                                  framework_2_outliers[outlier2_position].replace("-", ""), 1, -1, -1,
                                                  0, score_only=True)
            if shortest_length - alignments <= 2:
                position_merge.append(outlier2_position)
        correct_merge.append(position_merge)
    return correct_merge


def merge_outliers(outliers):
    """
    Appends 2 columns comparing against the 2 other framework for matches. Each column contains a list of indexes from the outlier file
    :param outliers: The dataframe of the 3 frameworks.
    :return: The dataframe of each framework with appended column comparing against index of other framework for matches.
             Column contains list of indexes from other framework that match.
    """

    merged_outliers = []
    for framework_number_1 in range(3):
        merges = []
        merge_column_name = []
        for framework_number_2 in range(3):
            if framework_number_1 != framework_number_2:
                merge_column_name.append(
                    'framework_' + str(framework_number_1) + '_merge_with_framework_' + str(framework_number_2))
                merges.append(framework_merge(outliers[framework_number_1], outliers[framework_number_2]))
        outliers[framework_number_1][merge_column_name[0]] = merges[0]
        outliers[framework_number_1][merge_column_name[1]] = merges[1]
        merged_outliers.append(outliers[framework_number_1])
    return merged_outliers


def heatmap_matrix_builder(outliers1, outliers2=None):
    """
    Builds the matrix for the heatmap matrix. Input either takes in a set of 3 frameworks containing outliers and does an
    interallele merging and creates a heatmap or input takes 2 sets of 3 frameworks of outliers and does a replicate
    comparison on them.
    :param outliers1: Outliers of 3 frameworks combined.
    :param outliers2: Optional. Outliers from 3
    :return:
    """

    if outliers2 is None:
        outliers2 = outliers1
    outlier_length_1 = len(outliers1)
    outlier_length_2 = len(outliers2)
    matrix = np.zeros((outlier_length_1, outlier_length_2))
    for outlier1_position in range(outlier_length_1):
        for outlier2_position in range(outlier_length_2):
            shortest_length = min(len(outliers1[outlier1_position].replace("-", "")),
                                  len(outliers2[outlier2_position].replace("-", "")))
            alignments = pairwise2.align.globalms(outliers1[outlier1_position].replace("-", ""),
                                                  outliers2[outlier2_position].replace("-", ""), 1, -1, -1, 0,
                                                  score_only=True)
            matrix[outlier1_position][outlier2_position] = shortest_length - alignments
    return matrix


def build_heatmap_graph(heatmap_matrix, rep_count_1, rep_count_2):
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(heatmap_matrix, linewidth=0.5)
    plt.title("Pairwise alignment score difference between two sequences")
    if rep_count_1 is None:
        plt.xlabel("Replicate Sequence Indexes")
        plt.ylabel("Replicate Sequence Indexes")
    else:
        plt.xlabel("Replicate " + str(rep_count_1 + 1) + " Sequence Indexes")
        plt.ylabel("Replicate " + str(rep_count_2 + 1) + " Sequence Indexes")

    st.write(fig)
    sns.set()
    matrix_merge_threshold = np.where(heatmap_matrix >= 6, 1, 0)
    fig1, ax1 = plt.subplots()
    sns.heatmap(matrix_merge_threshold, linewidth=0.5)
    plt.title("Merges between two sequences")
    if rep_count_1 is None:
        plt.xlabel("Replicate Sequence Indexes")
        plt.ylabel("Replicate Sequence Indexes")
    else:
        plt.xlabel("Replicate " + str(rep_count_1 + 1) + " Sequence Indexes")
        plt.ylabel("Replicate " + str(rep_count_2 + 1) + " Sequence Indexes")
    st.write(fig1)
    sns.set()


def get_sorted_index(outliers_count):
    """
    Sorts sequences by outlier_count. This needs to be fixed...
    """
    list_index = sorted(range(len(outliers_count)), key=lambda k: outliers_count[k])
    return reversed(list_index)


def Kevinsmerge(outliers_inter_allele_1, outliers_inter_allele_2, rep_count_1, rep_count_2):
    matrix = heatmap_matrix_builder(outliers_inter_allele_1, outliers_inter_allele_2)
    build_heatmap_graph(matrix, rep_count_1, rep_count_2)


def generate_significant_tables(outliers_inter_allele_1, outliers_count, mother_clone_dist, v_genes, matrix):
    significant_sequences = []
    significant_sequences_length = []
    significant_sequences_count = []
    significant_sequences_mother_clone_dist = []
    significant_v_genes = []

    list_of_possible_mother_clone_dist = sorted(set(mother_clone_dist))
    all_merged_indexes = []
    for unique_mother_clone_dist in list_of_possible_mother_clone_dist:
        for index in range(len(outliers_inter_allele_1)):
            if index not in all_merged_indexes and mother_clone_dist[index] == unique_mother_clone_dist:
                significant_sequences.append(outliers_inter_allele_1[index])
                significant_sequences_length.append(len(outliers_inter_allele_1[index]))
                significant_sequences_mother_clone_dist.append(mother_clone_dist[index])
                significant_v_genes.append(v_genes[index])

                sum_of_count = 0
                row_of_merged_matrix = np.where(matrix[index] >= 6, 1, 0)
                for position in range(len(row_of_merged_matrix)):
                    if row_of_merged_matrix[position] == 0 and position not in all_merged_indexes:
                        sum_of_count += outliers_count[position]
                        all_merged_indexes.append(position)
                significant_sequences_count.append(sum_of_count)

    dataframe1 = pd.DataFrame({"Sequence": significant_sequences, "Count": significant_sequences_count,
                               "Mother clone distance": significant_sequences_mother_clone_dist,
                               "Length": significant_sequences_length, "V-genes": significant_v_genes})
    return dataframe1


def Kevinsfullmerge(outliers_inter_allele_1, outliers_inter_allele_2, outliers_count, mother_clone_dist, v_genes, fr_clones, row, save_path):
    matrix = heatmap_matrix_builder(outliers_inter_allele_1, outliers_inter_allele_2)
    build_heatmap_graph(matrix, None, None)
    st.write("Significant sequences")
    if len(outliers_inter_allele_1) > 0:
        dataframe1 = generate_significant_tables(outliers_inter_allele_1, outliers_count, mother_clone_dist, v_genes, matrix)
        fr_clone = max(np.array(fr_clones)[np.where([i is not None for i in fr_clones])], key=len)  # np.array(fr_clones)[np.where([i is not None for i in fr_clones])[0]][0]
        seqsuse = pd.DataFrame({'count': dataframe1['Count'], 'sequence_length': dataframe1['Length'], 'sequence': dataframe1['Sequence'],'found_substring': np.repeat(1, len(dataframe1['Sequence'])), 'found_independent_Kmer': np.repeat(True, len(dataframe1['Sequence']))})
        seqsuse.reset_index(inplace=True)
        if seqsuse.shape[0] > 1:
            subclone_all = subclone.Subclone(seq_file=seqsuse, fr=None, clone_kmer=None, independentK=None, fr_clone=fr_clone, must_be_identical_percent=0, completed_subclone_table=True, row=row, clone_key=row['Clone File'] + "_merge", save_path=save_path)
        st.write(dataframe1)
        for kept_sequence in dataframe1['Sequence'].to_numpy():
            st.write(kept_sequence)
        ## 04/24 we do not need streamlit interactions for the linux version
        ##ste.download_button(label='subclone file 1', data=convert_df(seqsuse), file_name='significant_merged_table_1.csv', mime='text/csv')
        return dataframe1

"""
def merge(all_subclone_tables, frameworks, replicate_count, fr_clones, row):
    # Merging Steps: Need information about number of frameworks and number of replicates
    number_of_frameworks = len(frameworks)
    number_of_replicates = replicate_count
    all_outliers = []
    all_outliers_count = []
    all_outliers_mother_clone_dist = []
    all_v_genes = []
    # Merge between replicates
    for rep_count_1 in range(0, number_of_replicates):
        for framework_count_all_outliers in range(number_of_frameworks):
            if len(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count_all_outliers]) != 0:
                all_outliers.extend(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count_all_outliers]['seq'].tolist())
                all_outliers_count.extend(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count_all_outliers]['sum of counts'].tolist())
                all_outliers_mother_clone_dist.extend(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count_all_outliers]['mother_clone_distance'].tolist())
                print(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count_all_outliers]['V-gene'])
                all_v_genes.extend(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count_all_outliers]['V-gene'].tolist())
                for rep_count_2 in range(rep_count_1, number_of_replicates):
                    st.write('Replicate ' + str(rep_count_1 + 1) + ' Merge With Replicate ' + str(rep_count_2 + 1))
                    outliers_inter_allele_1 = []
                    outliers_inter_allele_2 = []
                    outliers_inter_allele_1_count = []
                    outliers_inter_allele_2_count = []
                    outliers_inter_allele_1_mother_clone_dist = []
                    outliers_inter_allele_2_mother_clone_dist = []
                    outliers_inter_allele_1_v_genes = []
                    outliers_inter_allele_2_v_genes = []
                    for framework_count in range(number_of_frameworks):
                        outliers_inter_allele_1.extend(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count]['seq'].tolist())
                        outliers_inter_allele_2.extend(
                            all_subclone_tables[rep_count_2 + number_of_replicates * framework_count]['seq'].tolist())
                        outliers_inter_allele_1_count.extend(
                            all_subclone_tables[rep_count_1 + number_of_replicates * framework_count]['sum of counts'].tolist())
                        outliers_inter_allele_2_count.extend(
                            all_subclone_tables[rep_count_2 + number_of_replicates * framework_count]['sum of counts'].tolist())
                        outliers_inter_allele_1_mother_clone_dist.extend(
                            all_subclone_tables[rep_count_1 + number_of_replicates * framework_count][
                                'mother_clone_distance'].tolist())
                        outliers_inter_allele_2_mother_clone_dist.extend(all_subclone_tables[rep_count_2 + number_of_replicates * framework_count]['mother_clone_distance'].tolist())
                        outliers_inter_allele_1_v_genes.extend(all_subclone_tables[rep_count_1 + number_of_replicates * framework_count]['V-gene'].tolist())
                        outliers_inter_allele_2_v_genes.extend(all_subclone_tables[rep_count_2 + number_of_replicates * framework_count]['V-gene'].tolist())
                    if len(outliers_inter_allele_1) == 0 or len(outliers_inter_allele_2) == 0:
                        st.write("No Comparison Possible. At least one group does not have outliers found")
                    else:
                        Kevinsmerge(outliers_inter_allele_1, outliers_inter_allele_2, rep_count_1, rep_count_2)
                        data1 = {"Sequence": outliers_inter_allele_1, "Count": outliers_inter_allele_1_count,
                                 "Distance to Mother Clone": outliers_inter_allele_1_mother_clone_dist,
                                 "V-gene": outliers_inter_allele_1_v_genes}
                        data2 = {"Sequence": outliers_inter_allele_2, "Count": outliers_inter_allele_2_count,
                                 "Distance to Mother Clone": outliers_inter_allele_2_mother_clone_dist,
                                 "V-gene": outliers_inter_allele_2_v_genes}
                        dataframe1 = pd.DataFrame(data1)
                        dataframe2 = pd.DataFrame(data2)
                        st.write("Replicate " + str(rep_count_1 + 1) + " Table")
                        st.write(dataframe1)
                        st.write("Replicate " + str(rep_count_2 + 1) + " Table")
                        st.write(dataframe2)
"""
    #Kevinsfullmerge(all_outliers, all_outliers, all_outliers_count, all_outliers_mother_clone_dist, all_v_genes, fr_clones, row)


def get_all_values_across_all_sample_for_column(all_subclone_tables, column_name):
    return [item for sublist in [outlier[column_name].tolist() for outlier in all_subclone_tables if len(outlier) != 0] for item in sublist]


def merge(all_subclone_tables, fr_clones, row, save_path):
    all_outliers = get_all_values_across_all_sample_for_column(all_subclone_tables, 'seq')
    all_outliers_count = get_all_values_across_all_sample_for_column(all_subclone_tables, 'sum of counts')
    all_outliers_mother_clone_dist = get_all_values_across_all_sample_for_column(all_subclone_tables, 'mother_clone_distance')
    all_v_genes = get_all_values_across_all_sample_for_column(all_subclone_tables, 'V-gene')
    #  breakpoint()
    Kevinsfullmerge(all_outliers, all_outliers, all_outliers_count, all_outliers_mother_clone_dist, all_v_genes, fr_clones, row, save_path)


