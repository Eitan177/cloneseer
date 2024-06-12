import pdb
import time

from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import numpy as np
import re
import pandas as pd
import streamlit as st
import os
from view_alignments_0401 import view_alignment
from pdb import set_trace


def replace_char_at_index(org_str, index, replacement):
    """
    Replace character at index in string org_str with the
    given replacement character.
    """
    new_str = org_str
    if index < len(org_str):
        new_str = org_str[0:index] + replacement + org_str[index + 1:]
    return new_str


def remove_nonletter(object_align):
    """
    """
    # TODO: This seems to be failing to remove characters that are not letters
    # 01/10 changed seqA and seqB
    seqA = object_align[0][0]
    seqB = object_align[0][1]
    seqB_for_file = seqB
    ind_pass = list()
    skips = re.finditer('-', seqA)
    counter = 0
    for jj in skips:
        seqB = replace_char_at_index(seqB, jj.start() - counter, '')
        ind_pass.append(jj.start() - counter)
        counter += 1
    return seqB, seqB_for_file, ind_pass


def process_seqs(lymphotrack_read_file):
    """
    Takes in the uploaded file object as the lymphotrack read file and processes it in exact order of the file.
    The unique files alternate between read info and the sequences. Method parses a column out for the count and
    sequence length.
    Creates a dataframe with columns count, sequence length and sequence
    """

    lymphotrack_file = pd.read_table(lymphotrack_read_file, header=None)
    lymphotrack_read_info = lymphotrack_file.iloc[np.arange(0, lymphotrack_file.shape[0], 2)]
    lymphotrack_read_sequence = lymphotrack_file.iloc[np.arange(1, lymphotrack_file.shape[0], 2)]
    num_seen = lymphotrack_read_info[0].apply(lambda x: int(re.sub('_.+', '', re.sub('.+Tcount', '', str(x)))))
    len_seq = lymphotrack_read_info[0].apply(lambda x: int(re.sub('.+length', '', str(x))))
    
    lymphotrack_read_df = pd.DataFrame({'count': np.array(num_seen), 'sequence_length': np.array(len_seq), 'sequence': lymphotrack_read_sequence[0]})
    
    return lymphotrack_read_df


class Subclone:
    """
    This is the subclone creation class
    """

    def __init__(self, seq_file, fr, clone_kmer, independentK, fr_clone,input_nanograms, must_be_identical_percent, completed_subclone_table, row, clone_key, save_path,aggregate):
        """
        Initialization method.

        Assigns the following variables as initial part of initialization
            * self.cluster: Dataframe for the input file containing count, sequence length and the actual sequence. Currently always runs process_seqs
            * self.input_file_gui: Input file will exist only for a merge. Otherwise will be false
            * self.fr_total: Matches to the clone
            * self.cluster_clone: Dataframe but only with columns containing matches to a kmer.
            * self.input_nanograms: DNA input in nanograms for the sample
            * self.must_be_identical: [0.75, 0.85, 0.89]. Tinker with this number may be necessary in the main file.
            * self.clone_kmer: All of the kmer splits of the original clone. Was originally self.res
            * self.fr_clone: All of the clones inputted into the original box
            * self.fr: [3, 2, 1]
            * self.see_abundant: False
            * self.int(max(closest_match)+1)_seqs: list()
            * self.independent: independentK cursor
            * self.sample_id: sample ID
            * self.sample_replicate_number: sample replicate number
            * self.clone_id: The name of the clone
            * self.save_path: Path where item will be saved
        """
        if must_be_identical_percent is not None:
            self.must_be_identical_percent = must_be_identical_percent
        else:
            self.must_be_identical_percent = 0
        self.input_nanograms = input_nanograms    
        self.clone_kmer = clone_kmer
        self.fr_clone = fr_clone
        self.fr = fr
        self.see_abundant = False
        self.fr_subclone_seqs = list()
        self.fr_subclone_seqs_unfiltered = list()
        self.independentK = independentK
        self.sample_id = row['Sample ID']
        self.sample_replicate_number = row['Replicate Number']
        self.clone_id = clone_key
        self.save_path = save_path
        self.aggregate=aggregate

        if self.aggregate:
            self.cluster = process_seqs(seq_file)
            table_w_sequence_count_sorted_name= self.save_path + "/" + self.sample_id + "_Rep" + str(self.sample_replicate_number) + "_FR" + str(self.fr) + "_" + self.clone_id + "_Subclonetable_Sorted.xlsx"
            if os.path.exists(re.sub('[.]xlsx','filtered.xlsx',table_w_sequence_count_sorted_name)):
                self.fr_subclone_seqs =pd.read_excel(re.sub('[.]xlsx','filtered.xlsx',table_w_sequence_count_sorted_name))
            if os.path.exists(table_w_sequence_count_sorted_name):    
                self.fr_subclone_seqs_unfiltered=pd.read_excel(table_w_sequence_count_sorted_name) 
        else:
            if not completed_subclone_table:
                self.cluster = process_seqs(seq_file)
                self.input_file_gui = True
                self.fr_total, self.cluster_clone = self.process_match(self.cluster)
            else:
                self.already_processed(seq_file)
            self.first(fr, clone_kmer, independentK, fr_clone)

    def first(self, fr, clone_kmer, independentK, fr_clone):
        """
        Class is visible in order to override fr, clone_kmer, independentK, fr_clone
        self.see_abundant and self.fr_subclone_seqs are also reset.

        """
        self.clone_kmer = clone_kmer
        self.fr_clone = fr_clone
        self.fr = fr
        self.see_abundant = False
        self.fr_subclone_seqs = pd.DataFrame(columns=['seq', 'sum of counts', 'mother_clone_distance', 'v-genes'])
        self.fr_subclone_seqs_unfiltered = pd.DataFrame(columns=['seq', 'sum of counts', 'mother_clone_distance', 'v-genes'])
        self.independentK = independentK

        # st.write("Framework " + str(self.fr) + " kmer match in progress")
        # st.write("Formatting Framework " + str(self.fr) + " sequences that match kmer from clone. Currently processing " + str(self.cluster_clone.shape[0]) + " sequences.")
        self.clean_cluster()
        self.create_visual()
        
        # col_a, col_b, col_c = st.columns(3)
        # print(self.fr_total)
        # col_a.write('match to clone total ' + str(self.fr_total))
        # col_b.write('total ' + str(sum(self.cluster['count'])))
        # col_c.write('clone fraction ' + str(np.round(self.fr_total / sum(self.cluster['count']), 4)))

    def clean_cluster(self):
        sequence_from_cluster = self.cluster_clone['sequence']
        length_from_cluster = self.cluster_clone['sequence_length']
        count_from_cluster = self.cluster_clone['count']
        filtered_data = []
        start = time.time()
        for sequence_number in range(len(sequence_from_cluster)):
            ##local alignment for fr and clone discrepancies
            if (4-self.fr)==round(len(self.fr_clone)/100):
                alignment = pairwise2.align.globalms(self.fr_clone, sequence_from_cluster.iloc[sequence_number], 1, -1, -5, 0)
            else:
                alignment = pairwise2.align.localms(self.fr_clone, sequence_from_cluster.iloc[sequence_number],2, -3, -4, -2)  
            seqB = alignment[0][1]
            #print(time.time())
            if seqB.lstrip("-").find("------------------") == -1:
                filtered_data.append([seqB, length_from_cluster.iloc[sequence_number], count_from_cluster.iloc[sequence_number]])
        self.cluster_clone = pd.DataFrame(filtered_data, columns=['sequence', 'sequence_length', 'count'])

    def process_match(self, ses):
        """
        Creates two new columns called found_substring which is the number of times any clone kmer is seen and
        found_independent_Kmer is true when the found_substring value is greater than the independent K value set in the GUI

        ses_count_sum is the sum of the found substrings all found independent kmer that were true.

        Process performed by first finding all unique starting points for the kmer split to create multiples. self.clone_kmer[0] are the k lengths
        Ind_single_kmer is responsible for grabbing all the multiples starting from independent_k and this processes all the clone_kmers
        """
        ses['found_substring'] = 0  # np.array()
        ses['found_independent_Kmer'] = False  # np.array()
        # pro = st.progress(0)
        count_col = 0
        for independent_Ks in np.arange(0, len(self.clone_kmer[0])):  # This is an array from 0 to K
            ses['found_substring'] = 0
            for ind_single_kmer in np.arange(independent_Ks, len(self.clone_kmer), len(self.clone_kmer[0])):  # Multiples of K starting from independent K
                single_kmer = self.clone_kmer[ind_single_kmer]
                ses['found_substring'] = ses['found_substring'] + ses.sequence.str.contains(single_kmer)
                count_col += 1
            ses['found_independent_Kmer'] = ses['found_independent_Kmer'] + (self.independentK <= ses['found_substring'])
            # pro.progress(count_col / len(self.clone_kmer))
        
        if ses.index[0] != 0:
            ses.reset_index(inplace=True)
            ses.rename(columns={0: 'seq'}, inplace=True)
        ses_count_sum = sum(ses[ses['found_independent_Kmer']]['count'])
        ses = ses[ses['found_independent_Kmer']].sort_values('count', ascending=False)  # [0:5000]
        
        return ses_count_sum, ses

    def get_partclus_and_aln(self, partclus, clone):
        """
        Performs a pairwise alignment on the clone sequence and uploads into variable all_pairwise_alignment_objects which contains an array
        that creates [Alignment(seqA=, seqB=, score=, start=, end=)].

        Filters out by

        """
        abundant_seq_or_fr_clone = clone
        aligned_counter = 0
        all_pairwise_alignment_objects = list()
        # cols_show = st.columns(10)

        # Appends pairwise objects into the all_pairwise_alignment_objects.
        # Appears as [Alignment(seqA=, seqB=, score=, start=, end=)]
        for x in partclus['sequence']:
            if (4-self.fr)==round(len(self.fr_clone)/100):
                all_pairwise_alignment_objects.append(pairwise2.align.globalms(abundant_seq_or_fr_clone, x, 1, -1, -5, 0))
            else:
                all_pairwise_alignment_objects.append(pairwise2.align.localms(abundant_seq_or_fr_clone, x,1, -1, -4, -1)) 

            #all_pairwise_alignment_objects.append(pairwise2.align.globalms(abundant_seq_or_fr_clone, x, 1, -1, -5, 0))
            aligned_counter += 1

        # First tries to remove all dashes from sequences. Adds a new column called score from the alignment. Creates a list
        # of SeqRecord files that
        New_seqBs = [remove_nonletter(aaa) for aaa in all_pairwise_alignment_objects]
        if not New_seqBs:
            match_len = len(clone)
        else:
            match_len = np.min((len(clone), np.max([len(x[1]) for x in New_seqBs])))
        
        partclus['score'] = [ii[0].score for ii in all_pairwise_alignment_objects]  # Gathers the scores from the alignment
        aln = [SeqRecord(seq=jjj[0][0], name=str(jjj[2]), description=str(jjj[1]), id=str(jjj[3]), annotations={'seq': jjj[0][1], 'insert': jjj[0][2]}) for jjj in zip(New_seqBs, partclus['count'].tolist(), partclus['score'], partclus['sequence'].tolist())]
        aln = [aln[jjj] for jjj in np.argsort(partclus['score'])]
        partclus = partclus.sort_values(['score'])
        # 12/23 remove
        #score_to_use = [pairwise2.align.localms(gg.id, clone, 1, -1, -5, -1, score_only=True) for gg in aln]
        
        # 12/23 make this new score_to_use
        score_to_use = np.round([pairwise2.align.localms(gg.id, clone, 1, -1.01, -5.0001, -1.000001, score_only=True) for gg in aln], 7)
        sequence_lengths=[len(ff.id.replace('-','')) for ff in aln ]
        mismatch_penalty = np.floor(100 - (np.array(score_to_use) * 100) % 100) % 100
        gap_open_penalty = (100-np.floor((np.array(score_to_use)*10000) % 100)) % 100
        gap_extension = (100-np.ceil((np.array(score_to_use)*1000000) % 100)) % 100
        matching_nuc_len= np.floor(score_to_use)+mismatch_penalty + (gap_open_penalty*5) + gap_extension
        matching_nucleotide_num = match_len - mismatch_penalty - gap_open_penalty - gap_extension
        
        

        # aln_sort_by_score = [hh[0] for hh in zip(aln, score_to_use) if not hh[1] <= (len(clone) * self.must_be_identical_percent)]

        # 06/19 comment this out
        #aln_sort_by_score = [hh[0] for hh in zip(aln, matching_nucleotide_num) if not hh[1] < match_len * self.must_be_identical_percent]  # TODO: Delete after using as test

        # 06/19 try this formulat
        if aln == []:
            return pd.DataFrame(),[]
        else:
            if self.fr==3 and len(self.fr_clone)>200:
                print('fr is three and len clone is greater than 200')
                #breakpoint()
            #aln_sort_by_score = [hh[0] for hh in zip(aln, matching_nuc_len,sequence_lengths) if not hh[1] <  hh[2] * self.must_be_identical_percent]  # TODO: Delete after using as test
            aln_sort_by_score = [hh[0] for hh in zip(aln, matching_nuc_len,sequence_lengths) if not hh[1] <  min(hh[2],len(self.fr_clone)) * self.must_be_identical_percent]
            if aln_sort_by_score == []:
                return pd.DataFrame(),[]
            # st.write('keep ', str(len(aln_sort_by_score)) + ' kmer matches, filter ' + str(np.sum([s <= (len(clone) * self.must_be_identical_percent) for s in score_to_use])) + ' after pairwise alignment')
            else:
            
                #06/19 comment
                #partclus = partclus[matching_nucleotide_num >= match_len * self.must_be_identical_percent]
                
                partclus_ret =pd.DataFrame([hh[0] for hh in zip(partclus.values.tolist(), matching_nuc_len,sequence_lengths) if not hh[1] <  min(hh[2],len(self.fr_clone)) * self.must_be_identical_percent])
                partclus_ret.columns=partclus.columns
                
                #06/19 add
                aln = aln_sort_by_score
                #if self.fr==3 and len(self.fr_clone)>200:
                #    breakpoint()
                return partclus_ret, aln  # , clip left, clip right

    def create_visual(self):
        # with st.spinner('formatting framework sequences for ' + str(self.cluster_clone.shape[0]) + ' sequences and calculating alignment scores'):
        part1, aln1 = self.get_partclus_and_aln(self.cluster_clone, self.fr_clone)
        print("Framework " + str(self.fr))
        if len(aln1) > 0:
            tally_sum, subclone_seqs_fr_unfiltered,subclone_seqs_fr = view_alignment(aln1, self.fr_clone, self.fr, self.sample_id, self.sample_replicate_number, self.clone_id, self.save_path, self.cluster['count'].sum(),self.input_file_gui)
            self.fr_total = tally_sum
            self.fr_subclone_seqs_unfiltered = subclone_seqs_fr_unfiltered
            self.fr_subclone_seqs = subclone_seqs_fr

    def already_processed(self, seq_file):
        self.cluster = seq_file
        self.cluster_clone = self.cluster
        self.fr_total = self.cluster_clone['count'].sum()
        self.input_file_gui = False
