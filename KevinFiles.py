from bdb import set_trace
import streamlit as st
#import neatbio.sequtils as utils
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from collections import Counter
import base64, zlib
# data pkgs
import matplotlib.pyplot as plt
import matplotlib
from Bio.SeqRecord import SeqRecord
import numpy as np
import pdb
import gzip
import shutil
from itertools import islice
import pandas as pd
st.set_page_config(layout="wide")
from streamlit import caching

import altair as alt  
import bokeh
from scipy.cluster.hierarchy import dendrogram, linkage
import pickle
import re
from Bio.Seq import Seq
import sys
from view_alignments_0326 import *
import time


with st.form(key='prep'):

    seq_file = st.file_uploader("Upload File",type = [".tsv",".txt",".csv"])

    clone = st.text_input("Clone sequence ")
    reorder_sequences =st.checkbox('reorder the sequences')
    submit_button = st.form_submit_button(label='Submit')

if seq_file is not None:
    aln=list()
    aln_pd= pd.read_csv(seq_file)
    for ind, row in aln_pd.iterrows():

        aln.append(SeqRecord(seq=row['sequence'],name=str(row['closest match']),description=str(row['numberobserved']),id=str(row['outliers']),annotations={'seq':row['sequence'],'insert':row['inserts']}))

    tallysum=view_alignment(aln, clone,3,reorder=reorder_sequences)
    