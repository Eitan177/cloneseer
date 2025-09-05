from operator import sub

import matplotlib.pyplot
import numpy as np
from bokeh.io import export_png
from bokeh.plotting import figure,show,output_file, save
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
import bokeh
from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering
# from scipy.cluster.hierarchy import linkage
from fastcluster import linkage
from scipy.spatial.distance import pdist, squareform
import scipy
import streamlit as st
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis
from bokeh import io
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import pdb



@st.experimental_memo
def getlinkage(uniqdf, distancemetric):  # , optimal_ordering):
    return (linkage(uniqdf, method=distancemetric))  # ,optimal_ordering=optimal_ordering))


@st.cache()
def convert_df(fr_df):
    return fr_df.to_csv(index=False).encode('utf-8')


def make_alignmentplotwithcluster(aln_pd, useconsensus, sample_id, sample_replicate_number, clone_id, save_path, fr=1,
                                  fontsize="9pt", plot_width=750, see_seq=False, reorder=False, input_file_gui=True,full_plot=False):
    """Bokeh sequence alignment view"""
    aln = list()
    for ind, row in aln_pd.iterrows():
        aln.append(
            SeqRecord(seq=row['sequence'], name=str(row['closest match']), description=str(row['numberobserved']),
                      id=str(row['outliers']),
                      annotations={'seq': row['sequence'],'homology_against_mother_clone':row['homology_against_mother_clone'], 'insert': row['inserts'], 'alnscores': row['alnscores']}))

    
    consensus_sequence = useconsensus

    # seqs = [[rec.seq] * int(np.ceil(int(rec.description))) for rec in (aln)]
    forfileseqs = [rec.annotations['seq'] for rec in aln]
    forfileinsert = np.array([rec.annotations['insert'] for rec in aln])
    seqs = np.array([rec.seq for rec in aln])
    seqs=np.array(["".join(seq) for seq in seqs])
    unformatted_seqs=aln_pd['sequenceNotformatted'].tolist()
    counts = [int(np.ceil(int(rec.description))) for rec in aln]
    alnscores = np.array([rec.annotations['alnscores'] for rec in aln])
    homology_against_mother_clone=np.array([rec.annotations['homology_against_mother_clone'] for rec in aln])
    closest_match = np.array([float(rec.name) for rec in aln])
    subclonecall = np.array([float(rec.id) for rec in aln])
    
    #neworder=np.lexsort([-1*np.array(counts),homology_against_mother_clone])
    #if clone_id == 'PositiveControl_FR1_filtered' and fr==2:
    #    breakpoint()
    neworder = np.argsort(seqs)
    forfileseqs = np.array(forfileseqs)[neworder].tolist()
    forfileinsert = forfileinsert[neworder]
    seqs = seqs[neworder]

    unformatted_seqs=np.array(unformatted_seqs)[neworder]

    counts = np.array(counts)[neworder].tolist()
    alnscores = alnscores[neworder]
    homology_against_mother_clone=homology_against_mother_clone[neworder]
    closest_match = closest_match[neworder]
    subclonecall = subclonecall[neworder]

    
    
    if int(max(closest_match)+1) <= 11 and int(max(closest_match)+1) >= 3:
        subclonecolors = Viridis[int(max(closest_match)+1)]
    elif int(max(closest_match)+1) > 11 and int(max(closest_match)+1) < 14:

        subclonecolors = Viridis[11] + Cividis[3][0:int(max(closest_match)+1) - 11]
    elif int(max(closest_match)+1) >= 14 and int(max(closest_match)+1) <= 22:
        subclonecolors = Viridis[11] + Cividis[int(max(closest_match)+1) - 11]
    elif int(max(closest_match)+1) > 22 and int(max(closest_match)+1) <= 24:
        subclonecolors = Magma[11] + Viridis[11] + Cividis[3][0:int(max(closest_match)+1) - 22]

    elif int(max(closest_match)+1) >= 24 and int(max(closest_match)+1) <= 32:
        subclonecolors = Magma[11] + Viridis[11] + Cividis[int(max(closest_match)+1) - 22]
        
    elif int(max(closest_match)+1) >= 33 and int(max(closest_match)+1) < 240:
        cc = Magma[11] + Viridis[11] + Cividis[10] + Inferno[11] + Plasma[11] + Cividis[10] + Magma[11] + Viridis[11] + \
             Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + \
             Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11]
        subclonecolors = cc[0:int(max(closest_match)+1)]

    else:
        subclonecolors = np.tile(Cividis[11], 500)[0:int(max(closest_match)+1)]
    st.write('a maximum '+str(int(max(closest_match)+1)) + ' subclones identified')

    colors = get_colors_consensus(seqs, consensus_sequence)

    gg = np.array(colors).reshape((-1, len(seqs[0])))
    
    subclonecolormatch = np.repeat([subclonecolors[int(ii)] for ii in closest_match], 5).reshape(-1, 5)
    
    gg = np.hstack((gg, subclonecolormatch))


    
    for_ordering, countsub = np.unique(subclonecolormatch, return_counts=True)
    countsub_sort_ind = np.argsort(-countsub)
    for_ordering = for_ordering[countsub_sort_ind]

    fororder_list = []

    if reorder:
        for ind, jj in pd.DataFrame(for_ordering).iterrows():
            fororder_list.append(np.where(np.all(subclonecolormatch == np.array(jj), axis=1))[0])
        contiguous_color_order = np.hstack(fororder_list)
        gg = gg[contiguous_color_order]
        forfileinsert = forfileinsert[contiguous_color_order]
        alnscores = alnscores[contiguous_color_order]
        homology_against_mother_clone=homology_against_mother_clone[contiguous_color_order]
        closest_match = closest_match[contiguous_color_order]
        counts = np.array(counts)[contiguous_color_order]
        seqs = seqs[contiguous_color_order]
        unformatted_seqs=unformatted_seqs[contiguous_color_order]
        subclonecall = subclonecall[contiguous_color_order]
        subclonecolormatch = subclonecolormatch[contiguous_color_order]


   
    locofinsertion2=[]

    locofinsertion2 = [[int(y) for y in np.unique(np.array( jj.replace('[', '').replace(']', '').replace(' ','').split(',')))] for jj in forfileinsert[np.where(forfileinsert != '[]')]]
    for jj,hh in zip(np.where(forfileinsert != '[]')[0], locofinsertion2):
        gg[jj, hh] = 'gray'


    gg_subc = gg.copy()
    seqs_subc = seqs.copy()
    unformatted_seqs_subc=unformatted_seqs.copy()

    for ii in set(closest_match):
        gg_subc[closest_match == ii] = gg[np.logical_and(closest_match == ii, subclonecall > 0)]
        seqs_subc[closest_match == ii]= seqs[np.logical_and(closest_match == ii, subclonecall > 0)]
        unformatted_seqs_subc[closest_match == ii]=unformatted_seqs[np.logical_and(closest_match == ii, subclonecall > 0)]

    seqs_subc_unique = pd.DataFrame({'seq':seqs_subc}).drop_duplicates()

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

    col_score_kms = st.columns(11)

    # setting distance_threshold=0 ensures we compute the full tree.
    dendo, axis = plt.subplots(1, 1, figsize=(5, 15))

    uniq, counts = unique_nnn, countgg

    if uniq.shape[0] > 1:
        uniqdf = pd.DataFrame(uniq)

        uniqdf = np.array(uniqdf)
        countsforlink = counts.copy()

        if full_plot==False and input_file_gui:
             if homology_against_mother_clone[np.argmax(countsforlink)] == np.min(homology_against_mother_clone):# and np.max(countsforlink) > 100:   

                 countsforlink[np.argmax(countsforlink)] = 1#np.ceil((sum(counts) - max(counts)) / 100)




        dfwcountnomain = np.repeat(uniqdf, countsforlink, axis=0)
        mismatchmutatedspotsuse = np.sum(dfwcountnomain != 0, axis=0) / dfwcountnomain.shape[0] > 0.01
        indelmutatedspotsuse = np.sum(dfwcountnomain == 5, axis=0) / dfwcountnomain.shape[0] > 0.001
        indexfirstletternotindel = np.where(~indelmutatedspotsuse)[0][0]
        if indexfirstletternotindel > 0:
            indelmutatedspotsuse[0:indexfirstletternotindel] = False
        mutatedspotsuse = mismatchmutatedspotsuse + indelmutatedspotsuse



        countview_fordf = np.repeat(counts, countsforlink)
        dfwcounts = np.repeat(uniqdf, countsforlink, axis=0)

 
        uniqtocolors = np.repeat(ugg, countsforlink, axis=0)
        justsubc_tocolors = np.repeat(gg_subc, countsforlink, axis=0)
        seqs_subc = np.repeat(seqs_subc, countsforlink, axis=0)
        unformatted_seqs_subc = np.repeat(unformatted_seqs_subc, countsforlink, axis=0)
        counts=np.repeat(counts, countsforlink, axis=0)
        forfileseqs = np.repeat(forfileseqs, countsforlink)
        forfileinsert = np.repeat(forfileinsert, countsforlink)
        alnscores = np.repeat(alnscores, countsforlink)
        homology_against_mother_clone=np.repeat(homology_against_mother_clone, countsforlink)
        closest_match = np.repeat(closest_match, countsforlink)
        countview_alnscores = alnscores
        
        # uniqtocolors[np.where(countview_fordf==np.max(countview_fordf)),0:(uniqtocolors.shape[1]-5) ]=['purple']
        # justsubc_tocolors[np.where(countview_fordf==np.max(countview_fordf)),0:(justsubc_tocolors.shape[1]-5) ]=['purple']
        if full_plot==False:
            if homology_against_mother_clone[np.argmax(counts)] == np.min(homology_against_mother_clone):
                homology_min_purple_ind = np.argmax(counts)
            else:
                homology_min_purple_ind = np.argmin(homology_against_mother_clone)  
            uniqtocolors[homology_min_purple_ind, 0:(uniqtocolors.shape[1] - 5)] = ['purple'] 




        seqs_for_view = np.repeat(seqs, countsforlink)

    else:
        seqs_for_view = seqs
        uniqtocolors = ugg
        justsubc_tocolors = gg_subc
    ## 03/26
    incrementforview = int(np.round(seqs_for_view.shape[0] / 1000))
    ##incrementforview = 1
    ## 03/26

    if incrementforview > 1:
        indsview = np.arange(0, seqs_for_view.shape[0], incrementforview)

        print([np.where(justsubc_tocolors[:,-1]==jj)[0][0] for jj in np.unique(subclonecolormatch)])
       
        if full_plot==False:
            indsview = np.append(indsview, homology_min_purple_ind)
        indsview = np.append(indsview,np.array([np.where(justsubc_tocolors[:,-1]==jj)[0][0] for jj in np.unique(subclonecolormatch)]))#np.argmax(countview_alnscores))
        indsview = np.sort(np.unique(indsview))

        indsview = np.sort(indsview)
        seqs_for_view = seqs_for_view[indsview]
        uniqtocolors = uniqtocolors[indsview]
        alnscores = alnscores[indsview]
        closest_match = closest_match[indsview]
        justsubc_tocolors = justsubc_tocolors[indsview]
        seqs_subc=seqs_subc[indsview]
        unformatted_seqs_subc=unformatted_seqs_subc[indsview]

        forfileseqs = forfileseqs[indsview]
        forfileinsert = forfileinsert[indsview]

    #st.download_button(label='sequences in alignment shown', data=convert_df(pd.DataFrame(
    #    {'sequence': seqs_for_view, 'sequenceNotformatted': forfileseqs, 'inserts': forfileinsert.tolist(),
    #     'group': closest_match, 'alnscore': alnscores})), file_name='alignmenview_sequences_' + str(fr) + '.csv',
    #                   mime='text/csv', key=np.random.rand())

    colors_for_view = np.flip(uniqtocolors, axis=1).ravel().tolist()
    just_subc_colors_for_view = np.flip(justsubc_tocolors, axis=1).ravel().tolist()
    
    N = len(seqs_for_view[0])
    S = len(seqs_for_view)
    width = .4

    # 03/26
    # x = np.arange(1,N+1)
    x = np.arange(1, N + 6)
    y = np.arange(0, S, 1)

    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()

    # use recty for rect coords with an offset
    recty = gy + .5
    # recty2 = gy2+0.5
    # recty3 = gy3+0.5
    h = 1 / S

    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, colors=colors_for_view[::-1]))  # text=text
    subc_source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, colors=just_subc_colors_for_view[::-1]))  # text=text

    if see_seq:
        plot_height = len(seqs_for_view) * 10 + 50

    ## 03/26
    # x_range = Range1d(0,N+1, bounds='auto')
    x_range = Range1d(0, N + 6, bounds='auto')

    view_range = (0, N)

    tools = "xpan, xwheel_zoom, reset, save"

    # sequence text view with ability to scroll along x axis

    p1 = figure(title=None, width=plot_width, height=700,
                x_range=x_range, y_range=(0, S), tools=tools,  # "xpan,reset",
                min_border=0, toolbar_location='below')  # , lod_factor=1)
    glyph = Text(x="x", y="y", text="text", text_align='center', text_color="black",\
                 # text_font="monospace",\
                 text_font_size=fontsize)

    rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.4)

    p2 = figure(title=None, width=plot_width, height=700,
                x_range=x_range, y_range=(0, S), tools=tools,  # "xpan,reset",
                min_border=0, toolbar_location='below')  # , lod_factor=1)
    glyph2 = Text(x="x", y="y", text="text", text_align='center', text_color="black",
                  # text_font="monospace",
                  text_font_size=fontsize)

    rects2 = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                  line_color=None, fill_alpha=0.4)

    p2.add_glyph(source, glyph)
    p2.add_glyph(source, rects)

    p2.grid.visible = False
    p2.xaxis.major_label_text_font_style = "bold"
    p2.yaxis.minor_tick_line_width = 0
    p2.yaxis.major_tick_line_width = 0

    p1.add_glyph(subc_source, glyph2)
    p1.add_glyph(subc_source, rects2)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p1, p2]], toolbar_location='below')
    if full_plot==False:
        suffix="_Clean_V_Dirty.html"
    else:
        suffix="_withReferenceInProportion_Clean_V_Dirty.html"    
    if fr is None:

        output_file(filename=save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_" + str(clone_id) + suffix, title="clean html")
        save(p)
    else:
        output_file(filename=save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + str(clone_id) +suffix)
        
        save(p)
        filename=save_path + "/" + sample_id + "_Rep" + str(sample_replicate_number) + "_FR" + str(fr) + "_" + str(clone_id) +suffix

    unformatted_seqs_subc_unique = pd.DataFrame({'notformatted':unformatted_seqs_subc}).drop_duplicates()
    return unformatted_seqs_subc_unique #seqs_subc_unique


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
