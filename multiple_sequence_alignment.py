import os, io, random
import pdb
import string
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
#import igblast_parser
import panel as pn
import panel.widgets as pnw
pn.extension()

from io import StringIO
import sys
import subprocess

from bokeh.plotting import figure, output_file,save
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import export_png
from KevinsFunction import *    

def blastp_get_top_hits(input_fp, db_fp, organism='human', element='V',debug = True, top = '3'):
    """
    Function that parses the hits table from blastp (bin/igblastp) into pandas dataframe. The hits table shows the top 3 hits (default) of the germline seuqences that are with highest identity % to the query seuqence.
    :param input_fp: string, file path of the input fasta file
    :param db_fp:  string, file path of the serach database, i.e. the value for `-germline_db_V` argument
    :param organism: string, optimal, default = human
    :param debug: bool, optional, default = True. If True, it will print the raw data from cmd (before parsing into dataframe)
    :param top: int, optional default = 3. The top N hits to extract.
    :return: a tuple of (df, the best matching germline allele name, the best matching germline allele identity score)
    """

    # igblastp can only search V database
    if element == 'V':
        usedb='-germline_db_V'
    elif element == 'J':
        usedb='-germline_db_J'
    else:
        usedb='-germline_db_D'    

    cmd = ['/Users/eitanhalper-stromberg/Downloads/ncbi-igblast-1.20.0-src/c++/Clang1300-DebugMT64/bin/igblastp', usedb, db_fp, '-query', input_fp, '-organism', organism, '-outfmt', '7','-num_alignments_V', top]

    a = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    ## set igdata to ncbi directory
    ## export IGDATA=/Users/eitanhalper-stromberg/Documents/0215sublcone_simplified/subclone_simplified/ncbi-igblast-1.20.0/
    
    # show output in notebook
    #a.wait()
    #output = a.stdout.readlines()
    #print(output)

    # parse output into string
    b = StringIO(a.communicate()[0].decode('utf-8'))
    #(b, err) = a.communicate()
    all_data = [x.strip() for x in str(b.getvalue()).split('#')]
    if debug:
        print(all_data)

    top=int(top)    
    try:

        # the hit table is the 2nd from the last in the list above.
        hits=[]
        for line in all_data:
            if line.startswith(str(top) + ' hits found'):
                hits+=line.split('\n')[1:top+1]
        
        # default fields
        fields = 'chain type\tquery id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq.start\tq.ends\t s.start\ts.end\tevalue\tbit score'
        
        # add col header
        hits.insert(0, fields)
        # parse into df
        data = StringIO('\n'.join(hits))
        
        df = pd.read_csv(data, sep='\t')

        # get df, the best matching germline allele & identity
        
        return df, df['subject id'][0], df['% identity'][0]
    except:
        print('None of the seuqenced you provided in the fasta file ')
        return None, None, None

def run_igblast(file_path):
    cmd1 = ["bin/igblastn", "-query", file_path,"-germline_db_V", "database_clean/IGHV_clean.fasta", "-germline_db_D", "database_clean/IGHD_clean.fasta", "-germline_db_J", "database_clean/IGHJ_clean.fasta", "-organism", "human", "-out", "igblastout.txt"]
    #breakpoint()
    igblastoutput = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd='/content/ncbi-igblast-1.22.0')
    igblastoutput.wait()
    
    cmd2=["igblast-parser", "--in", "igblastout.txt", "--out","parsed_igblastout"]
    igblastparseroutput = subprocess.Popen(cmd2, stdout=subprocess.PIPE)
    igblastparseroutput.wait()
    with open("igblastout.txt", "r") as f: 
        lines = f.readlines()

    querysequencename=[]
    percentline=[]
    count_sequences=0
    count_percents=0
    for line in lines: 
        if(re.findall('^Query\=',line) != []):
            querysequencename.append(line.strip())
            count_sequences+=1
        elif(re.findall('^Total\t',line) != []):
            
            percentline.append(re.sub('^Total\t','',line).strip().split('\t'))
            
            count_percents+=1
        if count_sequences != count_percents and (re.findall('^V  ',line) != []):
            #breakpoint()
            percentline.append(['','',re.sub('^.+\/','',re.sub('\)','',line.split(' ')[3])),re.sub('\/.+','',re.sub('\(','',line.split(' ')[3])), '','',  re.sub('%','',line.split(' ')[2])])
            count_percents+=1
    querysequencename=[re.sub('^.+_','',re.sub('^Query\= ','',x)) for x in querysequencename]
    #percentline=[re.sub('^Total\t','',x).split('\t') for x in percentline]
    v_length=[]
    v_match_length=[]
    v_percent_match=[]
    for ii in percentline:
        v_length.append(ii[2])
        v_match_length.append(ii[3])
        v_percent_match.append(ii[6])
       
    igblast_parsed_output_vmatch_percent=pd.DataFrame({'querysequencename':querysequencename,'v_length':v_length,'v_match_length':v_match_length,'v_percent_match':v_percent_match})

    igblast_parsed_output=pd.read_csv("parsed_igblastout.csv")
    igblast_parsed_output['umi']=igblast_parsed_output['umi'].apply(lambda x: re.sub('^.+_','',x))
    
    igblast_parsed_output=pd.merge(igblast_parsed_output,igblast_parsed_output_vmatch_percent,left_on='umi',right_on='querysequencename')
    
    return igblast_parsed_output

def create_multiple_sequence_alignment(table, save_path, file_id):
    make_fasta_file(table, save_path, file_id)
    
    igblast_parsed_output=run_igblast(save_path + "/" + file_id + "_Subclonetable_Sorted_toblast.fasta")
    table['umi']=[str(y) for y in np.arange(len(table.index))]
    
    table=pd.merge(table,igblast_parsed_output)#,on='umi')
    aln = AlignIO.read(save_path + "/" + file_id + "_Subclonetable_Sorted.fasta", 'fasta')
    p = view_alignment(aln, plot_width=3000)
    pn.pane.Bokeh(p)
    output_file(filename=save_path + "/" + file_id + "_Subclonetable_Alignment.html")
    save(p)
    #export_png(p, filename=save_path + "/" + file_id + "_Subclonetable_Alignment.png")
    ## 02/21 get V element from the igblast
    
    table['outliers']=1
    table['sequence']=table['seq']
    
    #data=process_fasta_from_ebi("Vs.fasta", table,id_separator="|",id_separator_index=1,element_to_match="IGHV")
    
    #data=process_fasta_from_ebi("IGHJ.fasta", table,id_separator=">",id_separator_index=1,element_to_match="IGHJ")
    
    #data=process_fasta_from_ebi("IGHD.fasta", data,id_separator=">",id_separator_index=1,element_to_match="IGHD")
    
    #dfv, top_germ_allelev, top_identityv =blastp_get_top_hits(save_path + "/" + file_id + "_Subclonetable_Sorted.fasta", db_fp= '/Users/eitanhalper-stromberg/Documents/0215sublcone_simplified/subclone_simplified/ncbi-igblast-1.20.0/database/Homo_sapiens_clean/IG_dna/IGHV_clean', organism='human', element='V')
    ## 02/21 get J element from the igblast
    #dfj, top_germ_allelej, top_identityj =blastp_get_top_hits(save_path + "/" + file_id + "_Subclonetable_Sorted.fasta", db_fp= '/Users/eitanhalper-stromberg/Documents/0215sublcone_simplified/subclone_simplified/ncbi-igblast-1.20.0/database/Homo_sapiens_clean/IG_dna/IGHJ_clean', organism='human', element='V')
    ## 02/21 get D from the igblast
    #dfd, top_germ_alleled, top_identityd = blastp_get_top_hits(save_path + "/" + file_id + "_Subclonetable_Sorted.fasta", db_fp= '/Users/eitanhalper-stromberg/Documents/0215sublcone_simplified/subclone_simplified/ncbi-igblast-1.20.0/database/Homo_sapiens_clean/IG_dna/IGHD_clean', organism='human', element='V')
   
    #return dfv, dfj, dfd
    
    return table

def make_fasta_file(table, save_path, file_id):
    fasta_file_for_igblast = open(save_path + "/" + file_id + "_Subclonetable_Sorted_toblast.fasta", "w")
    fasta_file_for_view = open(save_path + "/" + file_id + "_Subclonetable_Sorted.fasta", "w")
    fasta_lines_for_igblast = []
    fasta_lines_for_view = []
    for row_count in range(len(table.index)):
        fasta_lines_for_igblast.append(">Seq_" + str(row_count) + "\n")
        fasta_lines_for_view.append(">Seq_" + str(row_count) + "\n")
        fasta_lines_for_igblast.append(str(re.sub('-+','',table["seq"].iloc[row_count])) + "\n")
        fasta_lines_for_view.append(str(table["seq"].iloc[row_count]) + "\n")
    fasta_file_for_igblast.writelines(fasta_lines_for_igblast)
    fasta_file_for_view.writelines(fasta_lines_for_view)
    fasta_file_for_igblast.close()
    fasta_file_for_view.close()



def view_alignment(aln, fontsize="9pt", plot_width=800):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)
    N = len(seqs[0])
    S = len(seqs)
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, width= plot_width, height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, width=plot_width, height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p], [p1]], toolbar_location='below')
    return p


def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs = {'A':'red','T':'green','G':'orange','C':'blue','-':'white','N':'brown'}
    colors = [clrs[i] for i in text]
    return colors


def muscle_alignment(seqs):
    """Align 2 sequences with muscle"""
    filename = 'temp.faa'
    SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=name+'.txt')
    stdout, stderr = cline()
    align = AlignIO.read(name+'.txt', 'fasta')
    return align