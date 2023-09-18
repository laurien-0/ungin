#required imports
import os
import pathlib
import re

import numpy as np
import pandas as pd

import Bio.SeqIO as SeqIO
import Bio.SeqUtils.ProtParam as bp

import matplotlib.pyplot as plt
import seaborn as sns

def write_df_to_fasta(df, output_file):
    """
    Writes sequences stored in a dataframe to FASTA.

    Parameters:
        df (dataframe)       : Dataframe with columns 'ID' and 'Sequence'.
        output_file (string) : Path to save .fasta file.
    """
    with open(output_file,"w") as g:
        for i in range(len(df)):
            gene = df.loc[i,"ID"]
            seq = df.loc[i,"Sequence"]
            g.write(f">{gene}\n{seq}\n")

def parse_uracil_genome(genome):
    """
    Extracts information on protein products of the genome and stores this in a dataframe. 

    Parameters:
        genome (string)   : Absolute or relative path to GenBank record of genome.

    Returns:
        df (dataframe)    : Stores the protein ID, protein sequence and the sequence length for all CDSs.
    """
    prot_id = []
    prot_seq = []
    prot_seqlen = []
    
    with open(genome):
        for record in SeqIO.parse(genome, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    prot_id.append(feature.qualifiers["gene"][0])
                    prot_seq.append(feature.qualifiers["translation"][0])
                    prot_seqlen.append(len(feature.qualifiers["translation"][0]))
    df = pd.DataFrame(list(zip(prot_id,prot_seq,prot_seqlen)),columns =['ID','Sequence','Seq Length']).sort_values(by=['ID'])
    df.index = np.arange(len(df))
    
    return df

def parse_blast_output(df, input_file):
    """
    Extracts information on protein products of the genome and stores this in a dataframe. 

    Parameters:
        df (dataframe)           : Dataframe containing same sequences as were submitted to BLAST in the same order.
        input_file (string)      : Path to BLAST output .txt file.

    Returns:
        df_blasthits (dataframe) : Stores the E-value and sequence identity of the first hit for all sequences with\
                                   at least one significant hit.
    """
    with open(input_file,"r") as g:
        scraped_text = g.readlines()
    p = re.compile(r'Query #\S+: (\S+)')
    q = re.compile(r'%\s+(.+?)\s+(.+?) ')
    i=0
    e_values=[]
    identities=[]
    for line in scraped_text:
        match = p.search(line)
        if match:
            match_2 = q.search(scraped_text[i+5])
            e_values.append(float(match_2.group(1)))
            identities.append(float(match_2.group(2)))
        i += 1
    df_blasthits=df
    df_blasthits['BLAST_evalues'] = e_values
    df_blasthits['BLAST_seqids'] = identities
    df_blasthits=df_blasthits[(df_blasthits["BLAST_evalues"] < 1e-05) & (df_blasthits["BLAST_seqids"] > 35)]
    df_blasthits=df_blasthits[["ID","BLAST_evalues","BLAST_seqids"]]
    
    return df_blasthits

def make_confidence_plots(output_dir,gene_code, plddt, pae, vmin, vmax, scale):
    """
    Creates side by side plot displaying the confidence data output by ESMFold Structure Prediction.\
    pLDDT values are shown as a line graph and PAE values as a heatmap.  Saves plot as a .png file.

    Parameters:
        output_dir (string) : Absolute or relative path to directory to store plots.
        gene_code (string)  : ID of gene which has undergone structure prediction
        plddt (numpy array) : Mean pLDDT value for each amino acid.
        pae (numpy matrix)  : PAE values for each residue compared to every other residue.
        vmin (integer)      : Minimum value of color scale.
        vmax (integer)      : Maximum value of color scale.
        scale (Boolean)     : True/False to indicate whether color scale bar should be included.
    """
    index = np.arange(1, len(plddt)+1)
    fig, axs = plt.subplots(1,2,figsize=(12,5))
    fig.suptitle(gene_code)
    axs[0].plot(index, plddt)
    axs[0].set_title('pLDDT scores')
    axs[0].set(ylim=(0, 100))
    axs[0].text(x=0,y=5,s=f"Mean pLDDT {round(plddt.mean(),2)}")
    sns.heatmap(pae, ax=axs[1], vmin=vmin, vmax=vmax, cbar=scale)
    axs[1].set_title('PAE matrix')
    plt.savefig(f'{output_dir}/{gene_code}.png')
    plt.close()
    
def analyse_confidence(input_dir, df, scale):
    """
    Data wrangling to allow confidence plots to be created for a set of ESM structure predictions.

    Parameters:
        input_dir (string): Absolute or relative path to directory containing ESM structure prediction output.
                            Should contain 2 directories:
                            /pae with pae matrices saved as id_pae.txt
                            /structures with pdb files saved id.pdb
        df (dataframe)    : Dataframe where sequences of proteins are stored.
        scale (Boolean)   : True/False to indicate whether color scale bar should be included.

    Returns:
        df (dataframe)    : Modified to include mean pLDDT scores
    """
    os.chdir(input_dir)
    if not os.path.exists('plots'):
        os.mkdir("plots")
    
    #establishing range of values for heatmap colour scale
    pae_min=1
    pae_max=1
    for file in os.listdir('pae/'):
        if file.endswith(".txt"):
            pae = np.loadtxt(f"pae/{file}")
            if pae.min() < pae_min:
                pae_min = pae.min()
            if pae.max() > pae_max:
                pae_max = pae.max()
    
    df["Mean pLLDT"] = None
    
    #looping over all genes that underwent structure prediction
    for file in os.listdir('structures/'):
        if file.endswith(".pdb"):
            gene_code = str(file[:4])
        
            with open(f"structures/{file}", "r") as f:
                res_index = []
                all_plddts = []
                for line in f:
                    line = line.rstrip().split()
                    if line[0] != "ATOM":
                        continue
                    all_plddts.append(float(line[-2]))
                    res_index.append(int(line[5]))
                
            mean_plddt = sum(all_plddts) / len(all_plddts)
            df.loc[df['ID'] == gene_code, 'Mean pLLDT'] = mean_plddt
        
            temp_df = pd.DataFrame(list(zip(res_index, all_plddts)),columns =['Res', 'pLDDT'])
            temp_df = temp_df.groupby("Res")["pLDDT"].mean()        
            plddt = temp_df.to_numpy()
        
            pae = np.loadtxt(f"pae/{gene_code}_pae.txt")

            make_confidence_plots(output_dir='plots',gene_code=gene_code,plddt=plddt,pae=pae,vmin=pae_min,vmax=pae_max,scale=scale)
            
def plot_tm(df, filename, threshold=0.5):
    """
    Generates a 2 by 2 grid of box plots to visualise TMAlign output.  Top row = all TMScores in df dataframe.
    Bottom row = TMScores in df dataframe > threshold.  Plots on right are grouped by template.  Saves plot as a .png file.

    Parameters:
        df (dataframe)       : Dataframe with columns 'TMScore' and 'Template'.
        threshold (float)    : Float between 0 and 1.
        filename (string)    : Filename to save plot under.
    """
    filt=df[df['TMScore'] > threshold].reset_index(drop=True)
    
    fig, axs = plt.subplots(2,2,figsize=(7,7),gridspec_kw={'width_ratios': [1, 4]})
    sns.set_palette('viridis')
    sns.boxplot(data=df['TMScore'],ax=axs[0,0])
    axs[0,0].set(ylabel='TMScore',xticklabels=['All'])
    axs[0,0].set_title('All scores',fontsize=8,loc='left')
    sns.boxplot(data=filt['TMScore'],ax=axs[1,0])
    axs[1,0].set(ylabel='TMScore',xticklabels=['All'])
    axs[1,0].set_title(f'Scores over {threshold} threshold',fontsize=8,loc='left')
    sns.set_palette('Blues')
    sns.boxplot(y=df['TMScore'],x=df['Template'],ax=axs[0,1])
    axs[0,1].set(xlabel=None,ylabel=None)
    sns.boxplot(y=filt['TMScore'],x=filt['Template'],ax=axs[1,1])
    axs[1,1].set(xlabel=None,ylabel=None)
    plt.savefig(f'monomer_prediction/plots/{filename}.png')
    plt.close()
    
def plot_template_comp(df,filename,templates=False):
    """
    Generates heatmap to visualise TM-scores of templates compared to themselves.  Saves as a .png file.

    Parameters:
        filename (string)    : Absolute or relative path to .csv file containing TM Scores.
    """
    if templates == False:
        templates = list(df['Template'].unique())
    
    template_matrix = []
    
    for query in templates:
        data = df[df["Query"] == query]
        query_vector = []
        for template in templates:
            query_vector.append(float(data[data["Template"]==template]["TMScore"]))
        template_matrix.append(query_vector)
    
    df=pd.DataFrame(template_matrix)
    df.columns = templates
    df.index = templates

    sns.heatmap(df, vmin=0, vmax=1, cbar=True, cmap="Blues")
    plt.savefig(f'monomer_prediction/plots/{filename}.png')
    plt.close()