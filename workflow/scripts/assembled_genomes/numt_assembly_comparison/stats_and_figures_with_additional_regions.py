import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path
import sys
from scipy.stats import pearsonr, spearmanr
sys.path.append(snakemake.params.functions_path)
import functions as custom_functions
#import workflow.scripts.functions as custom_functions

#df_paths = snakemake.input.numt_paths

#genomes_info = pd.read_table("config/genomes.tsv")
genomes_info = pd.read_table(snakemake.input.genome_info)

dfs=[]
dropped={}

for index,row in genomes_info.iterrows():
    if row['Species'] in ['Bos_taurus','Capra_hircus']:
        tempdf=pd.read_csv(f"results/numt_discovery/numts/{row['Species']}/{row['Assembly Name']}/{row['Assembly Accession']}_{row['Assembly Name']}_numts.csv")
    else:
        tempdf=pd.read_csv(f"data/genome_comparison/additional_regions_compatible_seqs/{row['Species']}/{row['Assembly Name']}/{row['Assembly Accession']}/{row['Genome Abbreviation']}/all_numts/{row['Genome Abbreviation']}.csv")
    tempdf['genome']=row['Genome Abbreviation']
    dfs.append(tempdf)

df=pd.concat(dfs)
df=df.reset_index()

table1 = pd.DataFrame(columns=['total_numts','assembly_numts','total_numt_size',
    'genome_size','genome_coverage','genome_n50','genome_id','smallest_numt', 'biggest_numt', 'mean_size','size_std',
    'median_size', 'size_mad', 'lowest_identity', 'highest_identity', 'mean_identity', 'identity_stdev',
    'regions_with_multiple_fragments', 'regions_with_single_fragments','avg_multifragment_region_size', 'total_regions'])

for name,group in df.groupby(by='genome'):
    assembly = group
    assembly_numts = len(assembly)
    total_numts = len(group)
    genome_size = genomes_info[genomes_info['Genome Abbreviation']==name]['Total Sequence Length'].iloc[0]
    total_numt_size = group['nuc_size'].sum()#+group['ref_nuc_size'].sum()
    genome_coverage_percent = round((total_numt_size/genome_size)*100, 3)
    genome_id = genomes_info[genomes_info['Genome Abbreviation']==name]['Assembly Accession'].iloc[0]
    genome_n50 = genomes_info[genomes_info['Genome Abbreviation']==name]['Contig N50'].iloc[0]
    smallest_numt = group['nuc_size'].min()
    biggest_numt = group['nuc_size'].max()
    mean_size = group['nuc_size'].mean()
    size_stdev = group['nuc_size'].std()
    median_size = group['nuc_size'].median()
    size_mad = group['nuc_size'].mad()
    lowest_identity = group['%_identity'].min()
    highest_identity = group['%_identity'].max()
    mean_identity = group['%_identity'].mean()
    identity_stdev = group['%_identity'].std()
    fragmented_regions = len([g for _, g in group.groupby("region_id") if len(g) > 1])
    single_insertion_regions = len([g for _, g in group.groupby("region_id") if len(g) == 1])
    multifragment_regions_sizes = [g['region_length'].iloc[0] for i,g in group.groupby("region_id") if len(g) >1]
    avg_multifragment_region_size =np.rint(np.mean(multifragment_regions_sizes))
    total_regions=len(group['region_id'].unique())
    table1.loc[name]=[int(total_numts), int(assembly_numts),
        #int(ngs_numts),
        int(total_numt_size),
        int(genome_size), genome_coverage_percent, genome_n50, genome_id, smallest_numt, biggest_numt, mean_size, size_stdev,
        median_size, size_mad, lowest_identity, highest_identity, mean_identity, identity_stdev,
        fragmented_regions, single_insertion_regions, int(avg_multifragment_region_size), total_regions]

table1=table1.loc[genomes_info['Genome Abbreviation'].tolist()]


table1.to_csv(snakemake.output[0])

Path(snakemake.params[0]).mkdir(parents=True, exist_ok=True)
Path(snakemake.params[1]).mkdir(parents=True, exist_ok=True)


size_ranges = list(range(2000))[0::100]
size_ranges+=list(range(2000, 6000))[0::1000]
size_ranges.append(15000)

size_ranges_larger = [0,100,250,500,15000]
identity_ranges = [50,60,70,80,90,100]

for name, group in df.groupby(by='genome'):
    assembly_name=genomes_info[genomes_info['Genome Abbreviation']==name]['Assembly Name'].iloc[0]
    Path(f"{snakemake.params[0]}/{name}").mkdir(parents=True, exist_ok=True)
    assembly=group#[group['data_origin']=='assembly']
    custom_functions.assign_range_and_plot(assembly, range=size_ranges, outcol=f'NUMT size range in {assembly_name} genome (bp)', y_axis_label='Count', outpath=f"{snakemake.params[0]}/{name}/numt_size_distribution_assembly.png")
    custom_functions.assign_range_and_plot(assembly, range=identity_ranges, outcol=f'NUMT sequence identity range in {assembly_name} genome', outpath=f"{snakemake.params[0]}/{name}/numt_identity_distribution_assembly.png")

for name, group in df.groupby(by='genome'):
    assembly=group#[group['data_origin']=='assembly']
    assembly=custom_functions.assign_range_and_plot(assembly, range=size_ranges_larger, outcol='nuc_size_range_larger')
    assembly['nuc_size_range_larger']=assembly['nuc_size_range_larger'].astype(str)+" (bp)"
    size_ranges_order = assembly['nuc_size_range_larger'].sort_values().unique()
    sns.set(rc={'figure.figsize':(20,10)})
    g = sns.FacetGrid(assembly, col="nuc_size_range_larger",col_order=size_ranges_order.tolist())
    g.map(sns.kdeplot, "%_identity")
    g.set_titles('{col_name}')
    g.set_axis_labels("Identity (%)", "Density")
    Path(f"{snakemake.params[1]}/{name}").mkdir(parents=True, exist_ok=True)
    g.savefig(f"{snakemake.params[1]}/{name}/identity_distribution_by_size_range.png", dpi=300, bbox_inches = "tight")

plt.clf()
sns.set(rc={'figure.figsize':(20,15)})
size_kde=sns.kdeplot(data=df, x='nuc_size', hue='genome', legend=False)
size_kde.set(xlabel='NUMT size')
#size_kde.legend.remove()
plt.savefig(f"{snakemake.params[1]}/size_distribution.png", dpi=300, bbox_inches = "tight")

plt.clf()
sns.set(rc={'figure.figsize':(20,15)})
id_kde = sns.kdeplot(data=df, x='%_identity', hue='genome', legend=False)
id_kde.set(xlabel='NUMT sequence identity (%)')
#id_kde.legend.remove()
plt.savefig(f"{snakemake.params[1]}/identity_distribution.png", dpi=300, bbox_inches = "tight")
