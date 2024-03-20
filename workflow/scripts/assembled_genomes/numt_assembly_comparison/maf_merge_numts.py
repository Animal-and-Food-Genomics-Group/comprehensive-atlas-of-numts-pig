import pandas as pd
from Bio import AlignIO
from pybedtools import BedTool
import numpy as np


#other_species='Sus_scrofa'
#other_assembly_name='Berkshire_pig_v1'
#other_assembly_accession='GCA_001700575.1'

#other_numts_df=pd.read_csv(f"results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv")
#main_numts_df=pd.read_csv("results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv")

main_numts_df=pd.read_csv(snakemake.input.main_numts)
other_numts_df=pd.read_csv(snakemake.input.other_numts)

def maf_to_csv(maf_file, score_threshold=None, csv_outpath=None):
    #this function takes as input a .maf file and returns a pandas DataFrame containing all of the file's information
    #maf_file = os.path.abspath(maf_file)
    maf_dict = {'score':[], 'main_chromosome':[], 'other_chromosome':[],'mismap':[],
                'main_start':[], 'main_strand':[], 'main_srcsize':[], 'main_size':[],
                'other_start':[], 'other_strand':[], 'other_srcsize':[], 'other_size':[],
                'main_seq':[], 'other_seq':[]}

    with open(maf_file) as f:
    #parse the maf file to get score, EG2 and E-value
        lines = f.readlines()
        for line in lines:
            if line.startswith("a score="):
                line = line.strip().split()
                line = line[1:]
                maf_dict['score'].append(line[0].split("=")[1])
                maf_dict['mismap'].append(line[1].split("=")[1])


    for multiple_alignment in AlignIO.parse(maf_file, 'maf'):
        maf_dict['main_chromosome'].append(multiple_alignment[0].id)
        maf_dict['other_chromosome'].append(multiple_alignment[1].id)

        maf_dict['main_start'].append(multiple_alignment[0].annotations['start'])
        maf_dict['main_strand'].append(multiple_alignment[0].annotations['strand'])
        maf_dict['main_srcsize'].append(multiple_alignment[0].annotations['srcSize'])
        maf_dict['main_size'].append(multiple_alignment[0].annotations['size'])
        maf_dict['main_seq'].append(str(multiple_alignment[0].seq))
        maf_dict['other_start'].append(multiple_alignment[1].annotations['start'])
        maf_dict['other_strand'].append(multiple_alignment[1].annotations['strand'])
        maf_dict['other_srcsize'].append(multiple_alignment[1].annotations['srcSize'])
        maf_dict['other_size'].append(multiple_alignment[1].annotations['size'])
        maf_dict['other_seq'].append(str(multiple_alignment[1].seq))

    maf_dataframe = pd.DataFrame.from_dict(maf_dict)

    if score_threshold is not None:
        maf_dataframe = maf_dataframe.loc[maf_dataframe['score'].astype(int) > int(score_threshold)]

    maf_dataframe['other_end']= maf_dataframe['other_start'] + maf_dataframe['other_size']
    maf_dataframe['main_end']= maf_dataframe['main_start'] + maf_dataframe['main_size']
    maf_dataframe['other_start_forward']=maf_dataframe['other_start'].where(maf_dataframe['other_strand']==1, (maf_dataframe['other_srcsize']-maf_dataframe['other_end'])+1)
    maf_dataframe['other_end_forward']=maf_dataframe['other_end'].where(maf_dataframe['other_strand']==1, (maf_dataframe['other_srcsize']-maf_dataframe['other_start']))
    #if private numt perfectly, other start and end are identical, and in negative strand conversion gives end=start-1 basically
    maf_dataframe['other_start_forward']=maf_dataframe['other_start'].where(maf_dataframe['other_start']==maf_dataframe['other_end'], maf_dataframe['other_start_forward'])
    maf_dataframe['other_end_forward']=maf_dataframe['other_end'].where(maf_dataframe['other_start']==maf_dataframe['other_end'], maf_dataframe['other_end_forward'])

    if csv_outpath is not None:
        maf_dataframe.to_csv(csv_outpath)

    return(maf_dataframe)

#df=maf_to_csv(f"data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf")
#df_blast=pd.read_table(f"data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab", header=None)


df=maf_to_csv(snakemake.input.maf)
df_blast=pd.read_table(snakemake.input.blasttab, header=None)

df_blast.columns=['other_chromosome','main_chromosome','%_identity', 'aln_length','mismatches','gap_opens','other_start','other_end','main_start','main_end']
df_blast=df_blast[['%_identity','aln_length','mismatches','gap_opens']]
df=pd.concat([df, df_blast], axis=1)
df['main_strand']=df['main_strand'].replace(1,"+")
df['main_strand']=df['main_strand'].replace(-1,"-")
df['other_strand']=df['other_strand'].replace(1,"+")
df['other_strand']=df['other_strand'].replace(-1,"-")

df=df.reset_index()

df_main_bed=BedTool.from_dataframe(df[['main_chromosome','main_start','main_end','index']])
df_other_bed=BedTool.from_dataframe(df[['other_chromosome','other_start_forward','other_end_forward','index']])

main_numts=pd.read_table(snakemake.input.main_regions, header=None)
#main_numts=pd.read_table("data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed", header=None)
main_numts_bed=BedTool.from_dataframe(main_numts)
other_numts=pd.read_table(snakemake.input.other_regions, header=None)
#other_numts=pd.read_table(f"data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed", header=None)
other_numts_bed=BedTool.from_dataframe(other_numts)
main_merged=df_main_bed.intersect(main_numts_bed, wao=True).to_dataframe()
other_merged=df_other_bed.intersect(other_numts_bed, wao=True).to_dataframe()

main_merged=main_merged[['name','thickEnd']]
main_merged.columns=['index','main_region_id']
df=df.merge(main_merged, on='index', how='outer')
main_numts=main_numts[[1,2,3]]
main_numts.columns=['main_region_start','main_region_end','main_region_id']

df=df.merge(main_numts, on='main_region_id', how='outer')
df['main_region_size']=df['main_region_end']-df['main_region_start']
other_merged=other_merged[['name','thickEnd']]
other_merged.columns=['index','other_region_id']
other_numts=other_numts[[1,2,3]]
other_numts.columns=['other_region_start','other_region_end','other_region_id']
df=df.merge(other_merged, on='index', how='outer')
df=df.merge(other_numts, on='other_region_id',how='outer')
df['other_region_size']=df['other_region_end']-df['other_region_start']

df=df[~df['other_region_id'].isna()]
df=df[~df['main_region_id'].isna()]

df['main_gaps']=df['main_seq'].apply(lambda x:x.count("-"))
df['other_gaps']=df['other_seq'].apply(lambda x:x.count("-"))
df['other_gaps/aln_length']=df['other_gaps']/df['aln_length']
df['main_gaps/aln_length']=df['main_gaps']/df['aln_length']
df['other_size/main_size']=(df['other_size']/df['main_size'])
df['aln_length/main_region_size']=df['aln_length']/df['main_region_size']
df['aln_length/other_region_size']=df['aln_length']/df['other_region_size']
df['score']=df['score'].astype(int)

df=df.sort_values(by=['main_region_id','score'], ascending=False)
df['main_start']=df['main_start'].astype(int)
df['main_end']=df['main_end'].astype(int)

corr_dict={}

multimatches=[]

for name,group in df.groupby(by='main_region_id'):
    if len(group['other_region_id'].unique())==1:
        corr_dict[name]=group['other_region_id'].iloc[0]
    else:
        multimatches.append(group)

for i in multimatches:
    biggest=None
    biggest_sum=0
    for name,group in i.groupby(by='other_region_id'):
        if group['aln_length'].sum()>biggest_sum:
            biggest=name
            biggest_sum=group['aln_length'].sum()
    corr_dict[i['main_region_id'].iloc[0]]=biggest


newdf=[]

for k,v in corr_dict.items():
    newdf.append(df[(df['main_region_id']==k)&(df['other_region_id']==v)])

newdf=pd.concat(newdf)

newdf=newdf.drop(columns='index').drop_duplicates()


for name,group in newdf.groupby(by='other_region_id'):
    if name!='.':
        if len(group['main_region_id'].unique())>1:
            biggest=None
            biggest_sum=0
            newdf=newdf.drop(group.index)
            for name2,group2 in group.groupby(by='main_region_id'):
                if group2['aln_length'].sum()>biggest_sum:
                    biggest=group2
                    biggest_sum=group2['aln_length'].sum()
            newdf=pd.concat([newdf, biggest])
    else:
        for name2,group2 in group.groupby(by='main_region_id'):
            newdf=newdf.drop(index=group2.index)
            group2=group2.sort_values(by='aln_length', ascending=False)
            group2=group2.drop_duplicates(subset='main_region_id', keep='first')
            newdf=pd.concat([newdf,group2])

newdf['other_start_forward']=newdf['other_start'].where(newdf['other_strand']=='+', newdf['other_srcsize']-newdf['other_end'])
newdf['other_end_forward']=newdf['other_end'].where(newdf['other_strand']=='+', newdf['other_srcsize']-newdf['other_start'])
newdf=newdf.reset_index(drop=True)
newdf=newdf.reset_index()

outdf=pd.DataFrame(columns=['main_region_id','other_region_id','other_chromosome','other_region_start','other_region_end','alignment_length','alignment_length_by_main_region_length','alignment_length_by_other_region_length','gaps_by_alignment_length'])

for name,group in newdf.groupby(by='main_region_id'):
    outdf.loc[len(outdf)]=(name,group['other_region_id'].unique()[0], group['other_chromosome'].unique()[0],
    group['other_start_forward'].astype(int).min(), group['other_end_forward'].astype(int).max(),
    group['aln_length'].astype(int).sum(), group['aln_length/main_region_size'].sum().round(2),group['aln_length/other_region_size'].sum().round(2), group['other_gaps/aln_length'].sum().round(2)
    )

status=[]
for index,row in outdf.iterrows():
    if (row['gaps_by_alignment_length']>0.25)&(row['other_region_id']=='.'):
        tempstatus='PRIVATE'
    elif (row['other_region_id']=='.')&(row['alignment_length_by_main_region_length']>0.8):
        tempstatus='COMPATIBLE'
    elif (row['other_region_id']=='.')&(row['alignment_length_by_main_region_length']<=0.8):
        tempstatus='DROP'
    elif (row['alignment_length_by_other_region_length']<0.5)|(row['alignment_length_by_main_region_length']<0.5):
        tempstatus='DROP'
    else:
        tempstatus='CORRESPONDING'
    status.append(tempstatus)

outdf['status']=status

outdf.to_csv(snakemake.output[0], index=False)
