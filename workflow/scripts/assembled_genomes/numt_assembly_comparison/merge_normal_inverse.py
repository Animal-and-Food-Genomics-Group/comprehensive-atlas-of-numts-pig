import pandas as pd
from pybedtools import BedTool

main_numts_df=pd.read_csv(snakemake.input.ref_numts)
other_numts_df=pd.read_csv(snakemake.input.other_numts)


#other_species="Sus_scrofa"
#other_assembly_name='Sscrofa10.2'
#other_assembly_accession='GCF_000003025.5'
#other_species="Sus_cebifrons"
#other_assembly_name='Sus_cebifrons.v1'
#other_assembly_accession='GCA_905335845.1'
#other_species='Sus_scrofa'
#other_assembly_name='Berkshire_pig_v1'
#other_assembly_accession='GCA_001700575.1'
#
#main_numts_df=pd.read_csv("results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv")
#other_numts_df=pd.read_csv(f"results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv")

#other_numts_df=pd.read_csv(f"results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv")


#a=pd.read_csv(f"test/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.merged.csv")
#b=pd.read_csv(f"test/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.merged.csv")
#c=pd.read_csv(f"test/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.merged.csv")
#d=pd.read_csv(f"test/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.merged.csv")

a=pd.read_csv(snakemake.input.a)
b=pd.read_csv(snakemake.input.b)
c=pd.read_csv(snakemake.input.c)
d=pd.read_csv(snakemake.input.d)


a=pd.concat([a,b])

b=pd.concat([c,d]).drop_duplicates()

b.columns=[i.replace("other","TEMP") for i in b.columns]
b.columns=[i.replace("main","other") for i in b.columns]
b.columns=[i.replace("TEMP","main") for i in b.columns]

c=pd.concat([a,b]).drop_duplicates()

non_corr=c[(c['main_region_id']=='.')|(c['other_region_id']=='.')]
corr=c[(c['main_region_id']!='.')&(c['other_region_id']!='.')]

corr=corr[corr['status']!='DROP']

non_corr=non_corr[non_corr['status']!='DROP']
non_corr=non_corr.reset_index(drop=True)

private=non_corr[non_corr['status']=='PRIVATE'].sort_values(by='alignment_length', ascending=False)

private_main=private[private['main_region_id']!='.']
private_other=private[private['other_region_id']!='.']

private_other=private_other.drop_duplicates(subset='other_region_id',keep='first')
private_main=private_main.drop_duplicates(subset='main_region_id',keep='first')

corr=corr.sort_values(by=['alignment_length','main_chromosome','other_chromosome'], ascending=False)

corr=corr.drop_duplicates(subset=['main_region_id'], keep='first')
corr=corr.drop_duplicates(subset=['other_region_id'], keep='first')

comp=non_corr[non_corr['status']=='COMPATIBLE']

comp=comp.sort_values(by='alignment_length', ascending=False)
comp_main=comp[comp['main_region_id']!='.']
comp_other=comp[comp['other_region_id']!='.']

comp_main=comp_main.drop_duplicates(subset=['main_region_id'], keep='first')
comp_other=comp_other.drop_duplicates(subset=['other_region_id'], keep='first')


outdf=pd.DataFrame(columns=['main_region_id','main_chromosome','main_region_start','main_region_end','main_region_status','other_region_id','other_chromosome','other_region_start','other_region_end','other_region_status'])

for index,row in corr.iterrows():
    tempmain=main_numts_df[main_numts_df['region_id']==row['main_region_id']].iloc[0]
    tempother=other_numts_df[other_numts_df['region_id']==row['other_region_id']].iloc[0]
    outdf.loc[len(outdf)]=[tempmain['region_id'],tempmain['chromosome_code'], tempmain['region_start'], tempmain['region_end'],'CORRESPONDING',
                            tempother['region_id'],tempother['chromosome_code'], tempother['region_start'], tempother['region_end'],'CORRESPONDING',]


for index,row in comp_main.iterrows():
    tempmain=main_numts_df[main_numts_df['region_id']==row['main_region_id']].iloc[0]
    outdf.loc[len(outdf)]=[tempmain['region_id'],tempmain['chromosome_code'], tempmain['region_start'], tempmain['region_end'],'CORRESPONDING',
                            '.',row['other_chromosome'],row['other_region_start'],row['other_region_end'],'COMPATIBLE_SEQUENCE']


for index,row in comp_other.iterrows():
    tempother=other_numts_df[other_numts_df['region_id']==row['other_region_id']].iloc[0]
    outdf.loc[len(outdf)]=['.',row['main_chromosome'],row['main_region_start'],row['main_region_end'],'COMPATIBLE_SEQUENCE',
                            tempother['region_id'],tempother['chromosome_code'], tempother['region_start'], tempother['region_end'],'CORRESPONDING']


for index,row in private_main.iterrows():
    tempmain=main_numts_df[main_numts_df['region_id']==row['main_region_id']].iloc[0]
    outdf.loc[len(outdf)]=[tempmain['region_id'],tempmain['chromosome_code'], tempmain['region_start'], tempmain['region_end'],'CORRESPONDING',
                            '.',row['other_chromosome'],row['other_region_start'],row['other_region_end'],'PRIVATE']


for index,row in private_other.iterrows():
    tempother=other_numts_df[other_numts_df['region_id']==row['other_region_id']].iloc[0]
    outdf.loc[len(outdf)]=['.',row['main_chromosome'],row['main_region_start'],row['main_region_end'],'PRIVATE',
                            tempother['region_id'],tempother['chromosome_code'], tempother['region_start'], tempother['region_end'],'CORRESPONDING']

outdf=outdf.sort_values(by='main_region_id')

outdf.to_csv(snakemake.output[0], index=False)
