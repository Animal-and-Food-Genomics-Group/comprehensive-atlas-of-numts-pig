import pandas as pd
import sys
sys.path.append(snakemake.params.functions_path)
import functions as custom_functions

df = custom_functions.merge_maf_and_blasttab(snakemake.input[0], snakemake.input[1])
linear_df = custom_functions.merge_maf_and_blasttab(snakemake.input[2], snakemake.input[3])

consensus_length=linear_df['mtc_srcsize'].iloc[0]

df_before=df[df['mtc_start']<consensus_length]
df_over = df_before[df_before['mtc_end']>consensus_length]

samestart=linear_df[linear_df['mtc_start'].isin(df_over['mtc_start'])]
shorterend=samestart[samestart['mtc_end']<=consensus_length]

over_sets = {}
over_chroms = []
for index,row in df_over.iterrows():
    over_sets[index] = set(range(row['nuc_start'],row['nuc_end']+1))
    over_chroms.append(row['chromosome'])

linear_sets = {}
for index,row in linear_df.iterrows():
    if row['chromosome'] in over_chroms:
        linear_sets[index] = set(range(row['nuc_start'],row['nuc_end']+1))

todrop=set()
toadd=set()
for i,j in over_sets.items():
    for k,v in linear_sets.items():
        if df_over.loc[i]['chromosome']==linear_df.loc[k]['chromosome']:
            if v.issubset(j):
                todrop.add(k)
                toadd.add(i)

rows_to_drop=linear_df.loc[list(todrop)]
rows_to_drop.to_csv(snakemake.output[0])

rows_to_add=df.loc[list(toadd)]
rows_to_add.to_csv(snakemake.output[1])

linear_df=linear_df.drop(list(todrop))

outdf = pd.concat([linear_df, rows_to_add])

outdf.to_csv(snakemake.output[2], index=False)
