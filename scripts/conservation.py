import subprocess
import sys

import pandas as pd 
import polyclonal

# import gzip
# import shutil

sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

# snakemake variables
input_pdb_file = snakemake.input.input_pdb_file
env_chains = snakemake.params.env_chains
TRO11_site_map = snakemake.input.TRO11_site_map
BF520_site_map = snakemake.input.BF520_site_map
output_csv_file_name = snakemake.output.output_csv_file_name
output_pdb_file_name = snakemake.output.output_pdb_file_name
output_file_sub_name = snakemake.params.output_file_sub_name

# read pdb, as suggested here:
# https://gist.github.com/tbrittoborges/929ced78855945f3e296
colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
            (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
            (78, 80)]

names = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
         'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']

pdb = pd.read_fwf(input_pdb_file, names=names, colspecs=colspecs)

TRO11_df = pd.read_csv(TRO11_site_map)
BF520_df = pd.read_csv(BF520_site_map)

merged_df = TRO11_df.merge(BF520_df, on=['reference_site'])

merged_df = merged_df.assign(conserved=lambda x: x['sequential_wt_x']==x['sequential_wt_y'])
merged_df["conserved"] = merged_df["conserved"].astype(int)

merged_df['wildtype'] = merged_df['sequential_wt_x']
merged_df['mutant'] = 'A'
merged_df.to_csv(output_csv_file_name)

chains = env_chains.split(' ')

# make a dataframe with site escape averages for the eptiope for each chain
dfs = []
for chain in chains:
    chain_df = merged_df[['reference_site', 'conserved']].copy()
    chain_df['chain'] = chain
    dfs.append(chain_df)
df = pd.concat(dfs)

# polyclonal.pdb_utils.reassign_b_factor can only currently take int inputs
df = df[~df['reference_site'].str.contains("a|b|c|d|e|f|g")]
df['reference_site'] = df['reference_site'].astype(int)

df['site'] = df['reference_site']

# save the new PDB file
polyclonal.pdb_utils.reassign_b_factor(input_pdbfile=input_pdb_file, 
                                       output_pdbfile=f"{output_file_sub_name}.pdb", 
                                       df=df, 
                                       metric_col='conserved')
    
    # # compress the new PDB file for the docs
    # with open(f"{output_file_sub_name}_epitope_{epitope}.pdb", 'rb') as f_in:
    # with gzip.open(f"{output_file_sub_name}_epitope_{epitope}.pdb.gz", 'wb') as f_out:
    #     shutil.copyfileobj(f_in, f_out)
