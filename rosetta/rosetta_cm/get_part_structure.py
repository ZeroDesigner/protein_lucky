import pandas as pd
import Bio
from Bio.PDB import *
import os
import multiprocessing as mp
pdbl = PDBList() 
df=pd.read_csv('target.csv')
id_list=list(df['Accession  '])
pdb_list=[i.split('_')for i  in id_list]
def get_pdb_part(i):

    os.mkdir('direct_dir/'+i[0])
    os.chdir('direct_dir/'+i[0])
    print(i[0])
    pdbl.retrieve_pdb_file(i[0], pdir = '.', file_format = 'pdb')
    parser = PDBParser(PERMISSIVE=1) 
    pdb_file=i[0].lower()
    pdb_data = parser.get_structure(i[0],'pdb'+pdb_file+'.ent')
    chain_id=i[1]
    print(len(pdb_data[0][chain_id].get_unpacked_list()))
    start_num=pdb_data[0][chain_id].get_list()[0].get_id()[1]
    Dice.extract(pdb_data[0],chain_id,start_num,start_num+len(pdb_data[0][chain_id].get_unpacked_list()),'direct_dir/part_structure/'+i[0]+'_'+i[1]+'.pdb')

for i in pdb_list:
    try:
        get_pdb_part(i)
    except:
        True
pool = mp.Pool(processes = 30)
for i in pdb_list:
    pool.apply_async(get_pdb_part, (i, ))
pool.close()
pool.join()
