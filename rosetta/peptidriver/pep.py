import os
import shutil
import multiprocessing
import numpy as np
direct_path='./'
original_dir='./eg'
original_file=[original_dir+'/'+i for i in os.listdir(original_dir)]
pep_list=list(np.linspace(3,111,109)[11:])
for i in pep_list:
    new_dir=direct_path+'/'+str(int(i))+'_pep'
    os.mkdir(new_dir)
    for x in original_file:
        shutil.copy(x,new_dir)
    f=open(new_dir+'/new_peptide.flag','r')
    f_str=f.read()
    x=f_str.replace('num',str(int(i)))
    f1=open(new_dir+'/new_peptide.flag','w+')
    f1.write(x)
    f1.close()

new_dir_list=[direct_path+'/'+str(int(i))+'_pep' for i in pep_list]

def mulit_pep(i):
    os.chdir(i)
    os.system('PeptideDeriver.static.linuxgccrelease @new_peptide.flag')

pool = multiprocessing.Pool(processes = 20)
for i in new_dir_list:
    pool.apply_async(mulit_pep, (i, )) 
pool.close()
pool.join() 
