#not test
import os
import shutil
import multiprocessing as mp
fasta_dir=''
psipred_path=''
fasta_list=os.listdir(fasta_dir)
def run_psipred(i):
    os.chdir(fasta_dir)
    os.mkdir(i+'_psipred')
    shutil.copy(i,i+'_psipred')
    os.chdir(fasta_dir+'/'+i+'_psipred')
    os.system('nohup '+ psipred_path+' '+i)
pool = multiprocessing.Pool(processes = 20)
for i in fasta_list:
    pool.apply_async(run_psipred, (i, )) 
pool.close()
pool.join() 
