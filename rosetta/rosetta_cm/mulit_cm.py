import os,shutil
import multiprocessing
direct_dir='direct_dir/run_cm'
direct_dir_list=[]
for i in range(15):
    os.mkdir(direct_dir+'/'+str(i))
    direct_dir_list.append(direct_dir+'/'+str(i))
    shutil.copy('direct_dir/flags',direct_dir+'/'+str(i))

def rosetta_script(i):
    os.chdir(i) 
    os.system('rosetta_scripts.static.linuxgccrelease @flags')
    return 


pool = multiprocessing.Pool(processes = 15)
for i in direct_dir_list:
    pool.apply_async(rosetta_script, (i, )) 
pool.close()
pool.join() 
