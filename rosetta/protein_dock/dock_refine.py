#rosetta_local_docking
#删除ZN
#sed -i -e '/ZN/d' *.pdb
#rosetta 准备pdb文件

score_jd2.static.linuxgccrelease -in:file:l list_file  -out:pdb  -ignore_zero_occupancy false -overwrite
#对接前准备
nohup mpirun -np 20 docking_prepack_protocol.mpi.linuxgccrelease -s *.pdb -docking:partners A_B -docking::dock_rtmin true -docking::sc_min true -out:pdb -overwrite &

#local dock
import os
direct=''
dir_new=[direct+'/'+i for i in os.listdir(direct)]
for i in os.listdir(direct):
    os.chdir(direct)
    os.mkdir(i.split('.')[0])
    shutil.copy(i,i.split('.')[0])
    os.chdir(i.split('.')[0])
    os.system('nohup mpiexec -np  20 docking_protocol.mpi.linuxgccrelease -s '+i+' -partners A_B -dock_pert 3 8 -ex1 -ex2aro  -use_input_sc -dock_mcm_trans_magnitude 0.7 -dock_mcm_rot_magnitude 5.0 -nstruct 1000')
#dock refine
#find 
import shutil
df=pd.read_csv('score.sc',sep='\s+',skiprows=[0])
df.sort_values(by='total_score',inplace=True)
shutil.copy(list(df['description'])[0]+'.pdb','../')

#local dock refine
import os,shutil
direct=''
dir_new=[direct+'/'+i for i in os.listdir(direct)]
for i in os.listdir(direct):
    os.chdir(direct)
    os.mkdir('../dock_refine/'+i.split('.')[0])
    shutil.copy(i,'../dock_refine/'+i.split('.')[0])
    os.chdir('../dock_refine/'+i.split('.')[0])
    os.system('nohup mpiexec -np  20 docking_protocol.mpi.linuxgccrelease -s '+i+' -partners B_A  -ex1 -ex2aro -use_input_sc -dock_mcm_trans_magnitude 0.7 -dock_mcm_rot_magnitude 5.0 -nstruct 1000 -docking_local_refine')
 



