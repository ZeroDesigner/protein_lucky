#multiprocessing to run gromacs
import os,shutil
import multiprocessing
#path to protein
original_dir='./protein'
original_file=[original_dir+'/'+i for i in os.listdir(original_dir)]
#path to 
direct_path='./md'
direct_dir_list=[direct_path+'/'+i.split('.')[0] for i in os.listdir(original_dir)]
#path to mdp file
eg_mdp_dir='./eg_mdp'
for i in range(len(direct_dir_list)):
    os.mkdir(direct_dir_list[i])
    shutil.copy(original_file[i],direct_dir_list[i]+'/'+'original.pdb')
    for x in os.listdir(eg_mdp_dir):
        shutil.copy(eg_mdp_dir+'/'+x,direct_dir_list[i])
#add string as you need
str1='gmx pdb2gmx -ignh -ff amber99sb-ildn -f original.pdb -o original.gro -p topol.top -water tip3p'
str2='gmx editconf -f original.gro -o original-PBC.gro -c -d 1.0 -bt cubic'
str3='gmx solvate -cp original-PBC.gro -cs spc216.gro -o original_solv.gro -p topol.top'
str4='gmx grompp -f ions.mdp -c original_solv.gro -p topol.top -o ions.tpr'

str5='echo "SOL"|gmx genion -s ions.tpr -o original_solv.gro -p topol.top -pname NA -nname CL -neutral'
str6='gmx grompp -f minim.mdp -c original_solv.gro -p topol.top -o em.tpr'

str7='gmx mdrun -v -deffnm em'

str8='gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr'
str9='gmx mdrun -deffnm nvt -v'

str10='gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr'
str11='gmx mdrun -deffnm npt -v'

str12='gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr'
str13='gmx mdrun -deffnm md_0_1'
#simply run 
def gmx_gro(i):
    os.chdir(i)
    os.system(str1)
    os.system(str2)
    os.system(str3)
    os.system(str4)
    os.system(str5)
    os.system(str6)
    os.system(str7)
    os.system(str8)
    os.system(str9)
    os.system(str10)
    os.system(str11)
    os.system(str12)
    os.system(str13)
pool = multiprocessing.Pool(processes = 1)
for i in direct_dir_list:
    pool.apply_async(gmx_gro, (i, ))
pool.close()
pool.join()
