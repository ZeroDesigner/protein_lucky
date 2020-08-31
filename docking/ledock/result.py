import os
path='/Users/sujiaqi/Desktop/2bsm'
f_result=open('result.txt','w+')
f_result.write('ligand score(kcal/mol)\n')
for i in os.listdir(path):
    if 'dok' in i :
        f=open(i,'r')
        f_list=f.readlines()
        f_result.write(i.split('.')[0]+' '+f_list[1].split()[-2]+'\n')
f_result.close()

