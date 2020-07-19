import os
import sys
import shutil
import re

#这个函数是计算sdf文件的总分子数目
def countmol(sdf_file_name):
        sdf_file=open(sdf_file_name,'r')
        sdf_str=sdf_file.read()
        mol_num=sdf_str.count('$$$$') 
        return mol_num;
        
#函数计算分子数目分割后，返回的每个part的分子运算量
def sim_mol(mol_num,core):
        #计算除数(a)，余数(b)
        a=int(mol_num/core)
        b=int(mol_num%core)
        #生成列表，分割后的分子数目
        #若有余数
        if b != 0 :
            mol_list=[[(a*x)+1,a*(x+1)] for x in range(core)]           
            mol_list.append([(a*core)+1,mol_num])
        #若无
        else:
            mol_list=[[(a*x)+1,a*(x+1)] for x in range(core)]
        return mol_list
#可以自定义初始分子与终端分子
def sim_mol_diy(core,fir,last):
        #计算除数(a)，余数(b)
        a=int((last-fir)/core)
        b=int((last-fir)%core)
        #生成列表，分割后的分子数目
        #若有余数
        if b != 0 :
              test_mol_list=[[fir+(x*a)+1,fir+a*(x+1)] for x in range(core)]
              #替换第一个元素，最末尾追加一个元素
              test_mol_list[0]=[fir,fir+a]
              test_mol_list.append([last-b+1,last])
        #若无
        else:
              test_mol_list=[[fir+(x*a)+1,fir+a*(x+1)] for x in range(core)]
              #替换第一个元素以及最末尾的元素
              test_mol_list[0]=[fir,fir+a]
              test_mol_list[-1]=[last-a+1,last]
        return test_mol_list
        
#创建conf文件，子文件夹，修改相关的参数
#并获取创建的文件的相关列表
def make_conf(conf,mol_list):
        path=os.getcwd()
        list_path=[]
        conf_file=open(conf,'r')
        conf_str=conf_file.read()
        for i in range(len(mol_list)):
            #创建临时conf文件,以及其相关文件夹
            num=i+1
            tmp_conf_dir='GOLD_'+str(num)
            tmp_conf_file='GOLD_'+str(num)+'.conf'
            os.mkdir(tmp_conf_dir)
            tmp_conf=open(tmp_conf_file,'w+')
            #创建相关列表[tmp_conf_dir,tmp_conf_file]
            list_path.append([path +'/' + tmp_conf_dir , path + '/' + tmp_conf_dir + '/gold.conf'])
            #替换字符串，替换字符串为:'start_at_ligand 1\n'
            f_str=conf_str
            replace_str='start_at_ligand '+ str(mol_list[i][0]) + ' finish_at_ligand ' + str(mol_list[i][1])
            old_str='start_at_ligand 1'
            f_str=f_str.replace(old_str,replace_str)
            replace_str2='/data/home/pengruchao/sujiaqi/DOCK/GOLD/NS3/' + tmp_conf_dir
            old_str2='directory_target'
            f_str=f_str.replace(old_str2,replace_str2)
            #写入临时文件，移动临时文件
            tmp_conf.write(f_str)
            tmp_conf.close()
            shutil.move(tmp_conf_file,tmp_conf_dir+'/gold.conf')
        conf_file.close()
        return list_path
            
def pub(template_pbs,list_path):
        #创建一个文件夹，存放qsub文件
        path=os.getcwd()
        qsub_path=path+'/qsub_dir'
        os.mkdir(qsub_path)
        #创建新的templat_pbs文件，并且进行移动和复制
        pbs_file=open(template_pbs,'r')
        pbs_str=pbs_file.read()
        for i in range(len(list_path)):
                #创建临时pbs以及相关的out，err文件
                num=i+1
                tmp_pbs_name='GOLD_'+str(i+1)+'.pbs'
                tmp_out_name=qsub_path+'/GOLD_'+str(num)+'.out'
                tmp_err_name=qsub_path+'/GOLD_'+str(num)+'.err'
                tmp_pbs=open(tmp_pbs_name,'w+')
                #修改PBS文本内容，存储
                f_str=pbs_str
                f_str=f_str.replace('#PBS -N gold_test','#PBS -N '+'GOLD'+'_'+str(num))
                f_str=f_str.replace('#PBS -o out','#PBS -o '+tmp_out_name)
                f_str=f_str.replace('#PBS -e err','#PBS -e '+tmp_err_name)
                f_str=f_str.replace('gold.conf',str(list_path[i][1]))
                tmp_pbs.write(f_str)
                tmp_pbs.close()
                shutil.move(tmp_pbs_name,qsub_path)
        pbs_file.close()
        return 
        
def qsub(qsub_path):
        list_qsub=os.listdir(qsub_path)
        for i in range(len(list_qsub)):
                gold_pbs_path=qsub_path+'/'+list_qsub[i]
                os.system('nohup qsub '+gold_pbs_path+' &')
        return


if __name__ == '__main__':
        x=countmol('last.sdf')
        mol_list=sim_mol(x,2)
        mol_list=sim_mol(20,10,500)
        list_path=make_conf('gold.conf',mol_list)
        pub('gold.pbs',list_path)
        qsub('/data/home/pengruchao/sujiaqi/DOCK/GOLD/NS3/qsub_dir')
