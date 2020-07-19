#--coding:utf-8--
import sys
import math
from collections import defaultdict
import pandas as pd
if __name__ == "__main__":

    #IO,文件转列表
    f=open(infile,'r')
    content = [l.strip().split(',') for l in f]
	
    #产生字典，键mur：值score
    del content[0]
    dic = defaultdict(list)
    for line in content:
        key = line[3]
        dic[key].append(line[0])
        
    #临时列表，tmp_list=[[Mur,max_score]]
    tmp_list = []
    for key,value in dic.items():
        tmp_value=[float(i) for i in value]
        max_value= max(tmp_value)
        tmp_list.append([key,max_value])
    
    #终列表，result_list=[[Score,ligand,SMI,Mur]]
    result_list=[]
    for line in content:
        for i in range(len(tmp_list)):
            if (float(line[0])==tmp_list[i][1] and line[3]==tmp_list[i][0]):
                result_list.append(line)
    #导入数据库，赋值表头，排序，输出
    data=pd.DataFrame(result_list)
    data.columns=['Score','ligand','SMI','Mur']
    data.sort_values('Score',ascending=False) 
    data.to_csv('fin_score_lig_SMI_Mur.csv',index=False)
