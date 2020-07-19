import pandas as pd
import os
import shutil
abs_path='direct/run_cm/'
list_dir=os.listdir(os.getcwd())
list_dir=[abs_path+i for i in list_dir]
score_pdb=[]
for i in list_dir:
    df=pd.read_csv(i+'/score.sc',sep='\s+',skiprows=[0])
    df.sort_values(by='total_score',ascending=False,inplace=True)
    score_pdb.append([i+'/'+df['description'][df['total_score'].idxmax()]+'.pdb',df['total_score'][df['total_score'].idxmax()]])
df2=pd.DataFrame(score_pdb)
df2.columns=['location','score']
df2['location'][df2['score'].idxmax()]
df2.to_csv('sum_cm.csv',index=None)
tmp_list=df2['location'][df2['score'].idxmax()].split('/')
name=tmp_list[-2]+'_'+tmp_list[-1]
shutil.copy(df2['location'][df2['score'].idxmax()],abs_path+name)
