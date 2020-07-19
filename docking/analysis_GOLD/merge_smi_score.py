import pandas as pd

data1=pd.read_csv('test1.csv')    
data2=pd.read_csv('test2.csv')  

data1.columns=['Score','ligand'] 
data2.columns=['SMILES','ligand'] 

data3=pd.merge(data1,data2,on='ligand')    
data3.to_csv('Score_Ligand_SMI.csv',index=False)   

data4=data3[['SMILES','ligand']]
data4.to_csv('10000.smi',sep=' ',index=False,header=None)            
