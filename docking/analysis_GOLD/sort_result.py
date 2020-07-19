import pandas as pd
data = pd.read_csv(file_name, sep='\s+',header=None)   
data.sort_values(0,ascending=False) 
data.columns = ['Score','S(PLP)','S(hbond)','cho','metal','clash','tors','intcor','time','file','ligand']  
data.to_csv(out_file,index=False)    
data2=data[['Score','ligand']]
data2.to_csv(out_file2,index=False) 
