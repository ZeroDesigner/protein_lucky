import os
import pandas as pd
input='./'
input_file_list=[input+'/'+i+'/'+'input_protein.peptiderive.txt' for  i in  os.listdir(input) if  os.path.isdir(i) ]
information=[]
file_colounm=['Receptor', 'Partner', 'Peptide length', 'Position', 'Interface score', 'Relative interface score (%)', 'Sequence']
final_info_list=[]
for i in input_file_list:
	f=open(i,'r')
	f_lines=f.readlines()
	for x in range(len(f_lines)):
		if f_lines[x] == '## Best linear peptides for all chain pairs\n':
			tmp=f_lines[x+4].split('|')[1:-1]
			new=[i.strip() for i in tmp]
	final_info_list.append(new)
	f.close()
df=pd.DataFrame(final_info_list)
df.columns=file_colounm
df[["Peptide length"]] = df[["Peptide length"]].astype(int)
df.sort_values(by=["Peptide length"],inplace=True)
df.to_csv('summary_peptidriver.csv',index=None)
