purpose:use rosettta-peptidriver to find a peptide which contribut lots to delta G
step:
(1) add the protein complex structure to 'eg' dir, just two chain (A,B) , PDB format, A is the receptor ,B is the partner
(2) change the num in 'pep_list=list(np.linspace(3,111,109)[11:])' in the file-'pep.py'
3 is the minimal lenth of the peptide
111 is the largest lenth of the peptide
109 = the largest lenth of the peptide - the minimal lenth of the peptide
(3)'pool = multiprocessing.Pool(processes = 20)'
chang the num '20' to fit your own cpu
(4) run pep.py
(5) when you  finsihed , change the 'input_protein.peptiderive.txt' in the str -
input_file_list=[input+'/'+i+'/'+'input_protein.peptiderive.txt' for  i in  os.listdir(input) if  os.path.isdir(i) ]
in the file 'summary_pep_result.py'
(6)run summary_pep_result.py