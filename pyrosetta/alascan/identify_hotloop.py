import sys
import numpy as np
import os
def extract_ddG(filename):
	ddG_total = []
	index_total = []
	chain_total = []

	ddG_files = []
	with open(filename) as fi:
		lines = fi.readlines()
	read = False
	ddG_individual = []
	index_individual = []
	chain_individual = []
	for line in lines:
		if line[0:2] == "=>":
			ddG_files.append(line)
		if read:
			if line[0:3] == "===":
				read = False
				ddG_total.append(ddG_individual)
				index_total.append(index_individual)
				chain_total.append(chain_individual)
				ddG_individual = []
				index_individual = []
				chain_individual = []
			else:
				ddG = float(line.split()[3])
				ddG_individual.append(ddG)
				index = int(line.split()[1])
				index_individual.append(index)
				chain = line.split()[0]
				chain_individual.append(chain)
		else:
			if line[0:3] == "All":
				read = True
			else:
				continue
	return ddG_files, ddG_total, index_total, chain_total

def is_heat_nine(ddG_list, residues, chain_list):
	assert len(ddG_list) == len(residues)
	assert len(residues) == len(chain_list)
	if len(ddG_list) < 3:
		return False
	ddG_sum = np.sum(np.array(ddG_list))
	ddG_cutoff = 0.5 * ddG_sum
	for i in range(len(ddG_list) - 3):
		total_res_above_cutoff = 0
		start_resid = residues[i]
		if (residues[i + 3] - start_resid > 10) or (chain_list[i + 3] != chain_list[i]):
			continue
		j = i + 3
		while j < len(ddG_list) and residues[j] - start_resid < 10 and chain_list[j] == chain_list[i]:
			j += 1
		j = j - 1
		ddG_hotloop = np.sum(np.array(ddG_list[i:j+1]))
		if ddG_hotloop >= ddG_cutoff:
			if ddG_hotloop / (residues[j] - residues[i] + 1) >= 0.6:
				for k in range(i, j + 1):
					if ddG_list[k] >= 1.0:
						total_res_above_cutoff += 1
				if total_res_above_cutoff >= 3:
					return True
	return False



if __name__ == '__main__':
	filename = sys.argv[1]
	ddG_files, ddG_total, index_total, chain_total = extract_ddG(filename)

	assert len(ddG_files) == len(ddG_total)
	assert len(ddG_files) == len(index_total)
	assert len(ddG_files) == len(chain_total)

	for i in range(len(ddG_files)):
		ddG_file = ddG_files[i]
		ddG_list = ddG_total[i]
		residues = index_total[i]
		chain_list = chain_total[i]
		if is_heat_nine(ddG_list, residues, chain_list):
			print(ddG_file)

