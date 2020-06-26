import sys
import numpy as np
from os import path 
import os
import shutil
import pandas as pd
import copy 

min_residue = 3
max_distance = 8

def sort_function(a):
    subchains = a.split("/")[-1].split("_")[0:2]
    for i in range(len(subchains)):
        subchains[i] = int(subchains[i][1:])
    return subchains

def sort_function_chain(a):
    subchains = a.split("_")
    for i in range(len(subchains)):
        subchains[i] = int(subchains[i])
    return subchains

def sort_function_zip(a):
    return sort_function_chain(a[0])

def extract_head_directory():
    start = os.getcwd()
    root = []
    for root_temp in os.listdir(start):
        if "chain" in root_temp:
            root.append(root_temp)
    return root

def extract_files(head_dir):
    start = "ala_scan/PyRosettaResults/"
    start = "/".join([head_dir, start])
    root = []
    cur_dir = os.getcwd()
    for root_temp, dummy, dummy in os.walk(start):
        if root_temp != start:
            root.append("/".join([cur_dir, root_temp]))
    root.sort(key = sort_function)
    return root

def extract_hotspots(result_file):
    hotspots = []
    result_file = result_file + "/ddG_out_total.txt"
    read = False
    with open(result_file) as fo:
        for line in fo:
            if read:
                hotspots.append(line.strip())
            elif line == "Likely Hotspot Residues\n":
                read = True
    return hotspots

def if_less_than_cutoff(resid_array):
    assert len(resid_array) >= min_residue
    if resid_array[-1] - resid_array[0] <= max_distance:
        return True
    for i in range(len(resid_array) - 2):
        if resid_array[i + 2] - resid_array[i] <= max_distance:
            return True
    return False 

def identify_hotloop(hotspots):
    if len(hotspots) < 3:
        return False
    chainA_count = 0
    chainB_count = 0
    resid_A = []
    resid_B = []
    for info in hotspots:
        chain = info.split()[0]
        resid = int(info.split()[1])
        if chain == "A":
            chainA_count += 1
            resid_A.append(resid)
        else:
            chainB_count += 1
            resid_B.append(resid)

    if chainA_count < min_residue and chainB_count < min_residue:
        return False 

    if chainA_count >= min_residue:
        if if_less_than_cutoff(resid_A):
            return True

    elif chainB_count >= min_residue:
        return if_less_than_cutoff(resid_B)

    return False

def assign_heat(hotspots):
    if len(hotspots) == 0:
        return -1
    elif identify_hotloop(hotspots):
        return 1
    return 0


if __name__ == '__main__':
    reverse = sys.argv[1]
    if reverse == "True":
        reverse = True
    else:
        reverse = False
    root = extract_head_directory()
    d_inter_chain = {}
    d_intra_chain = {}
    inter_chain_index = []
    intra_chain_index = []
    for head_dir in root:
        temp_index = []
        temp_heat = []
        for file in extract_files(head_dir):
            chains = file.split("/")[-4].split("_")[-1]
            subchains = file.split("/")[-1].split("_")[0:2]
            for i in range(len(subchains)):
                subchains[i] = subchains[i][1:]
            if len(chains) == 2:
                if reverse:
                    subchains.reverse()
            subchains = "_".join(subchains)
            temp_index.append(subchains)
            temp_heat.append(assign_heat(extract_hotspots(file)))
        temp_heat = [x for _, x in sorted(zip(temp_index, temp_heat), key=sort_function_zip)]
        temp_index.sort(key=sort_function_chain)
        if len(chains) == 1:
            d_intra_chain.update({chains: temp_heat})
            if len(intra_chain_index) == 0:
                intra_chain_index = copy.deepcopy(temp_index)
            else:
                for i in range(len(temp_index)):
                    assert intra_chain_index[i] == temp_index[i]
        else:
            d_inter_chain.update({chains: temp_heat})
            if len(inter_chain_index) == 0:
                inter_chain_index = copy.deepcopy(temp_index)
            else:
                for i in range(len(temp_index)):
                    assert inter_chain_index[i] == temp_index[i]

    df_inter_chain = pd.DataFrame(data = d_inter_chain, index = list(inter_chain_index))
    df_intra_chain = pd.DataFrame(data = d_intra_chain, index = list(intra_chain_index))

    df_inter_chain.to_csv("inter.csv")
    df_intra_chain.to_csv("intra.csv")

