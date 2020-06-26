import os
from ala_scan import *
from pickout_hotloop import *
from os import path
import shutil 
import sys

def main():
    pdb_files = os.listdir("../")
    for pdb_filename in pdb_files:
        if pdb_filename == "test.pdb" or pdb_filename[-3:] != "pdb":
            continue
        print("Processing %s ..." % pdb_filename)
        partners = "A_B"
        mutant_aa = "A"
        neighbor_cutoff = 8.0
        interface_cutoff = 4.0
        repack_cutoff = 6.5
        trials = 20
        trial_output = "ddG_out"
        output = False
        pdb_filename = "../" + pdb_filename
        scanning(pdb_filename, partners, mutant_aa,
            neighbor_cutoff, interface_cutoff, repack_cutoff,
            trials, trial_output, output)

if __name__ == '__main__':
    # main()
    if path.exists("results_filtered.txt"):
        os.remove("results_filtered.txt")
    hotloop_found = False
    for directory in extract_files():
        file_directory_head = "Results_Clean"
        hotspots = extract_hotspots(directory)
        if identify_hotloop(hotspots):
            pdb_filename = "../" + directory.split("/")[1][0:-4] + ".pdb"
            assert os.path.exists(pdb_filename)
            file_directory = "/".join((file_directory_head, directory.split("/")[1][0:-4]))
            if not os.path.exists(file_directory):
                os.makedirs(file_directory) 
            pdb_copy = "/".join((file_directory, pdb_filename.split("/")[1]))
            shutil.copyfile(pdb_filename, pdb_copy)
            file_name = file_directory + "/out_hotspots.txt"
            hotloop_found = True
            write = False
            result_file = directory + "/ddG_out_total.txt"
            with open(result_file, "r") as fi:
                lines = fi.readlines()
                with open("results_filtered.txt", "a") as fo:
                    with open(file_name, "w+") as fo_separate:
                        header = "=> " + result_file + " <=\n"
                        fo.write(header)
                        for line in lines:
                            fo.write(line)
                            if write:
                                fo_separate.write(line)
                            elif line == "Likely Hotspot Residues\n":
                                write = True
                        fo.write("\n\n")
    if not hotloop_found:
        print("No Potential Hotloop Found")