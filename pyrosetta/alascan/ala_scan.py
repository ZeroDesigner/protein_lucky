#!usr/bin/env python

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.docking import calc_interaction_energy
from pyrosetta.rosetta.protocols.scoring import Interface

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import *

import re
import os
import sys
import optparse    
from decimal import Decimal

init(extra_options = ["-restore_pre_talaris_2013_behavior", "-mute all", "-constant_seed", "-ignore_zero_occupancy false"])

'''
Parameters:

    pose: (pyrosetta.rosetta.core.pose) protein structure

    partners: (str) partner chain, separated by '_'

Returns:

    pose: (pyrosetta.rosetta.core.pose) protein structure

Does:
    
    Deletes all chains other than the ones stated in partners.

'''

def get_interface_residue(pose, num_res, partners, interface_cutoff, neighbor_cutoff=8.0):
    pose.dump_pdb("test.pdb")
    exlucde_position = [] #the residues to exclude from the mutation list
    atom_position_arr = [] #the array to store xyz position of all atoms
    CB_pos = np.array([]) #the array for xyz position of CB's
    interface_mask = np.zeros(num_res, dtype=bool) #mask for interface residues
    for residue_i in range(1, num_res + 1): #note that PyRosetta is 1 indexing
        if pose.residue(residue_i).name() in ["GLY"]: #PRO, GLY, CYS should be excluded
            exlucde_position.append(residue_i)
        num_atoms_i = pose.residue(residue_i).natoms() #get number of atoms on res i
        residue_position_arr = np.array([-1, -1, -1]) #the array to store xyz pos of all atoms on res i
        if "GLY" in pose.residue(residue_i).name(): #write a dummy position for GLY "CB"
            CB_pos = np.append(CB_pos, [0, 0, 0])
        for atom_i in range(1, num_atoms_i + 1):
            atom_name = re.sub(r'\s+\d+', '', pose.residue(residue_i).atom_name(atom_i)).lstrip() #clean the atom name
            if atom_name.startswith("H"): #ignore H's when looking for close by atoms
                continue
            pos_atom_i = np.array(list(pose.xyz(pyrosetta.rosetta.core.id.AtomID(atom_i, residue_i)))) #get xyz position
            if "CB" in atom_name:
                CB_pos = np.append(CB_pos, pos_atom_i) #get xyz position of CB's
            residue_position_arr = np.vstack((residue_position_arr, pos_atom_i)) #update residue position array
        atom_position_arr.append(residue_position_arr[1:]) #update atom position array

    CB_pos = CB_pos.reshape(-1, 3)

    #find the shortest distance between two residues. 
    for i in range(len(atom_position_arr) - 1): 
        for j in range(i + 1, len(atom_position_arr)):
            if pose.pdb_info().chain(i + 1) == pose.pdb_info().chain(j + 1): #skip if on the same chain 
                continue
            min_dist = np.amin(cdist(atom_position_arr[i], atom_position_arr[j])) #find minimum distance between res i and res j 
            if min_dist < interface_cutoff:
                interface_mask[i] = True
                interface_mask[j] = True
                
    for i in range(num_res - 1): 
        for j in range(i + 1, num_res):
            if "GLY" in [pose.residue(i + 1).name(), pose.residue(j + 1).name()]:
                continue
            CB_dist = euclidean(CB_pos[i], CB_pos[j])
            if CB_dist < neighbor_cutoff:
                if pose.pdb_info().chain(i + 1) == pose.pdb_info().chain(j + 1): #skip if on the same chain
                    continue
                else:
                    interface_mask[i] = True
                    interface_mask[j] = True

    #remove excluded res
    for i in exlucde_position:
        interface_mask[i - 1] = False

    return interface_mask

'''
Parameters:

    pdb_filename: (str) the name of the pdb file to load pose from

    partners: (str) partner chain, separated by '_'

    mutant_aa: (str) one letter code for the mutant residue

    neighbor_cutoff: (float) the radius from CB on residue A within which 
                     to look for CB's on another residue

    interface_cutoff: (float) the radius from any heavy atom on residue A within
                      which to look for any other heavy atoms on another residue

    repack_cutoff: (float) the radius from a mutant CB within which to repack other
                   residues (CB-CB distance)

    trials: (int) number of trials to perform
    
    trial_output: (str) prefix of the energy output file
    
    output: (bool) True to output every mutant


Returns:

    None

Does:
    
    Performs alanine scanning by mutating all interface residues and calculate 
    binding energy.

'''
def scanning(pdb_filename, partners, mutant_aa = 'A',
        neighbor_cutoff = 8.0, interface_cutoff = 4.0, repack_cutoff = 6.5,
        trials = 1, trial_output = '', output = False):

    if pdb_filename == "NO_INPUT":
        pdb_filename = input("No PDB file name was read, please enter the PDB file name below:\n")

    if partners == "NO_INPUT":
        partners = input("No partner chains were read, please enter the partner chains below separated by \'_\'\n")

    pdb_code = pdb_filename[3:].split(".")[0]
    global file_directory

    file_directory = 'PyRosettaResults/' + pdb_code + '_' + partners + '/'
    if not os.path.exists(file_directory):
        os.makedirs(file_directory)

    pose = pose_from_file(pdb_filename) #load pose from file


    chains_all = [i[1] for i in list(pyrosetta.rosetta.core.pose.conf2pdb_chain(pose).items())]
    if len(chains_all) != 2:
        sys.exit("ERROR: More than two chains are found in the PDB file\n")
        
    ddG_scorefxn = create_score_function('pre_talaris_2013_standard','score12') #sf setup
    ddG_scorefxn.set_weight(core.scoring.fa_atr, 0.44)
    ddG_scorefxn.set_weight(core.scoring.fa_rep, 0.07)
    ddG_scorefxn.set_weight(core.scoring.fa_sol, 0.32)
    ddG_scorefxn.set_weight(core.scoring.hbond_bb_sc, 0.49) 

    num_res = pose.size()
    interface_mask = get_interface_residue(pose, num_res, partners, interface_cutoff, neighbor_cutoff) #get interface res


    for trial in range( trials ): #perform ala scanning
        ddG_mutants = {}
        for i in range(1, pose.total_residue() + 1):
            if interface_mask[i - 1]:
                filename = ''
                if output:
                    filename = pose.pdb_info().name()[:-4] + '_' +\
                        pose.sequence()[i-1] +\
                        str(pose.pdb_info().number(i)) + '->' + mutant_aa #get the file name for mutants
                ddG_mutants[i] = interface_ddG(pose, i, mutant_aa, 
                    ddG_scorefxn, repack_cutoff, filename ) #get binding energy

        #output results
        residues = list( ddG_mutants.keys() )  
        display = [pose.pdb_info().chain(i) + " " + str(pose.pdb_info().number(i)) 
        + " " + pose.sequence()[i - 1] + '\t' + str(ddG_mutants[i]) + '\n' for i in residues]

        f = open(file_directory + trial_output + '_' + str(trial + 1) + '.txt' , 'w' )
        f.writelines(display)
        f.close()

    scanning_analysis(trial_output)

'''
Parameters:
    
    pose: (pyrosetta.rosetta.core.pose) protein structure

    mutant_position: (int) position of mutation (1 indexing)

    mutant_aa: (str) one letter code of the target mutation amino acid

    scorefxn: (pyrosetta.rosetta.core.scoring.ScoreFunction) score function

    cutoff: (float) repack cutoff

    out_filename: (str) name of the result output file

Return: 
    
    ddg: (float) binding energy

Does:

    Mutate the residue at mutant_position to mutant_aa and calculate the 
    binding energy difference

'''
def interface_ddG( pose, mutant_position, mutant_aa, scorefxn,
        cutoff, out_filename = ''):
    # 1. create a reference copy of the pose
    wt = Pose()    # the "wild-type"
    wt.assign(pose)

    # 2. create a copy of the pose for mutation
    mutant = Pose()
    mutant.assign(pose)

    # 3. mutate the desired residue
    # the pack_radius argument of mutate_residue (see below) is redundant
    #    for this application since the area around the mutation is already
    #    repacked

    pyrosetta.toolbox.mutants.mutate_residue(mutant, mutant_position, 
                                                 mutant_aa, 0.0, scorefxn)
    # 4. calculate the "interaction energy"
    # the method calc_interaction_energy is exposed in PyRosetta however it
    #    does not alter the protein conformation after translation and may miss
    #    significant interactions
    # an alternate method for manually separating and scoring is provided called
    #    calc_binding_energy (see Interaction Energy vs. Binding Energy below)

    wt_score = calc_binding_energy(wt, scorefxn,
        mutant_position, cutoff)
    mut_score = calc_binding_energy(mutant, scorefxn,
        mutant_position, cutoff)

    ddg = mut_score - wt_score

    # -write the mutant structure to a PDB file
    if out_filename:
        mutant.dump_pdb(out_filename)

    return ddg

'''
Parameters:
    
    pose: (pyrosetta.rosetta.core.pose) protein structure

    scorefxn: (pyrosetta.rosetta.core.scoring.ScoreFunction) score function

    center: (int) position of mutation (1 indexing)

    cutoff: (float) repack cutoff

Return: 
    
    bindin energy : (float) binding energy

Does:

    Calculate binding energy

'''

def calc_binding_energy(pose, scorefxn, center, cutoff):
    # create a copy of the pose for manipulation
    test_pose = Pose()
    test_pose.assign(pose)

    # setup packer options
    # the sidechain conformations of residues "near the interface", defined as
    #    within  <cutoff>  Angstroms of an interface residue, may change and
    #    must be repacked, if all residues are repacked, aberrant sidechain
    #    conformations near the interface, but independent of complex
    #    interactions, will be repacked for the mutant and wild-type structures
    #    preventing them from adding noise to the score difference
    # this method of setting up a PackerTask is different from packer_task.py
    tf = standard_task_factory()    # create a TaskFactory
    tf.push_back(core.pack.task.operation.RestrictToRepacking())    # restrict it to repacking

    # this object contains repacking options, instead of turning the residues
    #    "On" or "Off" directly, this will create an object for these options
    #    and assign it to the TaskFactory
    prevent_repacking = core.pack.task.operation.PreventRepacking()

    # the "center" (nbr_atom) of the mutant residue, for distance calculation
    center = test_pose.residue(center).nbr_atom_xyz()
    for i in range(1, test_pose.total_residue() + 1):
        # the .distance_squared method is (a little) lighter than .norm
        # if the residue is further than <cutoff> Angstroms away, do not repack
        if center.distance_squared(
                test_pose.residue(i).nbr_atom_xyz()) > cutoff**2:
            prevent_repacking.include_residue(i)

    # apply these settings to the TaskFactory
    tf.push_back(prevent_repacking)

    # setup a PackRotamersMover
    packer = protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)

    # repack the test_pose
    packer.apply(test_pose)


    # score this structure
    before = scorefxn(test_pose)

    #
    xyz = rosetta.numeric.xyzVector_double_t()    # a Vector for coordinates
    xyz.x = 500.0    # arbitrary separation magnitude, in the x direction
    xyz.y = 0.0    #...I didn't have this and it defaulted to 1e251...?
    xyz.z = 0.0    #...btw thats like 1e225 light years,
                   #    over 5e245 yrs at Warp Factor 9.999 (thanks M. Pacella)


    chain2starts = len(pose.chain_sequence(1)) + 1
    for r in range(chain2starts, test_pose.total_residue() + 1):
        for a in range(1, test_pose.residue(r).natoms() + 1):
            test_pose.residue(r).set_xyz(a,
                test_pose.residue(r).xyz(a) + xyz)

    # repack the test_pose after separation
    packer.apply(test_pose)

    # return the change in score
    return before - scorefxn(test_pose)


'''
Parameters:
    
    trial_output: (str) output file name

    hot_cutoff: (float) threshold for a res to be a hot spot

Returns:

    None

Does:
    
    Performs results analysis. Correct numbers and output final results
'''

def scanning_analysis(trial_output, hot_cutoff=0.6):
    global file_directory
    filenames = os.listdir(file_directory)
    filenames =[(file_directory + i) for i in filenames if trial_output in i]

    # perform an initial reading, to setup lists
    filename = filenames[0]
    f = open(filename , 'r')
    data = f.readlines()
    data = [i.strip() for i in data]    # remove "\n"
    f.close()

    # list of mutations, should be identical for all output files
    mutants = [i.split('\t')[0] for i in data]
    ddg = [float(i.split('\t')[1]) for i in data]

    # for all files beyond the first, add the "ddG" values
    for filename in filenames[1:]:
        f = open( filename , 'r' )
        data = f.readlines()
        data = [i.strip() for i in data]
        f.close()

        ddg = [float(data[i].split('\t')[1])+ddg[i] for i in range(len(data))]

    # average by dividing by the number of files
    ddg = [i/len(filenames) for i in ddg]

    for index, i in enumerate(ddg):
        if i >= 4.5: #cap at 4.5
            ddg[index] = 4.5
        elif i <= -2: #large negative, change to positive
            ddg[index] = 1
        elif i <= 0: #small negative, change to zero
            ddg[index] = 0

    # extract list elements (for ddg and thus mutants) with a ddg value more
    #    than 1 standard deviation away
    significant = [i for i in range(len(ddg)) if ddg[i] > hot_cutoff]
    
    # these are the hotspots
    hotspots = [mutants[i] for i in significant]

    with open(file_directory + trial_output + "_total.txt", "w+") as fo:
        fo.write("="* 80)
        fo.write("\n")
        fo.write( 'All Interface Residues\n' )
        for index, mutant in enumerate(mutants):
            fo.write("%5s%10.2f\n" % (mutant, Decimal(ddg[index])))
        fo.write("="* 80)
        fo.write("\n")
        fo.write( 'Likely Hotspot Residues\n' )
        for hotspot in hotspots:
            fo.write( hotspot + "\n" )

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('--pdb_filename', dest = 'pdb_filename',
        default = 'NO_INPUT',   
        help = 'the PDB file containing the protein')
    parser.add_option('--partners', dest = 'partners',
        default = 'NO_INPUT',    
        help = 'the relative chain partners for docking')
    parser.add_option('--mutant_aa', dest = 'mutant_aa',
        default = 'A',   
        help = 'the amino acid to mutate all residues to')
    parser.add_option('--neighbor_cutoff', dest = 'neighbor_cutoff',
        default = '8.0',   
        help = 'the distance (in Angstroms) to detect neighbors')
    parser.add_option('--trials', dest='trials',
        default = '20',    
        help = 'the number of trials to perform')
    parser.add_option('--trial_output', dest = 'trial_output',
        default = 'ddG_out',    
        help = 'the name preceding all output files')
    parser.add_option('--interface_cutoff', dest = 'interface_cutoff',
        default = '4.0',    
        help = 'the distance (in Angstroms) to detect interface residues')
    parser.add_option('--repack_cutoff', dest = 'repack_cutoff',
        default = '6.5',    
        help = 'the distance (in Angstroms) to detect residues for repacking\
            near the interface')
    parser.add_option('--output', dest = 'output',
        default = '',    
        help = 'if True, mutant structures are written to PDB files')

    (options,args) = parser.parse_args()


    # PDB file option
    pdb_filename = options.pdb_filename
    partners = options.partners
    # scanning options
    mutant_aa = options.mutant_aa
    neighbor_cutoff = float(options.neighbor_cutoff)
    interface_cutoff = float(options.interface_cutoff)
    repack_cutoff = float(options.repack_cutoff)
    # trials options
    trials = int(options.trials)
    trial_output = options.trial_output
    output = bool(options.output)

    scanning(pdb_filename, partners, mutant_aa,
            neighbor_cutoff, interface_cutoff, repack_cutoff,
            trials, trial_output, output)

