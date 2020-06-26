# ala_scan
Created by Tim Ling @ YSL lab

## Dependencies

### Python
* numpy
* matplotlib
* scipy

### PyRosetta
* [pyrosetta 4](http://www.pyrosetta.org/dow) installed through conda (pyrosetta-2020.02+release.22ef835b4a2-py37_0)


## Reference
	
	Gavenonis, Jason, et al. “Comprehensive Analysis of Loops at Protein-Protein Interfaces for Macrocycle Design.” 
		Nature Chemical Biology, vol. 10, no. 9, 2014, pp. 716–722., doi:10.1038/nchembio.1580.
	
	Kortemme, T., and D. Baker. “A Simple Physical Model for Binding Energy Hot Spots in Protein-Protein Complexes.” 
		Proceedings of the National Academy of Sciences, vol. 99, no. 22, 2002, pp. 14116–14121., doi:10.1073/pnas.202485799.
	
	Kortemme, T., et al. “Computational Alanine Scanning of Protein-Protein Interfaces.” 
		Science Signaling, vol. 2004, no. 219, 2004, doi:10.1126/stke.2192004pl2.
	
	Siegert, Timothy R., et al. “Identifying Loop-Mediated Protein–Protein Interactions Using LoopFinder.” 
		Methods in Molecular Biology Modeling Peptide-Protein Interactions, 2017, pp. 255–277., doi:10.1007/978-1-4939-6798-8_15.	

## Protocol Description

	1. Identify the interface residues. A resiude is on the interface if:
		a. A residue has a side chain having at least one atom 
		   within a sphere with 4Å radius of an atom belonging 
		   to the other partner in the complex.

		b. A residue that becomes significantly buried upon 
		   complex formation, as measured by an increase in the
		   number of CB atoms within a sphere with a radius of 
		   8Å around the CB atom of the residue of interest.
	
	2. Mutate the interface residue one by one to a target amino acid (alanine by defulat).

	3. Calculate the binding energy of both the wt and mt; the binding energy can be expressed by:
	
		E(A*B*) - E(A + B)

	   where A*B* is the complex and A + B is the refolded monomers.
	   All side chains of residues within 6.5A from the site of mutation is refolded. 
	   Note that both the wt and the mt are refolded.

	4. Calculate ddG which is the binding energy difference between the mt and wt.

## Quickstart
	
	To perform computational Alanine scan, do 

		python ala_scan.py --pdb_filename test.pdb --partners chainA_chainB 
		[--mutant_aa A] [--neighbor_cutoff 8.0] 
		[--trials 20] [--trial_output ddG_out] 
		[--interface_cutoff 4] [--repack_cutoff 6.5] [--output 0]

	Description of the flags:

		--pdb_filename: the name of the PDB file to load pose from. 

		--partners: the interacting chains. (Note: there can only be two chains.)

		--mutant_aa: the target mutation amino acid

		--neighbor_cutoff: the cutoff for interface criteria 1b. 

		--trials: the number of trials to perform.

		--trial_output: the name of the files to write output in.

		--interface_cutoff: the cutoff for interface criteria 1a. 

		--repack_cutoff: the cutoff for residue repacking from the site of mutation.

		--output: if to write each mutant into PDB files.

## Conducting Large-Scale Computational Alanine Scan
	
	Because ala_scan.py only executes alanine scan for one PDB structure,
	repeating the same commands for many structures can be tedious.
	Four helper scripts are provided for convience.
	
	1. parser.py:
	   Simply git clone the entire ala_scan folder to the directory where multiple 
	   PDB files are located. Then go to the ala_scan folder and use
	
		python parser.py 

	  The script will search for all the PDB files in the upper directory and conducts 
	  alanine scan for all the structures.

	2. comparison.py:
	   This script looks at all the results files and generates a spreadsheet of
	   results. -1 indicates that no hot spot residues were found, 
	   0 indicates that there are hot spot residues but no hot loop was found,
	   and 1 indicates that hot loops were found.

	
	3. pickout_hotloop.py and identify_hotloop.py:
	   These two scripts filters the results and writes the important hot loops into a file.
## Acknowledgments
The implementation utilizes code from the following:
* [Evan H. Baugh's D090_Ala_scan.py](https://graylab.jhu.edu/pyrosetta/downloads/scripts/demo/D090_Ala_scan.py)





