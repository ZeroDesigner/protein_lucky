simple way to use  rosetta_loop_model
1:creat the model_1_miss.pdb,
model_1_miss.pdb is the loop_miss protein which you need to model
then run command ' score_jd2.static.linuxgccrelease -in:file:s model_1_miss.pdb -renumber_pdb true -out:pdb -overwrite
'
OK,now you get the pdb_file --'model_1_miss_0001.pdb '
then make a file named miss_loop.file
2.  Create the loop file.

    Create a new file, called miss_loop.file, and open it in your text editor.
    Insert the following line:

        LOOP 5 8 0 5 1

    This excerpt from the loopmodel documentation describes the meaning of the
    6 columns in that line:

        column1  "LOOP":     Literally the string LOOP, identifying this line as a loop
                             In the future loop specification files may take other data.
        column2  "integer":  Loop start residue number
        column3  "integer":  Loop end residue number
        column4  "integer":  Cut point residue number, >=startRes, <=endRes.
        column5  "float":    Skip rate. default - never skip (0)
        column6  "boolean":  Extend loop. Set to 1 
3:loopmodel
run the command 'loopmodel.static.linuxgccrelease @flags'
you will get model_1_miss_0001_0001.pdb