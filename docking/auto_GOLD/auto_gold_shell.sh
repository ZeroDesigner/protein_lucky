#!/bin/sh
conffile=gold.conf
dir=`pwd`
# need to get this here to work out how many ligands it's got
file=`awk '/ligand_data_file/ { print $2 }' $conffile`
# Calculate number of molecules to be docked
n_mols=`grep '$$$$' $file | wc -l`
echo " $n_mols  $file "
# process you want
n_proc=25
echo " $n_proc "
#计算每个conf文件的ligands数目如果有余数
split=` expr $n_mols / $n_proc `
echo " $split "
#进入循环
proc=1 # docking counter
lig1=1 # start ligand
lig2=$split # end ligand
echo " $proc $lig1 $lig2 "
for((proc=1 ; $proc <= $n_proc ; )) ; do
        if [ $proc -eq $n_procs ]; then
        lig2=$n_mols # last process gets whatever's left
        fi
        echo " $lig1 $lig2 $proc "
        mkdir ligand-$lig1-$lig2
        cp gold.conf ligand-$lig1-$lig2/
        cd ligand-$lig1-$lig2/
        sed -i "s/start_at_ligand 1/start_at_ligand $lig1 finish_at_ligand $lig2/g" gold.conf
        echo " $lig1 $lig2  "
        setsid gold_auto  &
        cd ..
        proc=`expr $proc + 1`
        lig1=`expr $lig2 + 1`
        lig2=`expr $split \* $proc`
done
echo 'done'
