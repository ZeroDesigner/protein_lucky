function bestrank(){
	for file in `ls $input`
		do
		if [ -d  $input"/"$file ]; then			
        		if [ ! -f "bestranking.lst" ]; then
           			echo 'find a  bestrank.lst'
	   			pwd
	   			cp $input/$file/bestranking.lst $input/collection-of-bestrank/bestrank-$file.lst
			fi
		fi
	
			done
}
function deallst() {
	for lst in `ls $col `
		do 
		sed -i '1,6d' ./collection-of-bestrank/$lst
	done
	cat  $input/collection-of-bestrank/*.lst > result.lst
 }

mkdir "collection-of-bestrank"
echo "input the rootdir"
read input;
bestrank $input
col=$input"/"collection-of-bestrank
deallst $col
