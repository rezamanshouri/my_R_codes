while read line ; do
	
	gene=$(echo $line | cut -d ' ' -f 2)	
	echo -n $gene": "; grep $gene 3_removed_SR | cut -d "," -f 2 | sort | xargs | tr ' ' ,

done <5_non_uniq_frequencies




