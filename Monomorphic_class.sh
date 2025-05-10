invs="$(cat invappend2_v2.in)" #store in an object the inversions from the file invappend2_v2.in

for inv in $invs
do
    inv_col="$(head -1 data/Porubsky_invs_callset_GT_no_lowconf_v4.csv | tr "\t" "\n" | cat -n | grep "${inv}" | cut -f1 | tr -d " ")" #find the column number of the inversion in the genotypes file. 
    chr="$(cat data/Inversions_coordinates_porubsky_v3.csv | grep "$inv[[:space:]]" | awk '{print $2}' | sed -e 's/chr//g' | head -n 1)" #find the chromosome of the inversion in the genome coordinates file. 
    type="$(cat data/Inversions_coordinates_porubsky_v3.csv | grep "$inv[[:space:]]" | awk '{print $7}' | sed -e 's/chr//g' | head -n 1)" #find the type of the inversion in the genome coordinates file.
    if [ $chr == "Y" ]
    then
    
        echo "$inv,Chr Y inversion,$type" >> Genotyping_info.txt #if the inversion is on the Y chromosome, write it to the file Genotyping_info.txt
    else
        only_genotyped_rows="$(cut -f$inv_col data/Porubsky_invs_callset_GT_no_lowconf_v4.csv | head -n -1 | tail -n+2 | awk '$1!="NA"' | awk '$1!="ND"' | awk '$1!="NER"' | awk '$1!="./."')" #for each inversion extract only the genotyped individuals. 
        count_only_genotyped_rows="$(echo "$only_genotyped_rows" | wc -l)" #count the number of genotyped individuals.
        if [ $count_only_genotyped_rows -eq 1 ] 
        then
            echo "$inv,only 1 sample genotyped,$type" >> Genotyping_info.txt # if only one individual is genotyped, write it to the file Genotyping_info.txt
        else
            genotypes="$(echo "$only_genotyped_rows"| sed 's/\//\n/g' | sort | uniq)" #if more then one individual is genotyped, extract all alleles for that inversion. 
            count_genotypes="$(echo "$genotypes" | wc -l)" #count the number of alleles (Std or Inv only possible). 
            if [ $count_genotypes -eq 1 ] 
            then
            echo "$inv,monomorphic,$type" >> Genotyping_info.txt # if only one allele is present, write it to the file Genotyping_info.txt as a monomoprhic inversion.
            else
                echo "$inv,polymorphic,$type" >> Genotyping_info.txt #if none of the above conditions is true, the inversion is polymorphic. 
            fi
        fi
    fi
done 


