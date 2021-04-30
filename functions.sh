function func_process_fasta {
	echo -e "\ngunzip on $i"
        gunzip -k $i
        echo -e "\nmake fasta "
	base=`basename $i | cut -d "." -f1`
        seqtk seq -a $dir/trimmed_data/${base}.fq > $dir/trimmed_data/${base}.fasta
	rm $dir/trimmed_data/${base}.fq
	#Add _ in names fasta to keep info about F/R during blast:
	sed -e "s/ /_/g" $dir/trimmed_data/${base}.fasta > $dir/trimmed_data/${base}.fasta2 
	mv $dir/trimmed_data/${base}.fasta2 $dir/trimmed_data/${base}.fasta
}

function func_blastn {

	#If species2 is not specified, map all the reads on species
	if [  -z $species2 ];
	then
		echo -e "\nblastn of all the reads of ${sample} on ${species}"
        	blastn -query $dir/trimmed_data/${sample}_trimmed_cat.fasta \
                	-db $path_species \
                	-outfmt 6 \
                	-max_target_seqs 2 \
                	-out $dir/blastn/${sample}_vs_${species}.txt \
                	-num_threads ${threads_blastn}
	else
		#If species2 isspecified, map on species2 first, then take only the reads that mapped to map on species
		echo -e "\nblastn of all the reads of ${sample} on ${species2}"
		blastn -query $dir/trimmed_data/${sample}_trimmed_cat.fasta \
			-db $path_species2 \
                	-outfmt 6 \
                	-max_target_seqs 2 \
                	-out $dir/blastn/${sample}_vs_${species2}.txt \
                	-num_threads ${threads_blastn}

		#extract 1st colomn where name sequences to keep | keep each read once
		cut -f1 $dir/blastn/${sample}_vs_${species2}.txt | uniq > $dir/blastn/${sample}_vs_${species2}.lst
		#make fasta in which keep only reads that blasted against species2
		seqtk subseq $dir/trimmed_data/${sample}_trimmed_cat.fasta $dir/blastn/${sample}_vs_${species2}.lst > $dir/blastn/${sample}_trimmed_vs_${species2}.fasta
		echo -e "\nblastn of the reads of ${sample} that mapped on ${species} on ${species}"
		blastn -query $dir/blastn/${sample}_trimmed_vs_${species2}.fasta \
			-db $path_species \
			-outfmt 6 \
			-max_target_seqs 2 \
			-out $dir/blastn/${sample}_trimmed_vs_${species2}_vs_${species}.txt \
			-num_threads ${threads_blastn}
	fi

}

function func_bowtie2 {
	if [ ! -f ${path_species}.bt2 ];
        then
                echo -e "\nIndex $species ..."
                /opt/bowtie2-2.3.4.2/bowtie2-build --threads $threads $path_species $path_species
        fi

        echo -e "\nRunning bowtie2 on $species ..."
        /opt/bowtie2-2.3.4.2/bowtie2 -x $path_species \
                -1 $dir/trimmed_data/${sample}_trimmed_1P.fq.gz \
                -2 $dir/trimmed_data/${sample}_trimmed_2P.fq.gz \
                -S $dir/Bowtie2/${sample}_trimmed_vs_${species}_paired.sam \
                --threads ${threads_bowtie2}

        /opt/bowtie2-2.3.4.2/bowtie2 -x $path_species \
                        -U $dir/trimmed_data/${sample}_trimmed_1U.fq.gz,$dir/trimmed_data/${sample}_trimmed_2U.fq.gz \
                        -S $dir/Bowtie2/${sample}_trimmed_vs_${species}_unpaired.sam \
                        --threads ${threads_bowtie2}


        if [ ! -z $species2 ];
        then
                if [ ! -f ${species2}.bt2 ];
                then
                        echo -e "\nIndex $path_species2 ..."
                        /opt/bowtie2-2.3.4.2/bowtie2-build --threads $threads $path_species2 $path_species2
                fi

                echo -e "\nRunning bowtie2 on $species2 ..."
                /opt/bowtie2-2.3.4.2/bowtie2 -x $path_species2 \
                -1 $dir/trimmed_data/${sample}_trimmed_1P.fq.gz \
                -2 $dir/trimmed_data/${sample}_trimmed_2P.fq.gz \
                -S $dir/Bowtie2/${sample}_trimmed_${species2}_paired.sam \
                --threads ${threads_bowtie2}

                /opt/bowtie2-2.3.4.2/bowt${sample}_ie2 -x $path_species2 \
                        -U $dir/trimmed_data/${sample}_trimmed_1U.fq.gz,$dir/trimmed_data/${sample}_trimmed_2U.fq.gz \
                        -S $dir/Bowtie2/${sample}_trimmed_${species2}_unpaired.sam \
                        --threads ${threads_bowtie2}
        fi

	echo -e "\nTurn sam files into sort'ed bam files"
	for i in $dir/Bowtie2/${sample}*.sam
	do
		base=`basename $i | cut -d "." -f1`
		samtools view -bh $i | samtools sort -o $dir/Bowtie2/${base}_sorted.bam
		rm $i
	done

	echo -e "\nMerge paired and unpaired"
	samtools merge $dir/Bowtie2/${samples}_trimmed_vs_${species}_merged_sorted.bam $dir/Bowtie2/${samples}_trimmed_vs_${species}_*_sorted.bam
	if [ ! -z $species2 ];
        then
	        samtools merge $dir/Bowtie2/${samples}_trimmed_vs_${species2}_merged_sorted.bam $dir/Bowtie2/${samples}_trimmed_vs_${species2}_*_sorted.bam
	fi
	rm $dir/Bowtie2/${samples}_trimmed_vs_${species}_*paired_sorted.bam
	
	if [ $coverage = TRUE ]:
	then
		func_cov
	fi

}

function func_cov {

mkdir -p $dir/Coverage

echo -e "\nCalculing coverage with bedtools..."
bedtools genomecov -ibam $dir/Bowtie2/${samples}_trimmed_vs_${species}_merged_sorted.bam \
                       -g $path_species -d > $dir/Coverage/${samples}_trimmed_vs_${species}coverage_positions
length_species=`awk '!/^>/{l+=length($0)}END{print l}' $path_species`
cat $dir/Coverage/${samples}_trimmed_vs_${species}_coverage_positions | awk '{sum+=$3} END {print "Average coverage of ${samples} on ${species} = ",sum/NR}' 
cat $dir/Coverage/${samples}_trimmed_vs_${species}_coverage_positions | awk '$3!=0' | awk 'END {print "Proportion of ${species}'s genome covered by ${sample} =",NR/$length_species}' 

if [ ! -z $species2 ];
then
	bedtools genomecov -ibam $dir/Bowtie2/${samples}_trimmed_vs_${species2}_merged_sorted.bam \
                       -g $path_species2 -d > $dir/Coverage/${samples}_trimmed_vs_${species2}_coverage_positions
	length_species2=`awk '!/^>/{l+=length($0)}END{print l}' $path_species2`
	cat $dir/Coverage/${samples}_trimmed_vs_${species2}_coverage_positions | awk '{sum+=$3} END {print "Average coverage of ${samples} on ${species2} = ",sum/NR}' 
	cat $dir/Coverage/${samples}_trimmed_vs_${species2}_coverage_positions | awk '$3!=0' | awk 'END {print "Proportion of ${species2}'s genome covered by ${sample} =",NR/$length_species2}' 
fi

}
