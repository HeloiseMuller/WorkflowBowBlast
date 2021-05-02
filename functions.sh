function func_check {
	set -e
	if [ $integrity = TRUE ];
	then
		$path_fastqc -version
	fi
	if [ $trimmomatic = TRUE ];
	then
		echo -e "\nTrimmomatic v`java -jar $path_trimmomatic PE -version`"
	fi
	if [ $bowtie2 = TRUE ];
	then
		$path_bowtie2/bowtie2 --version | grep "bowtie2"
		$path_samtools --version | grep "samtools"
	fi
	if [ $coverage = TRUE ];
	then
		$path_bedtools --version
	fi
	if [ $blastn = TRUE ];
	then
		$path_blastn -version
	fi
	set +e
}


function func_process_fasta {
	echo -e "\ngunzip on $i"
        gunzip -k $i
        echo -e "\nmake fasta from $i"
	base=`basename $i | cut -d "." -f1`
        seqtk seq -a $dir/trimmed_data/${base}.fq > $dir/trimmed_data/${base}.fasta
	rm $dir/trimmed_data/${base}.fq
	#Add _ in names fasta to keep info about F/R during blast:
	sed -e "s/ /_/g" $dir/trimmed_data/${base}.fasta > $dir/trimmed_data/${base}.fasta2 
	mv $dir/trimmed_data/${base}.fasta2 $dir/trimmed_data/${base}.fasta
}

function func_blastn {

	if [ ! -f ${path_species}.nhr ];
	then
	echo -e "\nRunning blast indexing on $species ..."
	makeblastdb -in $path_species \
		-dbtype nucl \
		-parse_seqids \
		-out $path_species \
		-logfile $dir/IndexingBlast_$species.log
	fi

	#If species2 is not specified, map all the reads on species
	if [  -z $species2 ];
	then
		echo -e "\nblastn of all the reads of ${sample} on ${species}"
        	$path_blastn -query $dir/trimmed_data/${sample}_trimmed_cat.fasta \
                	-db $path_species \
                	-outfmt 6 \
                	-max_target_seqs 2 \
                	-out $dir/blastn/${sample}_trimmed_vs_${species}.txt \
                	-num_threads ${threads_blastn}
	
	#If species2 is specified, map on species2 first, then take only the reads that mapped to map on species
	else	
	        if [ ! -f ${path_species2}.nhr ];
	        then
        	echo -e "\nRunning blast indexing on $species2 ..."
        	makeblastdb -in $path_species2 \
                	-dbtype nucl \
                	-parse_seqids \
                	-out $path_species \
                	-logfile $dir/IndexingBlast_$species2.log
        	fi

		
		echo -e "\nblastn of all the reads of ${sample} on ${species2}"
		$path_blastn -query $dir/trimmed_data/${sample}_trimmed_cat.fasta \
			-db $path_species2 \
                	-outfmt 6 \
                	-max_target_seqs 2 \
                	-out $dir/blastn/${sample}_trimmed_vs_${species2}.txt \
                	-num_threads ${threads_blastn}

		#extract 1st colomn where name sequences to keep | keep each read once
		cut -f1 $dir/blastn/${sample}_trimmed_vs_${species2}.txt | sort | uniq > $dir/blastn/${sample}_trimmed_vs_${species2}.lst
		#make fasta in which keep only reads that blasted against species2
		seqtk subseq $dir/trimmed_data/${sample}_trimmed_cat.fasta $dir/blastn/${sample}_trimmed_vs_${species2}.lst > $dir/blastn/${sample}_trimmed_vs_${species2}.fasta
		echo -e "\nblastn of the reads of ${sample} that mapped on ${species} on ${species}"
		$path_blastn -query $dir/blastn/${sample}_trimmed_vs_${species2}.fasta \
			-db $path_species \
			-outfmt 6 \
			-max_target_seqs 2 \
			-out $dir/blastn/${sample}_trimmed_vs_${species2}_vs_${species}.txt \
			-num_threads ${threads_blastn}
	fi

}

function func_bowtie2 {
	if [ ! -f ${path_species}.1.bt2* ];
        then
                echo -e "\nBowtie2 indexing on $species ..."
                $path_bowtie2/bowtie2-build --threads ${threads_bowtie2} $path_species $path_species >> $dir/bowtie2Indexing.out
        fi

        echo -e "\nRunning bowtie2 on $species ..."
        $path_bowtie2/bowtie2 -x $path_species \
                -1 $dir/trimmed_data/${sample}_trimmed_1P.fq.gz \
                -2 $dir/trimmed_data/${sample}_trimmed_2P.fq.gz \
                --threads ${threads_bowtie2} | \
		$path_samtools view -bS -@ ${threads_bowtie2} > $dir/Bowtie2/${sample}_trimmed_vs_${species}_paired.bam

        $path_bowtie2/bowtie2 -x $path_species \
                        -U $dir/trimmed_data/${sample}_trimmed_1U.fq.gz,$dir/trimmed_data/${sample}_trimmed_2U.fq.gz \
                        --threads ${threads_bowtie2} | \
                	$path_samtools view -bS -@ ${threads_bowtie2} > $dir/Bowtie2/${sample}_trimmed_vs_${species}_unpaired.bam


        if [ ! -z $species2 ];
        then
                if [ ! -f ${path_species2}.1.bt2* ];
                then
                        echo -e "\nBowtie2 indexing on $species2 ..."
                        $path_bowtie2/bowtie2-build --threads $threads $path_species2 $path_species2 >> $dir/bowtie2Indexing.out
                fi

                echo -e "\nRunning bowtie2 on $species2 ..."
                $path_bowtie2/bowtie2 -x $path_species2 \
                -1 $dir/trimmed_data/${sample}_trimmed_1P.fq.gz \
                -2 $dir/trimmed_data/${sample}_trimmed_2P.fq.gz \
                --threads ${threads_bowtie2}  | \
                $path_samtools view -bS -@ ${threads_bowtie2} > $dir/Bowtie2/${sample}_trimmed_vs_${species2}_paired.bam

                $path_bowtie2/bowtie2 -x $path_species2 \
                        -U $dir/trimmed_data/${sample}_trimmed_1U.fq.gz,$dir/trimmed_data/${sample}_trimmed_2U.fq.gz \
                        --threads ${threads_bowtie2} | \
                $path_samtools view -bS -@ ${threads_bowtie2} > $dir/Bowtie2/${sample}_trimmed_vs_${species2}_unpaired.bam
        fi 

	echo -e "\nSorting bam files ..."
	unsortedBams=`ls $dir/Bowtie2/${sample}_trimmed_vs_*.bam | grep -v '_sorted.bam'`
	for i in $unsortedBams
	do
		base=`basename $i | cut -d "." -f1`
		$path_samtools sort -o $dir/Bowtie2/${base}_sorted.bam -@ ${threads_bowtie2} $i
		rm  $i
	done

	echo -e "\nMerge paired and unpaired"
	$path_samtools merge $dir/Bowtie2/${sample}_trimmed_vs_${species}_merged_sorted.bam $dir/Bowtie2/${sample}_trimmed_vs_${species}_*_sorted.bam
	if [ ! -z $species2 ];
        then
	        $path_samtools merge $dir/Bowtie2/${sample}_trimmed_vs_${species2}_merged_sorted.bam $dir/Bowtie2/${sample}_trimmed_vs_${species2}_*_sorted.bam
	fi
	rm $dir/Bowtie2/${sample}_trimmed_vs_*paired_sorted.bam
	
	if [ $coverage = TRUE ];
	then
		$path_bedtools --version
		func_cov
	fi

}

function func_cov {

mkdir -p $dir/Coverage

echo -e "\nCalculing coverage with bedtools..."
$path_bedtools genomecov -ibam $dir/Bowtie2/${sample}_trimmed_vs_${species}_merged_sorted.bam \
                       -g $path_species -d > $dir/Coverage/${sample}_trimmed_vs_${species}_coverage_positions
length_species=`awk '!/^>/{l+=length($0)}END{print l}' $path_species`
cat $dir/Coverage/${sample}_trimmed_vs_${species}_coverage_positions | awk '{sum+=$3} END {print "Average coverage of ${samples} on ${species} = ",sum/NR}' 
cat $dir/Coverage/${sample}_trimmed_vs_${species}_coverage_positions | awk '$3!=0' | awk 'END {print "Proportion of ${species} genome covered by ${sample} =",NR/$length_species}' 

if [ ! -z $species2 ];
then
	$path_bedtools genomecov -ibam $dir/Bowtie2/${sample}_trimmed_vs_${species2}_merged_sorted.bam \
                       -g $path_species2 -d > $dir/Coverage/${sample}_trimmed_vs_${species2}_coverage_positions
	length_species2=`awk '!/^>/{l+=length($0)}END{print l}' $path_species2`
	cat $dir/Coverage/${sample}_trimmed_vs_${species2}_coverage_positions | awk '{sum+=$3} END {print "Average coverage of ${samples} on ${species2} = ",sum/NR}' 
	cat $dir/Coverage/${sample}_trimmed_vs_${species2}_coverage_positions | awk '$3!=0' | awk 'END {print "Proportion of ${species2} s genome covered by ${sample} =",NR/$length_species2}' 
fi

}
