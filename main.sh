#Get varaibles from config file
. $1

#Get functions
. functions.sh

#Check integrity if needed
if [ $integrity = TRUE ];
then
	echo -e "\nCalculate the md5"
	for i in $dir/$runIllumina/raw_data/$sample/*.fq.gz
	do
		md5sum $i >> $dir/md5.txt
	done
	
	echo -e "\nFasqc on raw data"
	mkdir -p $dir/fastqc
	for i in $dir/$runIllumina/raw_data/$sample/*.fq.gz
	do
		/opt/FastQC-0.11.8/fastqc $i -o $dir/fastqc
	done
fi

#Run trimmomatic if needed
if [ $trimmomatic = TRUE ];
then
	echo -e "Running trimmomatic..."
	mkdir -p $dir/trimmed_data
	java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
			-threads $threads -trimlog $dir/trimmed_data/${sample}_trimmings.log \
			-basein $dir/$runIllumina/raw_data/${sample}*1.fq.gz \
			-baseout $dir/trimmed_data/${sample}_trimmed.fq.gz \
			ILLUMINACLIP:$dir/trimmed_data/adapters.fa:2:30:10 \
			LEADING:20 TRAILING:20 \
			SLIDINGWINDOW:4:15 \
			MINLEN:36
fi

#Split threads between bowtie2 and blastn if use both
if [ $bowtie2 = TRUE ] &&  [ $blastn = TRUE ];
then    
        threads_bowtie2=`expr $threads / 2`
        threads_blastn=`expr $threads / 2`
elif [ $bowtie2 = TRUE ] && [ $blastn = FALSE ];
then
        threads_bowtie2=$threads
else
        threads_blastn=$threads
fi 


#Run bowtie2 if needed
if [ $bowtie2 = TRUE ];
then
	mkdir -p $dir/Coverage
	func_bowtie2 & PIDbowtie2=$!
fi

#Run blastn if needed
if [ $blastn = TRUE ];
then
	if [ ! -f $dir/trimmed_data/${sample}_cat.fasta ];
	then
		#Get fasta files
		for i in $dir/trimmed_data/${sample}_*.fq.gz
		do
			func_process_fasta $i & PIDfasta=$!
	
		done

		wait $PIDfasta

		#Add "UNPAIRED" at the end of unpaired read names to keep info during blastn
		for i in $dir/trimmed_data/${sample}_*U.fasta
        	do
                	base=`basename $i | cut -d "." -f1`
                	sed '/^>/ s/$/_UNPAIRED/g' $i > $dir/trimmed_data/${base}.fasta3
                	mv $dir/trimmed_data/${base}.fasta3 $dir/trimmed_data/${base}.fasta
        	done


		#Concatenate fasta files
		cat $dir/trimmed_data/${sample}_*.fasta > $dir/trimmed_data/${sample}_trimmed_cat.fasta
		rm $dir/trimmed_data/${sample}_*.fasta
	fi

	#Run blasn (on both species if there is a second one, keeping only reads in commom)
	func_blastn $dir/trimmed_data/${sample}_trimmed_cat.fasta & PIDblastn=$!

fi

echo -e "\nEnd of the workflow"
