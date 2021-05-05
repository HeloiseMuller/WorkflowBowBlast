echo -e "\n"

#Get varaibles from config file
. $1

#Get functions
path_WorkflowBowBlast=`dirname $(realpath $0)`
echo $path_WorkflowBowBlast
. $path_WorkflowBowBlast/functions.sh

#Check executables
echo -e "\nVersions of the executables:"
func_check

#Sample
echo -e "\nSample is $sample"

#Check integrity if needed
if [ $integrity = TRUE ];
then
	echo -e "\nCalculate the md5"
	for i in $dir/$runIllumina/raw_data/$sample/*.fq.gz
	do
		md5sum $i >> $dir/md5.txt
	done
	
	echo -e "\nFasqc on raw data"
	$path_fastqc -version
	mkdir -p $dir/fastqc
	for i in $dir/$runIllumina/raw_data/$sample/*.fq.gz
	do
		$path_fastqc $i -o $dir/fastqc
	done
fi

#Run trimmomatic if needed
if [ $trimmomatic = TRUE ];
then
	echo -e "Running trimmomatic..."
	mkdir -p $dir/trimmed_data
	java -jar $path_trimmomatic PE \
			-threads $threads -trimlog $dir/trimmed_data/${sample}_trimmings.log \
			-basein $rawIllumina \
			-baseout $dir/trimmed_data/${sample}_trimmed.fq.gz \
			LEADING:20 TRAILING:20 \
			SLIDINGWINDOW:4:15 \
			MINLEN:36
fi

#Split threads between bowtie2 and blastn if use parallel
if [ $parallel = TRUE ] && [ $bowtie2 = TRUE ] &&  [ $blastn = TRUE ];
then
	threads_bowtie2=`expr $threads / 2`
        threads_blastn=`expr $threads / 2`	
else
	threads_bowtie2=$threads
	threads_blastn=$threads
fi


#Run bowtie2 if needed. IF Coverage also true, will run it automaticly
if [ $bowtie2 = TRUE ];
then
	echo -e "\n"
	mkdir -p $dir/Bowtie2
	func_bowtie2 & PIDbowtie2=$!
	if [ $parallel = FALSE ];
	then    
	        wait $PIDbowtie2
	fi
fi

#If bowtie2 has been run independly, run coverage
if [ $bowtie2 = FALSE ] && [ $coverage = TRUE ];
then
	echo -e "\n"
	func_cov & PIDcov=$! 
fi

#Run blastn if needed
if [ $blastn = TRUE ];
then
	if [ ! -f $dir/trimmed_data/${sample}_trimmed_cat.fasta ];
	then
		echo -e "\n"
		PIDsfasta+= #empty array to put the PID of the loop
		#Get fasta files

		for i in $dir/trimmed_data/${sample}_*.fq.gz
		do
			func_process_fasta $i &
		       	PIDsfasta+=" $!"
	
		done

		wait $PIDsfasta

		#Add "UNPAIRED" at the end of unpaired read names to keep info during blastn
		for i in $dir/trimmed_data/${sample}_*U.fasta
        	do
			base=`basename $i | cut -d "." -f1`
                	sed '/^>/ s/$/_UNPAIRED/g' $i > $dir/trimmed_data/${base}.fasta3
                	mv $dir/trimmed_data/${base}.fasta3 $dir/trimmed_data/${base}.fasta
        	done


		#Concatenate fasta files
		cat $dir/trimmed_data/${sample}_*.fasta > $dir/trimmed_data/${sample}_trimmed_cat.fasta
		rm $dir/trimmed_data/${sample}_*U.fasta $dir/trimmed_data/${sample}_*P.fasta 
	fi

	echo -e "\n"
	#Run blasn (on both species if there is a second one, keeping only reads in commom)
	mkdir -p $dir/blastn
	func_blastn $dir/trimmed_data/${sample}_trimmed_cat.fasta & PIDblastn=$!

fi

wait $PIDbowtie2 $PIDcov $PIDblastn
echo -e "\nEnd of the workflow"
