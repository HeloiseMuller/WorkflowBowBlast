# WorkflowBowBlast

WorkflowBowBlast runs automaticaly several tools to process raw Illumina reads, such as trimmomatic, bowtie2, bedtools genomecov and blastn.

# 1) Installation

git clone https://github.com/HeloiseMuller/WorkflowBowBlast

# 2) Running the workflow

## Overview

It is advised to create a new directory for each project, WorkflowBowBlast will create several directories for the outputs. 

To use WorkflowBowBlast, one needs to set the configuration file which contains the location of the data, the choice of the tools to run with their location and some other parameters. 

Then run the main script with
`./main.sh config.txt`

## Config file

Here is an example of config.txt
```
#Directory of the project
dir=/home/user/Project

#Name of the run of Illumina
rawIlumina=X22DC1_1.fq.gz

#Name of the sample
sample=data

#Reference genome of the sample
path_species=/home/user/ref/genome1.fasta 
species=genome1 #no dot in species name

#Reference species of a second species
path_species2=/home/user/ref/genome2.fasta  #Comment if does not need it
species2=genome2 #no dot i species name

#Scripts to run
integrity=TRUE
trimmomatic=TRUE
bowtie2=TRUE
coverage=TRUE
blastn=TRUE
parallele=TRUE

#Number of threads ​
60
the name of the Illumina run. It needs to be in `$dir`. The raw Illumina reads needs to be organized as follow: `$dir/$runIlumina/raw_data/`to use (even number if parallel=TRUE)
threads=20


#Path of softwares to run


```

`dir=/home/user/Project` is the directory of your project. All the outputs will be created in this directory. 

`rawIlumina=X22DC1_1.fq.gz` is the path of the raw Illumina forward reads. No need to specify the reverse reads, they need to be at the same location.

`sample=data` is the name of the sample on which one wants to run the workflow. The fq.gz files need to begin with `sample` and it cannot contain any dots.

`path_species=/home/user/ref/genome1.fasta` is the path of the reference genome of the Illumina reads. The name of the species used in the outputs is `$species`. `$species` cannot contains any dots

`path_species2=/home/user/ref/genome2.fasta` is the path of a reference genome of another species of interest. The name of this other species used in the outputs is `$species2`. `$species2` cannot contains any dots. This second species should be used if one wants to look for chimeric reads between `$species` and `$species2`. If not, comment these options.

`integrity`, `trimmomatic`, `bowtie2`, `coverage` and `blastn` are the tools WorkflowBowBlast can run. They are detailed into more details thereafter. Set TRUE to run, FALSE to ignore them. `parallele=TRUE` allows WorkflowBowBlast run blastn without waiting the end of bowtie2 and coverage.

`threads=20` set the total number of threads to use. If `bowtie2=TRUE` and `blastn=TRUE`, `threads` will be split between them if they run at the same time.

## Steps of the workflow

### Integrity

This option calculate the md5 of the raw sequencing reads and run fastqc

### Trimmomatic

This option runs trimmomatic on the raw sequencing reads. For the following, one needs trimmed reads, so set `trimmomatic=TRUE`. If trimmomatic has already been ran, one can set `trimmomatic=FALSE` but check that the trimmed reads are located in $dir/trimmed_data/ and that they are named as follow: `$sample_trimmed.1P.fq.gz`, `$sample_trimmed.2P.fq.gz`, `$sample_trimmed.1U.fq.gz` and `$sample_trimmed.2U.fq.gz.

### Bowtie2

This option runs bowtie2 on the trimmed reads against `genome`. If `genome2` is specified, it also run bowtie2 on this other genome independently. After bowtie2, the resulting sam files are turn into sorted bam files and merged together.

### Coverage

This option runs Bedtools genomecov on the output(s) of bowtie2. The output is a three colomns file, contening the sequencing depth at each position of the genome. Coverage also automatically calcluates the average sequencing depth over the genome and the percentage of genome covered at least once by the Illumina reads.

### Blastn

This option can be ran at the same time as bowtie2 if `parallele=TRUE`. In this case, `threads` will be equally shared between blastn and bowtie2.
This option firslty processes the trimmed fq.gz files to obtain one fasta file. If `genome2` is specified, blastn is ran a first time on the fasta file against this `genome2`. Then, only reads that mapped on `genome2` are mapped on `genome1`.If `genome2` is commented, all reads of the fasta file are mapped on `genome1`.

