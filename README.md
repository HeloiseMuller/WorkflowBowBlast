# WorkflowBowBlast

WorkflowBowBlast runs automaticaly several tools to process raw Illumina reads, such as trimmomatic, bowtie2 and blastn.

# 1) Installation

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
runIlumina=X204SC2102112-Z02-G001

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
blastn=TRUE

#Number of threads to use (even number)
threads=20


#Path of softwares to run


```

`dir=/home/user/Project` is the directory of your project. All the outputs will be created in this directory. 

`runIlumina=X204SC2102112-Z02-G001` is the name of the Illumina run. It needs to be in `$dir`. The raw Illumina reads needs to be organized as follow: `$dir/$runIlumina/raw_data/`

`sample=data` is the name of the sample on which one wants to run the workflow. The fq.gz files need to begin with `sample` and it cannot contain any dots.

`path_species=/home/user/ref/genome1.fasta` is the path of the reference genome of the Illumina reads. The name of the species used in the outputs is `$species`. `$species` cannot contains any dots

`path_species2=/home/user/ref/genome2.fasta` is the path of a reference genome of another species of interest. The name of this other species used in the outputs is `$species2`. `$species2` cannot contains any dots. This second species should be used if one wants to look for chimeric reads between `$species` and `$species2`. If not, comment these options.

`integrity`, `trimmomatic`, `bowtie2` and `blastn` are the tools WorkflowBowBlast can run. They are detailed with more details thereafter. Set TRUE to run, FALSE to ignore them. `parallele=TRUE` allows WorkflowBowBlast to run bowtie2 and blastn at the same time.

`threads=20` set the total number of threads to use. If `bowtie2=TRUE` and `blastn=TRUE`, `threads` will be split between them if they run at the same time.

## Steps of the workflow

### Integrity

This option calculate the md5 of the raw sequencing reads and run fastqc

### Trimmomatic

This option runs trimmomatic on the raw sequencing reads. For the following, one needs trimmed reads, so set `trimmomatic=TRUE`. If trimmomatic has alrady been ran, one can set `trimmomatic=FALSE` but check that the file names are as follow: `$sample_trimmed.1P.fq.gz`

### Bowtie2

This option runs bowtie2 on the trimmed reads against `genome`. If `genome2` is specified, it also run bowtie2 on this genome independently. After bowtie2, the resulting sam files are turn into sorted bam files and merged together. Bedtools genomecov is then ran in order to get the coverage at each position of `genome` and `genome2`, if specified.

### Blastn

This option can be ran at the same time as bowtie2 if `parallele=TRUE`. In this case, `threads` will be equally shared between blastn and bowtie2.
This option firslty processes the trimmed fq.gz to obtain a fasta file. If `genome2` is specified, blastn is ran on the fasta file against this `genome2`. Then, only reads that mapped on `genome2` are mapped on `genome1`. If `genome2` is commented, all reads of the fasta file are mapped on `genome1`.

