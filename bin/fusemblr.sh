#!/bin/bash
set -euo pipefail

version="v1"

##A wrapper for the fusarium assembly project = fusemblr (not the best name but hey, it was almost fusembly, like fusilli; the worst pasta invented)

##this pipeline has 5 main steps:
## STEP 1: Downsampling on raw ONT reads (Filtlong)
## STEP 2: Polishing of raw ONT reads with paired end illumina data (Ratatosk)
## STEP 3: Assembly of polished ONT reads (Flye; modified to allow for larger minimum overlap values)
## Step 4: Polishing of assembly using Pacbio (NextPolish2; optional)
## Step 5: Filtering, reordering and renaming (Seqkit)

#####################################################################################
############# STEP -1. CREATING THE ENVIRONMENT. NOT ACTUALLY USED NOW ##############
#####################################################################################

##environment 
##installing
## Ratatosk (raw read polisher)
## filtlong (kmer based long-read downsampling/filtering tool)
## seqkit (many uses)
## flye (assembler)
## nextpolish2 (long-read polisher for use with pacbio hifi is available)
## fastp (quality control of illumina data used for nextpolish2)
## minimap2 (read alignment for nextpolish2)
## samtools (read alignment for nextpolish2)

# mamba create -n fusemblr ratatosk bioconda::filtlong bioconda::flye bioconda::fastp nextpolish2 bioconda::seqkit bioconda::minimap2 bioconda::samtools
#conda activate fusemblr

##modify the maximum value for minoverlap in flye (making 200kb..the N95 of a read dataset shouldn't exceed that)
#sed -i 's/v, 1000, 10000)/v, 1000, 200000)/g' $CONDA_PREFIX/lib/python3.*/site-packages/flye/main.py


##############################################################
################ STEP 0a: SETTING VARIABLES ##################
##############################################################

#default values, unless denoted when running MUM&Co
nanopore=""
pair1=""
pair2=""
genomesize=""
hifi=""
threads="1"
minsize="5000"
coverage="100"
weight="5"
minovl=""
prefix=""
output="fusemblr_output"
help="nohelp"

## to clean up a bunch of output from the tools in order to reduce all the unnecessary output
cleanup="yes"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
	-n|--nanopore)
	nanopore="$2"
	shift
	shift
	;;
	-1|--pair1)
	pair1="$2"
	shift
	shift
	;;
	-2|--pair2)
	pair2="$2"
	shift
	shift
	;;
	-g|--genomesize)
	genomesize="$2"
	shift
	shift
	;;
	-h|--hifi)
	hifi="$2"
	shift
	shift
	;;
	-t|--threads)
	threads="$2"
	shift
	shift
	;;
	-m|--minsize)
	minsize="$2"
	shift
	shift
	;;	
	-x|--coverage)
	coverage="$2"
	shift
	shift
	;;
 	-w|--weight)
	weight="$2"
	shift
	shift
	;;
	-v|--minovl)
	minovl="$2"
	shift
	shift
	;;
	-p|--prefix)
	prefix="$2"
	shift
	shift
	;;
	-o|--output)
	output="$2"
	shift
	shift
	;;
	-c|--cleanup)
	cleanup="$2"
	shift
	shift
	;;
	-h|--help)
	echo "
	
	fusemblr (version: ${version})
 
	fusemblr.sh -n nanopore.fq.gz -1 illumina.R1.fq.gz -2 illumina.R2.fq.gz -g 70000000
	
	Required inputs:
	-n | --nanopore		Nanopore long reads used for assembly in fastq or fasta format (*.fastq / *.fq) and can be gzipped (*.gz)
	-1 | --pair1		Paired end illumina reads in fastq format; first pair. Used for Rataosk polishing. Can be gzipped (*.gz)
	-2 | --pair2		Paired end illumina reads in fastq format; second pair. Used for Rataosk polishing. Can be gzipped (*.gz)	
	-g | --genomesize	Estimation of genome size, required for downsampling and assembly

	Recommended inputs:
	-h | --hifi		Pacbio HiFi reads required for assembly polishing with NextPolish2 (Recommended if available)
	-t | --threads		Number of threads for tools that accept this option (default: 1)
	
	Optional parameters:
	-m | --minsize		Minimum size of reads to keep during downsampling (Default: 5000)
	-x | --coverage		The amount of coverage for downsampling (X), based on genome size, i.e. coverage*genomesize (Default: 100)
	-v | --minovl		Minimum overlap for Flye assembly (Default: Calculated during run as N95 of reads used for assembly)
 	-w | --weight		The weighting used by Filtlong for selecting reads; balancing the length vs the quality (Default: 5)
	-p | --prefix		Prefix for output (Default: name of nanopore reads file (-a) before the fastq suffix)
	-o | --output		Name of output folder for all results (Default: fusemblr_output)
	-c | --cleanup		Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help		Print this help message

	"
	exit
	;;
	esac
done


#creates error message and exits if these values are not/incorrectly assigned 
[[ $nanopore == "" ]] && echo "ERROR: Path to ONT long-reads not found, assign using -n" && exit
[[ $pair1 == "" || $pair2 == ""  ]] && echo "ERROR: Missing paired-end short-reads as input, assign using -1 and -2" && exit
[[ $genomesize == "" ]] && echo "ERROR: Missing genome size estimation, assign using -g (can be a very rough estimate)" && exit

[[ $hifi == "" ]] && echo "WARNING: Running without PacBio Hifi data; skipping final polishing step"

##############################################################
####################### 0b. SETTING UP #######################
##############################################################


##redefining some variables and checking
##if prefix wasn't defined, grab the prefix of the nanopore reads
[[ $prefix == "" ]] && prefix=$( echo ${nanopore} | awk -F "/" '{print $NF}' | sed 's/\.fastq\.gz$//' | sed 's/\.fq\.gz$//' | sed 's/\.fastq$//' | sed 's/\.fq$//' )

##paths to raw data
nanoporepath=$( realpath ${nanopore} )
pair1path=$( realpath ${pair1} )
pair2path=$( realpath ${pair2} )

#check if the files given actually exist
[ ! -f "${nanoporepath}" ] && echo "ERROR: Cannot find path to nanopore reads provided by -n; check path is correct and file exists" && exit
[ ! -f "${pair1path}" ] && echo "ERROR: Cannot find path to illumina reads provided by -1; check path is correct and file exists" && exit
[ ! -f "${pair2path}" ] && echo "ERROR: Cannot find path to illumina reads provided by -2; check path is correct and file exists" && exit

if [[  $hifi != "" ]]
then
hifipath=$( realpath ${hifi} )
[ ! -f "${hifipath}" ] && echo "ERROR: Cannot find path to Pacbio Hifi reads provided by -h; check path is correct and file exists" && exit
fi


##newly defining some variables for output writing
## leave as is just setting to kb
minsize2=$( echo ${minsize} | awk '{print $1/1000}' )
## number of bases desired by filtlong (think 70000000*coverage)
target=$( echo $genomesize | awk -v coverage="$coverage" '{print $1*coverage}' )
## set a variable string of all the filtlong settings to track the read dataset being used
readstats=$( echo "min${minsize2}kb_${coverage}X_weightlen${weight}" )
## converts size of genome to Mb
size2=$( echo $genomesize | awk '{print $1/1000000}' )

[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

mkdir ${output}
cd ${output}

##################################################################
###################### STEP 1. DOWNSAMPLING ######################
##################################################################

echo "#################################################################"
echo "################## fusemblr: Starting fusemblr ##################"
echo "#################################################################"
echo "################## fusemblr: Step 1: Downsampling ONT reads"


##make a directory for placing the downsampled output
mkdir 1.filtlong_ont

##now run filtlong with the settings
filtlong --min_length ${minsize} -t ${target} --length_weight ${weight} ${nanoporepath} | gzip > 1.filtlong_ont/${prefix}.${readstats}.fq.gz


###########################################################################
###################### 2. READ POLISHING WITH RATAOSK #####################
###########################################################################

echo "################## fusemblr: Step 2: Polishing ONT reads"


## create directory for output of reads
mkdir 2.ratatosk_ont
## increase base quality score minimum to 90 due to high quality reads (-Q)
Ratatosk correct -v -Q 90 -c ${threads} -G -s $( ls ${pair1path} ${pair2path} ) -l 1.filtlong_ont/${prefix}.${readstats}.fq.gz -o 2.ratatosk_ont/${prefix}.${readstats}.ratatosk > ratatosk.fusemblr.log
##modify output file name
mv 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.fastq.gz 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.fq.gz
##move log file into folder
mv ratatosk.fusemblr.log 2.ratatosk_ont/

##get some stats
seqkit stats -N 50,90,95 --threads ${threads} 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.fq.gz > 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.stats.tsv

###get the read N90 to set as a variables in flye
if [[ $minovl == "" ]]
then
minovl=$( tail -n1 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.stats.tsv | awk '{print $11}' | sed 's/,//g' )
minovl2=$( echo ${minovl} | awk '{print $1/1000}' | awk -F "." '{print $1}' )
else
minovl2=$( echo ${minovl} | awk '{print $1/1000}' | awk -F "." '{print $1}' )

fi

#################################################################
################## STEP 3a. ASSEMBLY WITH FLYE ##################
#################################################################

echo "################## fusemblr: Step 3a: Assembling ONT reads (Flye)"


## create a variable that saves these variable in the naming scheme
assembly="${prefix}.${readstats}.ratatosk.flye_nanocorr_${size2}Mb_${coverage}X_minovl${minovl2}k"

## run flye
flye --nano-corr 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.fq.gz -m ${minovl} --genome-size ${genomesize} --asm-coverage ${coverage} --threads ${threads} -o 3a.flye_assembly/ > flye.fusemblr.log
## make sure the assembly file is named is a simple format
cp 3a.flye_assembly/assembly.fasta 3a.flye_assembly/${assembly}.fa
##move log file into folder
mv flye.fusemblr.log 3a.flye_assembly/flye.fusemblr.log

#################################################################
################ STEP 3b. ASSEMBLY WITH HIFIASM #################
#################################################################

echo "################## fusemblr: Step 3b: Assembling ONT reads (Hifiasm)"


## create a variable that saves these variable in the naming scheme
assembly="${prefix}.${readstats}.ratatosk.hifiasm"


if [[  $hifi != "" ]]
then
## run hifiasm
hifiasm -o 3b.hifiasm/${prefix} -t ${threads} -l0 --ont 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.fq.gz
else
hifiasm -o 3b.hifiasm/${prefix} -t ${threads} -l0 --ul 2.ratatosk_ont/${prefix}.${readstats}.ratatosk.fq.gz ${hifipath}
fi
## convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' 3b.hifiasm/${prefix}.bp.p_ctg.gfa > 3b.hifiasm/${prefix}.hifiasm.fa


##################################################################
############ STEP 4. ASSEMBLY POLISHING WITH PB HIFI  ############
##################################################################

##flye polishing
if [[  $hifi != "" ]]
then

echo "################## fusemblr: Step 4: Polishing assembly with Hifi"

mkdir 4a.flye_assembly.nextpolish2/

## prepare long-read alignments to assembly
minimap2 -ax map-hifi -t ${threads} 3a.flye_assembly/${assembly}.fa ${hifipath} | samtools sort -@ 4 -o 4a.flye_assembly.nextpolish2/minimap_pacbio.sort.bam -
samtools index 4a.flye_assembly.nextpolish2/minimap_pacbio.sort.bam

## prepare illumina data
fastp -5 -3 -n 0 -f 5 -F 5 -t 5 -T 5 -q 20 -i ${pair1path} -I ${pair2path} -o ${prefix}.R1.clean.fq.gz -O ${prefix}.R2.clean.fq.gz
## produce a 21 and 31-mer dataset (needs a lot of memory)
yak count -o k21.yak -k 21 -b 37 <(zcat ${prefix}.R*.clean.fq.gz) <(zcat ${prefix}.R*.clean.fq.gz)
yak count -o k31.yak -k 31 -b 37 <(zcat ${prefix}.R*.clean.fq.gz) <(zcat ${prefix}.R*.clean.fq.gz)

## now run Nextpolish with the inputs generated above
nextPolish2 -t ${threads} 4a.flye_assembly.nextpolish2/minimap_pacbio.sort.bam 3a.flye_assembly/${assembly}.fa k21.yak k31.yak > 4a.flye_assembly.nextpolish2/${assembly}.nextpolish2.fa

## remove intermediate (alignment/yak) files that take up a lot of space
rm 4a.flye_assembly.nextpolish2/minimap_pacbio.sort.*

##convert the final assembly to a simple name
cp 4a.flye_assembly.nextpolish2/${assembly}.nextpolish2.fa ${prefix}.flye.fa

else

##convert the final assembly to a simple name
cp 3a.flye_assembly/${assembly}.fa ${prefix}.flye.fa

fi


##hifiasm polishing
if [[  $hifi != "" ]]
then

echo "################## fusemblr: Step 4: Polishing assembly with Hifi"

mkdir 4b.hifiasm_assembly.nextpolish2/

## prepare long-read alignments to assembly
minimap2 -ax map-hifi -t ${threads} 3b.hifiasm/${prefix}.hifiasm.fa ${hifipath} | samtools sort -@ 4 -o 4b.hifiasm_assembly.nextpolish2/minimap_pacbio.sort.bam -
samtools index 4b.hifiasm_assembly.nextpolish2/minimap_pacbio.sort.bam

## now run Nextpolish with the inputs generated above
nextPolish2 -t ${threads} 4b.hifiasm_assembly.nextpolish2/minimap_pacbio.sort.bam 3a.flye_assembly/${assembly}.fa 4b.hifiasm_assembly.nextpolish2/k21.yak 4b.hifiasm_assembly.nextpolish2/k31.yak > 4b.hifiasm_assembly.nextpolish2/${prefix}.hifiasm.nextpolish2.fa

## remove intermediate (alignment/yak) files that take up a lot of space
rm 4b.hifiasm_assembly.nextpolish2/minimap_pacbio.sort.*
rm *.yak
rm ${prefix}.R*.clean.fq.gz
rm fastp.html
rm fastp.json

##convert the final assembly to a simple name
cp 4b.hifiasm_assembly.nextpolish2/${prefix}.hifiasm.nextpolish2.fa ${prefix}.hifiasm.fa

else

##convert the final assembly to a simple name
cp 3b.hifiasm/${prefix}.hifiasm.fa ${prefix}.hifiasm.fa

fi

#####################################################################
############ STEP 5. ASSEMBLY GAP FILLING  ############
#####################################################################

echo "################## fusemblr: Step 5: Filling gaps in Flye assembly using Hifiasm"

ragtag.py patch -o 5.ragtag_patch -f 25000 ${prefix}.flye.fa ${prefix}.hifiasm.fa
cp 5.ragtag_patch/ragtag.patch.fasta ${prefix}.flye.rt_patch.fa


#####################################################################
############ STEP 6. FILTERING, REORDERING AND RENAMING  ############
#####################################################################


echo "################## fusemblr: Step 6: Filtering, ordering and renaming"

##filter out any sequences smaller than 10kb, sort by length and then rename as numbered contig in order of largest to smallest (1 being the largest)
seqkit seq -m 10000 ${prefix}.flye.rt_patch.fa | seqkit sort -l -r - | awk 'BEGIN{n=1} {if($1 ~ ">") {print ">contig_"n; n++} else{print}}'  > ${prefix}.fa

#####################################################################
############################# FINISHED  #############################
#####################################################################

echo "################## fusemblr: Thanks for using fusemblr"
