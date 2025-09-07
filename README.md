<p align="center" >
    <img src="https://github.com/SAMtoBAM/fusemblr/blob/main/logo/fusemblr.png" width=100%>
</p>

[![Zenodo DOI](https://zenodo.org/badge/963417762.svg)](https://doi.org/10.5281/zenodo.15190276)
[![Anaconda_version](https://anaconda.org/samtobam/fusemblr/badges/version.svg)](https://anaconda.org/samtobam/fusemblr)
[![Anaconda_platforms](https://anaconda.org/samtobam/fusemblr/badges/platforms.svg)](https://anaconda.org/samtobam/fusemblr)
[![Anaconda_downloads](https://anaconda.org/samtobam/fusemblr/badges/downloads.svg)](https://anaconda.org/samtobam/fusemblr)
[![Anaconda-Server Badge](https://anaconda.org/samtobam/fusemblr/badges/latest_release_date.svg)](https://anaconda.org/samtobam/fusemblr)

**_fusemblr_** is a pipeline wrapper designed for the assembly of complex genomes using nanopore reads and paired-end illumina

**_fusemblr_**  was designed for the <i>Fusarium oxysporum</i> assembly project (hence the name) <br/>
The pipeline only requires Nanopore reads (the longer and higher coverage the better) and an estimation of genome size <br/>
Paired-end illumina reads and PacBio is optional <br/>

<i>Notably: Providing illumina and/or PacBio Hifi had very little impact on the resulting assemblies using our _Fusarium oxysporum_ datasets as we used recent ONT basecalled data, had high coverage and a good subset of long reads.</i>

# Easy installation

	conda install samtobam::fusemblr

# Container image

	docker pull ghcr.io/samtobam/fusemblr:latest

# How to run

 	fusemblr.sh -n nanopore.fq.gz -g 70000000
	
	Required inputs:
	-n | --nanopore		Nanopore long reads used for assembly in fastq or fasta format (*.fastq / *.fq) and can be gzipped (*.gz)
	-g | --genomesize	Estimation of genome size, required for downsampling and assembly

	Recommended inputs:
	-1 | --pair1		Paired end illumina reads in fastq format; first pair. Used for Ratatosk polishing. Can be gzipped (*.gz)
	-2 | --pair2		Paired end illumina reads in fastq format; second pair. Used for Ratatosk polishing. Can be gzipped (*.gz)		
 	-h | --hifi		Pacbio HiFi reads required for assembly polishing with NextPolish2 (Recommended if available)
	-t | --threads		Number of threads for tools that accept this option (default: 1)
	
	Optional parameters:
	-m | --minsize		Minimum size of reads to keep during downsampling (Default: 5000)
	-x | --coverage		The amount of coverage for downsampling (X), based on genome size, i.e. coverage*genomesize (Default: 100)
	-v | --minovl		Minimum overlap for Flye assembly,  (Default: Calculated during run as N95 of reads used for assembly)
 	-w | --weight		The weighting used by Filtlong for selecting reads; balancing the length vs the quality (Default: 5)
	-p | --prefix		Prefix for output (default: name of assembly file (-a) before the fasta suffix)
	-o | --output		Name of output folder for all results (default: fusemblr_output)
	-c | --cleanup		Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help		Print this help message



# Pipeline in 6 steps: <br/>
### 1. Downsampling of reads to a designated coverage using ```Filtlong```
##### &nbsp; &nbsp; -default is set to 100X (-x); which provided better assemblies compared to the typical 30-50X 
### 2. _Optional_: Polishing of downsampled reads with the paired-end illumina reads using ```Ratatosk correct``` 
##### &nbsp; &nbsp; -uses a baseline quality score (-Q) of 90 and therefore assumes mildly recent ONT data (e.g. R10 or high-accuracy basecalling)
### 3. Genome Assembly
#### 3.a. Assembly with```Flye``` 
##### &nbsp; &nbsp; -removed the hard coded maximium value for the minimum overlap threshold (previously 10kb) 
##### &nbsp; &nbsp; -by default the minimum overlap value is automatically provided as the read N90 after polishing
#### 3.b. Assembly with ```Hifiasm```
##### &nbsp; &nbsp; -if Hifi reads are provided: uses the ```--ul``` option, with both polished ONT and Hifi reads
##### &nbsp; &nbsp; -without Hifi: uses the ```--ont``` option, with only the polished ONT reads
### 4. 'Patch' the Flye assembly (target) using the the Hifiasm assembly (query) with ```Ragtag patch```
##### &nbsp; &nbsp; -uses a minimum unique alignment length (-f) of 25000 to be conservative during patching
### 5. _Optional_: Polishing of assembly with PacBio Hifi and paired-end illumina reads using ```NextPolish2```
### 6. Filtering (minimum contig length 10kb), reordering and renaming using ```Seqkit``` and ```awk```

## Schematic

<p align="center" >
    <img src="https://github.com/SAMtoBAM/fusemblr/blob/main/figures/fusemblr_schematic.png" width=70%>
</p>



Following assembly it is recommended that you run [PAQman](https://github.com/SAMtoBAM/PAQman) on your resulting assembly to comprehensively check the quality <br/>
It is recommended to feed your resulting assembly to PAQman alongside the 1.filtlong/*.fz.gz set of reads <br/>
This can also help you compare any assemblies you have to check for the best.


