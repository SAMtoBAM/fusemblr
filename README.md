# fusembly
fusembly is a pipeline wrapper designed for the assembly of complex genomes using nanopore reads and paired-end illumina

Fusembly was designed for the Fusarium consortium assembly project (hence the name) in order to improve the automated steps during genome assembly




 
	fusembly -n nanopore.fq.gz -h hifi.fq.gz -1 illumina.R1.fq.gz -2 illumina.R2.fq.gz
	
	Required inputs:
	-n | --nanopore		Nanopore long reads used for assembly in fastq or fasta format (*.fastq / *.fq) and can be gzipped (*.gz)
	-1 | --pair1		Paired end illumina reads in fastq format; first pair. Used for Rataosk polishing and PAQman evaluation. Can be gzipped (*.gz)
	-2 | --pair2		Paired end illumina reads in fastq format; second pair. Used for Rataosk polishing and PAQman evaluation. Can be gzipped (*.gz)	
	-g | --genomesize	Estimation of genome size, required for downsampling and assembly

	Recommended inputs:
	-h | --hifi		Pacbio HiFi reads required for assembly polishing with NextPolish2
	-t | --threads		Number of threads for tools that accept this option (default: 1)
	
	Optional parameters:
	-m | --minsize		Minimum size of reads to keep during downsampling (Default: 5000)
	-x | --coverage		The amount of coverage for downsampling (X), based on genome size, i.e. coverage*genomesize (Default: 100)
	-v | --minovl		Minimum overlap for Flye assembly,  (Default: Calculated during run as N95 of reads used for assembly)
	-p | --prefix		Prefix for output (default: name of assembly file (-a) before the fasta suffix)
	-o | --output		Name of output folder for all results (default: fusembly_output)
	-c | --cleanup		Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help		Print this help message
