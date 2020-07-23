# STITCHIT
Linking regulatory elements to genes. STITCHIT uses a two-level paradigm to first define regulatory elements for a gene and then estimate its relevance. Thus, STITCHIT presents a one-stop solution to derive regions with potential gene regulatory 
role.

The first step is a supervised segmentation method, that uses epigenomics data around a gene and the expression of that gene in many samples as information to derive de novo gene regulatory regions.
The second step then uses these derived regions and performs a sparse regression using the epigenomics data (signal) and gene expression (response) to output regions of importance and assess their weights and significance of each gene regulatory region.

Example applications of STITCHIT are:

1. You have measured ATAC-seq and RNA-seq on a number of patient samples for a condition of interest. You want to understand which regulatory regions (such as enhancers) exist and influence some genes of interest. For example for overlapping with GWAS or eQTL data.
2. You want to study gene regulation by integrating many samples of paired epigenomics (DNAse1, histone ChIP-seq) and expression data, for example using data from ENCODE or TCGA etc. 
3. You are interested in studying transcription factor (TF) binding and want to get a better estimate of gene-specific regions where TFs may be interacting with the DNA sequence to regulate genes of interest.



## Building the tool
To build the **STITCHIT** *cmake*, a C++11 compiler, and the boost library must be available on your system.
We have tested STITCHIT with Linux operating systems. 

Make sure you are in the projects main directory. Then, generate a build folder by typing:
`mkdir build`

Change into that folder:
`cd build`

To run *cmake* use:
`cmake ..`

To finally build the project use:
`make`

You can speed up the build process by using the *-j* option, specifying how many cores should be used for building

To execute all tests use:
`make test`

Tests can also be executed individually. They are located in the folder
`build/test`

## Input to STITCHIT
STITCHIT needs three sources of data as input:

1. Epigenomics data (such as DNAse1-seq, ATAC-seq, FAIRE-seq or histone ChIP-seq) in form of wiggle files.
2. Gene expression table (such as RNA-seq or microarrays). The expression data should in addition be given in discretized form as this is used for the first step in STITCHIT.
3. Gene annotation in GTF format.

Data in 1. and 2. needs to be measured jointly on many samples (ideally 30 or more).

## How to run STITCHIT?
You can call STITCHIT it using the following command:

	./build/core/STITCHIT –b <DNase_bw> -a <gencode.v26.annotation.gtf> -d <Discretised_Expression> -o <Continuous_Expression>  -s <Chromosome_Size> -w <ExtensionSize> -c <cores> -p <p-value> -g <geneID> -z <Initial merge> -f <Outputpath> -r <maximum count side> -t <segmentsize>

I order to run STITCHIT with the data a number of parameters need to be set:


  - b: Path to big wig files holding epigenetic signal that should be used for the segmentation
  - a: Genome annotation file used to find the genomic position of the target gene
  - d: Discretised expression data
  - o: Continuous expression data
  - s: File holding the size of the individual chromosomes of the target organism
  - w: Extension up and downstream of the target gene
  - c: Number of cores used within STITCHIT
  - p: P-value threshold used to select a regulatory element
  - g: Target gene ID
  - z: Resolution used to merge the initial data to cut runtime (default 10)
  - f: Output path
  - r: Maximum size of the entire considered search space
  - t: Size of the segments used within STITCHIT (default 2000)
 

## Running alternative approaches
In our paper we compare with some alternative ways of estimating gene-specific regulatory elements. How to run these approaches is detailed below.

## How to run the Unified peaks approach?

In addition to learning novel regulatory elements, we provide the possibility to give previously defined regions as input to the second step of STITCHIT (the sparse regression step). Then the first step of defining the regions from the epigenomics signal is skipped. Using the data provided to STITCHIT it is decided which of these regions are associated with gene expression. 
	
	./build/core/UNIFIED_PEAKS -k <Merged peak file>

The new parameter here is *-k*, which provides the merged peak regions in bed file format <chr> <start> <end>. 
The file should be sorted.

## How to run GeneHancer?

Similarly to above we provide the possibility to give previously defined regions as input to the second step of STITCHIT (the sparse regression step), but without evaluation which of these regions are of interest for the gene. That means we are not filtering these regions before the regression. This mode makes sense, when you have obtained the gene-links of regulatory regions from a database such as GeneHancer.

	./build/core/GENEHANCER -k <GeneHancer File>

The new parameter here is *-k*, which provides the GeneHancer regions in bed file format <chr> <start> <end> <ENSGID>. 
The file must be sorted according to the geneIDs.
 
## FAQ
Q: I can not compile STITCHIT on my Mac! What is wrong?

A: Make sure that the symbolic links of your compilers are set properly. We have experienced that these are not updated properly in Mac OS updates.


Q: Where can I get the data from the manuscript?

A: The data is available online at zenodo (10.5281/zenodo.2547383)


Q: How should I cite STITCHIT?

A: You can cite our bioRxiv article: [doi.org/10.1101/585125](https://www.biorxiv.org/content/10.1101/585125v1.full)


## Acknowledgements
This project is being funded by the Bundesministerium für Bildung und Forschung (BMBF) with project number 01DP17005 under the acronym EPIREG.
