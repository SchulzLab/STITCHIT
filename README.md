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
 
## How to run the linear regression for REM filtering?

An R script is provided to run a linear model to filter REMs generated with either the MDL formulation of STITCHIT, UnifiedPeaks or GeneHancer mode. The most basic command to execute the script *Two_level_learning.R* is

    Rscript Two_level_learning.R --dataDir <Data directory of regulatory elements> --outDir <Output directory that will be generated> --reponse <Name of the response variable> --cores <#CPUs>

The script will learn a separate model for each gene for that segments are provided in the directory specified by *dataDir*. In the folder specified by *outDir* filtered segments will be stored. Similarly a performance overview file will be generated.

Further parameters are:
	
        - outDir Output directory (will be created if it does not exist)
        - dataDir Directory containing the data
        - response Name of the response variable
        - cores Number of cores to be use (default 1)
        - fixedAlpha Use a fixed value for the alpha parameter in elastic net regulatisation, do not perform a grid search
        - alpha Stepsize to optimise the alpha parameter in elastic net regularisation (default 0.05)
        - testsize Size of test data[%] (default 0.2)
        - regularisation L for Lasso, R for Ridge, and E for Elastic net (default E)
        - innerCV Number of folds for inner cross-validation (default 6)
        - outerCV Number of iterations of outer cross-validation to determine test error (default 3)
        - constraint Specifies a constraint on the coefficent sign, enter N for negative and P for positive constraint
        - performance Flag indiciating whether the performance of the model should be assessed (default TRUE)
        - seed Random seed used for random number generation (default random)
        - leaveOneOutCV Flag indicating whether a leave one out cross-validation should be used (default FALSE)
        - asRData Store feature coefficients as RData files (default FALSE)
        - randomise Randomise the feature matrix (default FALSE)
        - logResponse Flag indicating whether the response variable should be log transformed (default TRUE)
        - ftest Flag indicating whether partial F-test should be computed to assess the significance of each feature (default FALSE)
        - coefP p-value threshold for model coefficient (default 1, all OLS coefs will be returned)

The parameters *cores*, *innerCV*, *outerCV*, *ftest* and *alpha* can considerably affect runtime.


## How to run the Nested version of STITCHIT?

As explained [here](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab798/6368526), a nested execution of STITCHIT enables out-of-sample performance estimates. As this nested execution is more computationally involved, we provide a python script to manage the workflow: *NestedSTITCHIT.py*.
Note that this script exclusively works with STICHIT, not UnifiedPeaks or GeneHancer.

The script needs to be executed for each gene of interest as:

    python NestedSTITCHIT.py <Normalised BigWig files> <Gene Annotation> <Discretised Expression Data> <Original expression data> <OutputFolder> <Size File> <geneID> <#Folds>  <Folder for R output>

For each gene STITCHIT will be executed <#Folds> times on a randomly selected set of samples (Monte Carlo Cross validation) and a linear model will be learned for each subset and evaluated on hold out data not used for feature generation. Note that using the nested execution is very resource and time-intensive.
	
	
	
## FAQ
Q: I can not compile STITCHIT on my Mac! What is wrong?

A: Make sure that the symbolic links of your compilers are set properly. We have experienced that these are not updated properly in Mac OS updates.


Q: Where can I get the data from the manuscript?

A: The data is available online at zenodo (10.5281/zenodo.2547383). You can also access predicted REM-gene associations and other information through the (Epiregio database) [https://epiregio.de]. Please cite (Baumgarten et al. 2020)[https://academic.oup.com/nar/article/48/W1/W193/5847772] for use this way.


Q: How should I cite STITCHIT?

A: Please cite our NAR article (Schmidt et al., 2021) [10.1093/nar/gkab798](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab798/6368526).
   



## Acknowledgements
This project has been funded by the Bundesministerium für Bildung und Forschung (BMBF) with project number 01DP17005 under the acronym EPIREG and the German Center for Cardiovascular Research (DZHK).
