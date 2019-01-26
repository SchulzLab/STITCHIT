# STITCHIT
Linking regulatory elements to genes

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

## How to run STITCHIT?
You can call STITCHIT it using the following command:

	./build/core/STITCHIT –b <DNase_bw> -a <gencode.v26.annotation.gtf> -d <Discretised_Expression> -o <Continuous_Expression>  -s <Chromosome_Size> -w <ExtensionSize> -c <cores> -p <p-value> -g <geneID> -z <Initial merge> -f <Outputpath> -r <maximum count side> -t <segmentsize>

The parameters are:


  - b: Path to big wig files holding epigenetic signal that should be used for the segmentation
  - a: Genome annotation file used to find the genomic position of the target gene
  - d: Discretised expression data
  - o: Continuous expression data
  - s: File holding the size of the individual chromosomes of the target organism
  - w: Extension up and downstream of the target gene
  - c: Number of cores used within STITCHIT
  - p: P-value threshold used to select a regulatory element
  - g: Target gene ID
  - z: Resolution used to merge the initial data to cut runtime (we suggest 10)
  - f: Output path
  - r: Maximum size of the entire considered search space
  - t: Size of the segments used within STITCHIT (we suggest 2000)
 


## How to run Unified peaks?
	
	./build/core/UNIFIED_PEAKS -k <Merged peak file>

The new parameter here is *-k*, which provides the merged peak regions in bed file format <chr> <start> <end>. 
The file should be sorted.

## How to run GeneHancer?

	./build/core/GENEHANCER -k <GeneHancer File>

The new parameter here is *-k*, which provides the GeneHancer regions in bed file format <chr> <start> <end> <ENSGID>. 
The file must be sorted according to the geneIDs.
 
## FAQ
Q: I can not compile STITCHIT on my Mac! What is wrong?
A: Make sure that the symbolic links of your compilers are set properly. We have experienced that these are not updated properly in Mac OS updates.

## Acknowledgements
This project is being funded by the Bundesministerium für Bildung und Forschung (BMBF) with project number 01DP17005 under the acronym EPIREG.
