# SPAN
Binning Algorithm to group features based on binary classification.

# Execute
1) use make to compile
2) add '-h' flag to get more information on how to run the code

# Example
Sample data is provided in the data folder. The id files are not used. A sample command to rund the code would be:

./span.run -i ../data/Integrated_featureMatrix_ENSG00000107581_5000.tab -c -p 10 -o ../test_results.tab

"-i <path>" input file
"-o <path>" output file
"-c" use class labels
"-p 10" split data into 10 parts
