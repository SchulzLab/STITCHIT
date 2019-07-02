project (STITCHCore CXX)

add_executable (STITCH main.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp SPANInputGenerator.cpp CorComp.cpp BinSelection.cpp)

add_executable (UNIFIED_PEAKS peak_Combination.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp CorComp.cpp BinSelection.cpp PeakReader.cpp)

add_executable (GENEHANCER geneHancer_Annotation.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp CorComp.cpp BinSelection.cpp GeneHancerReader.cpp)

add_executable (REMSELECT coverage_Computation_from_bw.cpp CorComp.cpp BinSelection.cpp RemReader.cpp)

target_link_libraries (STITCH BigWig ${Boost_LIBRARIES} SPAN)
target_link_libraries (UNIFIED_PEAKS BigWig ${Boost_LIBRARIES})
target_link_libraries (GENEHANCER BigWig ${Boost_LIBRARIES})
target_link_libraries (REMSELECT BigWig ${Boost_LIBRARIES})
