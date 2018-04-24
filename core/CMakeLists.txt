project (STITCHCore CXX)

add_executable (STITCH main.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp SPANInputGenerator.cpp CorComp.cpp BinSelection.cpp)

add_executable (UNIFIED_PEAKS peak_Combination.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp CorComp.cpp BinSelection.cpp PeakReader.cpp)

target_link_libraries (STITCH BigWig ${Boost_LIBRARIES} SPAN)
target_link_libraries (UNIFIED_PEAKS BigWig ${Boost_LIBRARIES})
