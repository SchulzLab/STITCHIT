project (STITCHCore CXX)

add_executable (STITCH main.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp SPANInputGenerator.cpp CorComp.cpp BinSelection.cpp)

target_link_libraries (STITCH BigWig ${Boost_LIBRARIES} SPAN)
