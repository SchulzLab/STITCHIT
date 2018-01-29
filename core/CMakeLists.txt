project (STITCHCore CXX)

add_executable (LEG main.cpp GenomeSizeReader.cpp ExpressionReader.cpp GTFReader.cpp SPANInputGenerator.cpp CorComp.cpp BinSelection.cpp)

target_link_libraries (LEG BigWig ${Boost_LIBRARIES} SPAN)
