project(LinkingEnhancerToGenesBigWig C)

add_library(BigWig
	io.c
	bwValues.c
	bwRead.c
	bwStats.c
	bwWrite.c)

target_link_libraries(BigWig -lcurl -lm -lz)
