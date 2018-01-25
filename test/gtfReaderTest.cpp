#include "gtest/gtest.h"
#include "../core/GTFReader.cpp"
#include "../core/GenomeSizeReader.cpp"

#include "config.h"

class GTFReaderTest : public ::testing::Test{
	protected:
	GTFReaderTest()
		:g_(GenomeSizeReader(TEST_DATA_PATH("hg38_chrSize.txt")))
	{
		g_.loadGenomeSizeFile();

	}

	void virtual SetUp(){
	
	}

	public:
	GenomeSizeReader g_;
};

TEST_F(GTFReaderTest, Constructor){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	SUCCEED();
}

TEST_F(GTFReaderTest, FileName){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	ASSERT_EQ(TEST_DATA_PATH("gencode.sampleV26.gtf"),gtf_.getGTFfileName());
}

TEST_F(GTFReaderTest, FileDoesNotExist){
	GTFReader gtf_(TEST_DATA_PATH("gencode.v612.annotation.gtf"),g_.getGenomeSize(),5000);
	ASSERT_THROW(gtf_.findGenomicLocation("ENSG00000184990"),std::invalid_argument);
}

TEST_F(GTFReaderTest, PhonyFile1){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.Phony1.gtf"),g_.getGenomeSize(),5000);
	ASSERT_THROW(gtf_.findGenomicLocation("ENSG00000184990"),std::invalid_argument);
}

TEST_F(GTFReaderTest, PhonyFile2){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.Phony2.gtf"),g_.getGenomeSize(),5000);
	ASSERT_THROW(gtf_.findGenomicLocation("ENSG00000184990"),std::invalid_argument);
}

TEST_F(GTFReaderTest, GetWindow){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	ASSERT_EQ(gtf_.getWindow(),5000);
}

TEST_F(GTFReaderTest, SetWindow){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	gtf_.setWindow(0);
	ASSERT_EQ(gtf_.getWindow(),0);
}

TEST_F(GTFReaderTest, EmptyRegion){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	std::tuple<std::string,unsigned int, unsigned int,std::string> location = gtf_.getGenomicLocation();
	ASSERT_EQ(std::get<0>(location),"chr0");
	ASSERT_EQ(std::get<1>(location),0);
	ASSERT_EQ(std::get<2>(location),0);
}

TEST_F(GTFReaderTest, StandardGeneRetrievalNoShift){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),0);
	gtf_.findGenomicLocation("ENSG00000184990");
	std::tuple<std::string,unsigned int, unsigned int,std::string> location = gtf_.getGenomicLocation();
	ASSERT_EQ(std::get<0>(location),"chr14");
	ASSERT_EQ(std::get<1>(location),104753100);
	ASSERT_EQ(std::get<2>(location),104768494);
}

TEST_F(GTFReaderTest, StandardGeneRetrievalShift){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	gtf_.findGenomicLocation("ENSG00000184990");
	std::tuple<std::string,unsigned int, unsigned int,std::string> location = gtf_.getGenomicLocation();
	ASSERT_EQ(std::get<0>(location),"chr14");
	ASSERT_EQ(std::get<1>(location),104748100);
	ASSERT_EQ(std::get<2>(location),104773494);
}                                         
TEST_F(GTFReaderTest, StandardGeneRetrievalNegativePosition){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),104753110);
	gtf_.findGenomicLocation("ENSG00000184990");
	std::tuple<std::string,unsigned int, unsigned int,std::string> location = gtf_.getGenomicLocation();
	ASSERT_EQ(std::get<0>(location),"chr14");
	ASSERT_EQ(std::get<1>(location),0);
}

TEST_F(GTFReaderTest, StandardGeneRetrievalExceedingGenomeSizePosition){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),2500000);
	gtf_.findGenomicLocation("ENSG00000184990");
	std::tuple<std::string,unsigned int, unsigned int,std::string> location = gtf_.getGenomicLocation();
	ASSERT_EQ(std::get<0>(location),"chr14");
	ASSERT_EQ(std::get<2>(location),107043718);
	ASSERT_EQ(std::get<3>(location),"+");
}

TEST_F(GTFReaderTest, GeneDoesNotExistInGTF){
	GTFReader gtf_(TEST_DATA_PATH("gencode.sampleV26.gtf"),g_.getGenomeSize(),5000);
	ASSERT_THROW(gtf_.findGenomicLocation("ENSG11111111111"),std::invalid_argument);
}
