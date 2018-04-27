#include "gtest/gtest.h"
#include "../core/GeneHancerReader.cpp"
#include "config.h"

class GeneHancerReaderTest : public ::testing::Test{
	protected:
	GeneHancerReaderTest()
		:g_(GeneHancerReader(TEST_DATA_PATH("geneHancerSample.bed")))
	{
	}

	void virtual SetUp(){
	
	}
	public:
	GeneHancerReader g_;

};

TEST_F(GeneHancerReaderTest, Constructor){
	SUCCEED();
}

TEST_F(GeneHancerReaderTest, Initialization_Vector){
	ASSERT_EQ(0,g_.getEnhancers().size());
}

TEST_F(GeneHancerReaderTest, Initialization_Filename){
	ASSERT_EQ(TEST_DATA_PATH("geneHancerSample.bed"),g_.getFilename());
}

TEST_F(GeneHancerReaderTest, WrongGeneIDFormat){
	std::string geneID="14341";
	ASSERT_THROW(g_.findEnhancers(geneID),std::invalid_argument);
	ASSERT_EQ(0,g_.getEnhancers().size());
}

TEST_F(GeneHancerReaderTest, GeneIDNotPresent){
	std::string geneID="ENSG000344341";
	ASSERT_THROW(g_.findEnhancers(geneID),std::invalid_argument);
	ASSERT_EQ(0,g_.getEnhancers().size());
}

TEST_F(GeneHancerReaderTest, PresentOnChromosome2C){
	std::string geneID="ENSG00000000005";
	g_.findEnhancers(geneID);
	ASSERT_EQ(1,g_.getEnhancers().size());
}

TEST_F(GeneHancerReaderTest, PresentOnChromosome2R){
	std::string geneID="ENSG00000000419";
	g_.findEnhancers(geneID);
	auto peaks = g_.getEnhancers();
	ASSERT_EQ(3,peaks.size());

	ASSERT_EQ("chr20",std::get<0>(peaks[0]));
	ASSERT_EQ(49981884,std::get<1>(peaks[0]));
	ASSERT_EQ(49985177,std::get<2>(peaks[0]));

	ASSERT_EQ("chr20",std::get<0>(peaks[1]));
	ASSERT_EQ(49993491,std::get<1>(peaks[1]));
	ASSERT_EQ(49993751,std::get<2>(peaks[1]));

	ASSERT_EQ("chr20",std::get<0>(peaks[2]));
	ASSERT_EQ(50290893,std::get<1>(peaks[2]));
	ASSERT_EQ(50298200,std::get<2>(peaks[2]));
}
