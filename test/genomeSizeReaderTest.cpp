#include "gtest/gtest.h"
#include "../core/GenomeSizeReader.cpp"
#include "config.h"

class GenomeSizeReaderTest : public ::testing::Test{
	protected:
	GenomeSizeReaderTest()
		:g_(GenomeSizeReader(TEST_DATA_PATH("hg38_chrSize.txt"))),
		 g2_(GenomeSizeReader(TEST_DATA_PATH("hg19_chrSize.txt"))),
		 g3_(GenomeSizeReader(TEST_DATA_PATH("Phony_chrSize.txt")))
	{
	}

	void virtual SetUp(){
	
	}
	public:
	GenomeSizeReader g_;
	GenomeSizeReader g2_;
	GenomeSizeReader g3_;

};

TEST_F(GenomeSizeReaderTest, Constructor){
	SUCCEED();
}

TEST_F(GenomeSizeReaderTest, getfileName){
	ASSERT_EQ(TEST_DATA_PATH("hg38_chrSize.txt"),g_.getFilename());
}

TEST_F(GenomeSizeReaderTest, loadGenomeSizeFile){
	g_.loadGenomeSizeFile();
	SUCCEED();
}

TEST_F(GenomeSizeReaderTest, GenomeSizeFileDoesNotExist){
	ASSERT_THROW(g2_.loadGenomeSizeFile(),std::invalid_argument);
}

TEST_F(GenomeSizeReaderTest, GenomeSizeFileIsPhony){
	ASSERT_THROW(g3_.loadGenomeSizeFile(),std::invalid_argument);
}

TEST_F(GenomeSizeReaderTest, IsEmpty){
	ASSERT_TRUE(g_.getGenomeSize().empty());
}

TEST_F(GenomeSizeReaderTest, IsNotEmpty){
	g_.loadGenomeSizeFile();
	ASSERT_FALSE(g_.getGenomeSize().empty());
}

TEST_F(GenomeSizeReaderTest, chr1){
	g_.loadGenomeSizeFile();
	std::map<std::string, int> genomeSize = g_.getGenomeSize();
	ASSERT_EQ(genomeSize["chr1"],248956422);
}

TEST_F(GenomeSizeReaderTest, chrY){
	g_.loadGenomeSizeFile();
	std::map<std::string, int> genomeSize = g_.getGenomeSize();
	ASSERT_EQ(genomeSize["chrY"],57227415);
}

TEST_F(GenomeSizeReaderTest, chr21){
	g_.loadGenomeSizeFile();
	std::map<std::string, int> genomeSize = g_.getGenomeSize();
	ASSERT_EQ(genomeSize["chr21"],46709983);
}
