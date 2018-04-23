#include "gtest/gtest.h"
#include "../core/PeakReader.cpp"
#include "config.h"

class PeakReaderTest : public ::testing::Test{
	protected:
	PeakReaderTest()
		:g_(PeakReader(TEST_DATA_PATH("peakFile1.txt"))),
		 g2_(PeakReader(TEST_DATA_PATH("broken_peakFile.txt")))
	{
	}

	void virtual SetUp(){
	
	}
	public:
	PeakReader g_;
	PeakReader g2_;

};

TEST_F(PeakReaderTest, Constructor){
	SUCCEED();
}

TEST_F(PeakReaderTest, Initialization_Vector){
	ASSERT_EQ(0,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, Initialization_Filename){
	ASSERT_EQ(TEST_DATA_PATH("peakFile1.txt"),g_.getFilename());
}

TEST_F(PeakReaderTest, FileBroken){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",30,300);	
	ASSERT_THROW(g2_.findOverlappingPeaks(coordinates),std::invalid_argument);
}

TEST_F(PeakReaderTest, NotPresent){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr5",30,300);	
	ASSERT_THROW(g_.findOverlappingPeaks(coordinates),std::invalid_argument);
	ASSERT_EQ(0,g_.getOverlappingPeaks().size());
}
