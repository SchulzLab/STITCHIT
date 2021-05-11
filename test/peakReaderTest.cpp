#include "gtest/gtest.h"
#include "../core/PeakReader.cpp"
#include "config.h"

class PeakReaderTest : public ::testing::Test{
	protected:
	PeakReaderTest()
		:g_(PeakReader(TEST_DATA_PATH("peakFile1.txt"))),
		 g2_(PeakReader(TEST_DATA_PATH("broken_peakFile.txt"))),
		 g3_(PeakReader(TEST_DATA_PATH("peakFile2.txt")))
	{
	}

	void virtual SetUp(){
	
	}
	public:
	PeakReader g_;
	PeakReader g2_;
	PeakReader g3_;

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

TEST_F(PeakReaderTest, WrongChromosome){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr5",30,300);	
	ASSERT_THROW(g_.findOverlappingPeaks(coordinates),std::invalid_argument);
	ASSERT_EQ(0,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, NotPresentOnChromosome){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",3000,3100);	
	ASSERT_THROW(g_.findOverlappingPeaks(coordinates),std::invalid_argument);
	ASSERT_EQ(0,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, PeakFileNotSorted){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",30,420);	
	ASSERT_THROW(g3_.findOverlappingPeaks(coordinates),std::invalid_argument);
}

TEST_F(PeakReaderTest, InvalidGenomicCoordinate){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr5",300,30);	
	ASSERT_THROW(g_.findOverlappingPeaks(coordinates),std::invalid_argument);
	ASSERT_EQ(0,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, PresentOnChromosome){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",25,420);	
	g_.findOverlappingPeaks(coordinates);
	ASSERT_EQ(3,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, PresentOnChromosome2C){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",26,29);	
	g_.findOverlappingPeaks(coordinates);
	ASSERT_EQ(1,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, PresentOnChromosome3C){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",270,290);	
	g_.findOverlappingPeaks(coordinates);
	ASSERT_EQ(1,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, PresentOnChromosome4C){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",180,240);	
	g_.findOverlappingPeaks(coordinates);
	ASSERT_EQ(1,g_.getOverlappingPeaks().size());
}

TEST_F(PeakReaderTest, PresentOnChromosomeR){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",25,420);	
	g_.findOverlappingPeaks(coordinates);
	std::vector<std::pair<unsigned int, unsigned int> > peaks;
	peaks = g_.getOverlappingPeaks();
	ASSERT_EQ(24,std::get<0>(peaks[0]));
	ASSERT_EQ(30,std::get<1>(peaks[0]));
	ASSERT_EQ(220,std::get<0>(peaks[1]));
	ASSERT_EQ(280,std::get<1>(peaks[1]));
	ASSERT_EQ(400,std::get<0>(peaks[2]));
	ASSERT_EQ(499,std::get<1>(peaks[2]));
}

TEST_F(PeakReaderTest, PresentOnChromosome2R){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",26,29);	
	g_.findOverlappingPeaks(coordinates);
	std::vector<std::pair<unsigned int, unsigned int> > peaks;
	peaks = g_.getOverlappingPeaks();
	ASSERT_EQ(24,std::get<0>(peaks[0]));
	ASSERT_EQ(30,std::get<1>(peaks[0]));
}

TEST_F(PeakReaderTest, PresentOnChromosome3R){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",270,290);	
	g_.findOverlappingPeaks(coordinates);
	std::vector<std::pair<unsigned int, unsigned int> > peaks;
	peaks = g_.getOverlappingPeaks();
	ASSERT_EQ(220,std::get<0>(peaks[0]));
	ASSERT_EQ(280,std::get<1>(peaks[0]));
}

TEST_F(PeakReaderTest, PresentOnChromosome4R){
	std::tuple<std::string,unsigned int, unsigned int> coordinates = std::make_tuple("chr3",180,240);	
	g_.findOverlappingPeaks(coordinates);
	std::vector<std::pair<unsigned int, unsigned int> > peaks;
	peaks = g_.getOverlappingPeaks();
	ASSERT_EQ(220,std::get<0>(peaks[0]));
	ASSERT_EQ(280,std::get<1>(peaks[0]));
}
