#include <omp.h>
#include "gtest/gtest.h"
#include "../core/BinSelection.cpp"
#include "../core/ExpressionReader.cpp"
#include "../core/CorComp.cpp"

#include "config.h"

class BINSELECTIONTest : public ::testing::Test{
	protected:
	BINSELECTIONTest(){
	}

	void virtual SetUp(){
	}

	public:
};


//NOTE that coordinates are 1based, they are transformed to the zero base within the BinSelectionProcess, the intervals are both closed, i.e. [A,B].
TEST_F(BINSELECTIONTest, constructor){
	BinSelection bs = BinSelection(TEST_DATA_PATH("BigWigFiles"));
	ASSERT_EQ(bs.getMeanSignal().size(),0);
}


TEST_F(BINSELECTIONTest,matrixSize){
	BinSelection bs = BinSelection(TEST_DATA_PATH("BigWigFiles"));
	std::vector<std::pair<unsigned int, unsigned int> > segments;
	segments.push_back(std::make_pair(1,5));
	segments.push_back(std::make_pair(6,10));
	segments.push_back(std::make_pair(11,20));
	std::string chrom;
	chrom="chr3";
	bs.computeMeanSignal(chrom,segments);
	std::vector<std::vector<double> > signal;
	signal = bs.getMeanSignal();
	ASSERT_EQ(signal.size(),4);
	ASSERT_EQ(signal[0].size(),3);
	ASSERT_EQ(signal[1].size(),3);
	ASSERT_EQ(signal[2].size(),3);
}


TEST_F(BINSELECTIONTest,EntriesNotMatchingStepSize){
	BinSelection bs = BinSelection(TEST_DATA_PATH("BigWigFiles"));
	std::vector<std::pair<unsigned int, unsigned int> > segments;
	segments.push_back(std::make_pair(4,7));
	segments.push_back(std::make_pair(9,12));
	segments.push_back(std::make_pair(14,17));
	std::string chrom;
	chrom="chr3";
	bs.computeMeanSignal(chrom,segments);
	std::vector<std::vector<double> > signal;
	signal = bs.getMeanSignal();
	std::vector<std::string> sampleNames;
	sampleNames =  bs.getSampleNames();
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11")-sampleNames.begin()][0],0);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11")-sampleNames.begin()][1],0);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11")-sampleNames.begin()][2],0);

	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00DFM11")-sampleNames.begin()][0],10);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00DFM11")-sampleNames.begin()][1],10);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00DFM11")-sampleNames.begin()][2],10);

	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00Y4Y11")-sampleNames.begin()][0],12.5);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00Y4Y11")-sampleNames.begin()][1],17.5);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S00Y4Y11")-sampleNames.begin()][2],25);

	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S004BT12")-sampleNames.begin()][0],12.5);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S004BT12")-sampleNames.begin()][1],17.5);
	ASSERT_EQ(signal[std::find(sampleNames.begin(),sampleNames.end(),"B_S004BT12")-sampleNames.begin()][2],25);
}

TEST_F(BINSELECTIONTest,getSignalVectorBySegment){
	BinSelection bs = BinSelection(TEST_DATA_PATH("BigWigFiles"));
	std::vector<std::pair<unsigned int, unsigned int> > segments;
	segments.push_back(std::make_pair(4,7));
	segments.push_back(std::make_pair(9,12));
	segments.push_back(std::make_pair(14,17));
	std::string chrom;
	chrom="chr3";
	bs.computeMeanSignal(chrom,segments);
	std::vector<std::vector<double> > signal;
	signal = bs.getMeanSignal();
	std::vector<std::string> sN = bs.getSampleNames();

	std::vector<double> zS = bs.getSignalVectorBySegment(0);
	ASSERT_EQ(zS[std::find(sN.begin(),sN.end(),"B_S00CWT11")-sN.begin()],0);
	ASSERT_EQ(zS[std::find(sN.begin(),sN.end(),"B_S00DFM11")-sN.begin()],10);
	ASSERT_EQ(zS[std::find(sN.begin(),sN.end(),"B_S00Y4Y11")-sN.begin()],12.5);
	ASSERT_EQ(zS[std::find(sN.begin(),sN.end(),"B_S004BT12")-sN.begin()],12.5);

	std::vector<double> fS = bs.getSignalVectorBySegment(1);
	ASSERT_EQ(fS[std::find(sN.begin(),sN.end(),"B_S00CWT11")-sN.begin()],0);
	ASSERT_EQ(fS[std::find(sN.begin(),sN.end(),"B_S00DFM11")-sN.begin()],10);
	ASSERT_EQ(fS[std::find(sN.begin(),sN.end(),"B_S00Y4Y11")-sN.begin()],17.5);
	ASSERT_EQ(fS[std::find(sN.begin(),sN.end(),"B_S004BT12")-sN.begin()],17.5);

	std::vector<double> sS = bs.getSignalVectorBySegment(2);
	ASSERT_EQ(sS[std::find(sN.begin(),sN.end(),"B_S00CWT11")-sN.begin()],0);
	ASSERT_EQ(sS[std::find(sN.begin(),sN.end(),"B_S00DFM11")-sN.begin()],10);
	ASSERT_EQ(sS[std::find(sN.begin(),sN.end(),"B_S00Y4Y11")-sN.begin()],25);
	ASSERT_EQ(sS[std::find(sN.begin(),sN.end(),"B_S004BT12")-sN.begin()],25);

}

TEST_F(BINSELECTIONTest,getSignalVectorByExpression){
	BinSelection bs = BinSelection(TEST_DATA_PATH("BigWigFiles"));
	std::vector<std::pair<unsigned int, unsigned int> > segments;
	segments.push_back(std::make_pair(4,7));
	segments.push_back(std::make_pair(9,12));
	segments.push_back(std::make_pair(14,17));
	std::string chrom;
	chrom="chr3";
	bs.computeMeanSignal(chrom,segments);
	std::vector<std::vector<double> > signal;
	signal = bs.getMeanSignal();
	std::vector<std::string> sN = bs.getSampleNames();
	ExpressionReader exp (TEST_DATA_PATH("Expression_Data_Sample.txt"));
	exp.loadExpressionData("ENSG00000184990");
	std::vector<double> expression = bs.getExpressionVectorByNames(exp.getExpressionMap());
	ASSERT_EQ(expression[0],exp.getExpressionMap()[sN[0]]);
	ASSERT_EQ(expression[1],exp.getExpressionMap()[sN[1]]);
	ASSERT_EQ(expression[2],exp.getExpressionMap()[sN[2]]);	
	ASSERT_EQ(expression[3],exp.getExpressionMap()[sN[3]]);
}
