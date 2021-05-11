#include <omp.h>
#include "gtest/gtest.h"
#include "../SPAN/Data.cpp"
#include "../SPAN/Fraction.cpp"
#include "../SPAN/Span.cpp"
#include "../SPAN/Segment.cpp"
#include "../core/SPANInputGenerator.cpp"                                                                                                                                                                                                                                      
#include "../core/ExpressionReader.cpp"

#include "config.h"
class SPANTest : public ::testing::Test{
	protected:
	SPANTest(){
	}

	void virtual SetUp(){
	}

	public:
	Data d_;
	SPAN s_;
};

TEST_F(SPANTest, getElement){
     d_.read_file(TEST_DATA_PATH("Integrated_featureMatrix_ENSG00000107581_5000_V2.tab"),true,'g',1,false);
	SPAN s_ = SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > resultVector = s_.runSpan(d_,1,2);
	ASSERT_EQ(resultVector.size(),7);
}

TEST_F(SPANTest,setDat){
	ExpressionReader exp_ = ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample.txt"));
	exp_.loadExpressionData("ENSG00000184990");
	std::tuple<std::string, unsigned int, unsigned int,std::string> coordinates_=std::make_tuple("chr3",1,40,"+");
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	d_.setData(sig.getInputMatrix(),true,'g',1,false);
	SPAN s_ = SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > resultVector = s_.runSpan(d_);
	EXPECT_TRUE(resultVector.size()>7 && resultVector.size()<10);
}

TEST_F(SPANTest,setMC){
	ExpressionReader exp_ = ExpressionReader(TEST_DATA_PATH("Expression_Data_MultiClass.txt"));
	exp_.loadExpressionData("ENSG00000184990");
	std::tuple<std::string, unsigned int, unsigned int,std::string> coordinates_=std::make_tuple("chr3",1,40,"+");
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	d_.setData(sig.getInputMatrix(),true,'g',1,false);
	SPAN s_ = SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > resultVector = s_.runSpan(d_);
	ASSERT_EQ(resultVector.size(),9);
}

TEST_F(SPANTest,coordConversion){
	ExpressionReader exp_ = ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample.txt"));
	exp_.loadExpressionData("ENSG00000184990");
	std::tuple<std::string, unsigned int, unsigned int,std::string> coordinates_=std::make_tuple("chr3",4,41,"+");
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	d_.setData(sig.getInputMatrix(),true,'g',1,false);
	SPAN s_ = SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > resultVector = s_.runSpan(d_);
	std::vector<std::pair<unsigned int, unsigned int> > genomeConv= s_.convertSegmentationToGenomicCoordinates(resultVector,coordinates_);
	ASSERT_EQ(genomeConv[0].first,4);
	ASSERT_EQ(genomeConv[genomeConv.size()-1].second,41);
}
