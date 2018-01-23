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
          d_.read_file(TEST_DATA_PATH("Integrated_featureMatrix_ENSG00000107581_5000.tab"),true,false,'g',1,false);
	SPAN s_ = SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > resultVector = s_.runSpan(d_,1,2);
	ASSERT_EQ(resultVector.size(),7);
}

TEST_F(SPANTest,setDat){
	ExpressionReader exp_ = ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample.txt"));
	exp_.loadExpressionData("ENSG00000184990");
	std::tuple<std::string, unsigned int, unsigned int> coordinates_=std::make_tuple("chr3",1,49);
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	d_.setData(sig.getInputMatrix(),true,false,'g',1,false);
	SPAN s_ = SPAN();
	std::vector<std::pair<unsigned int, unsigned int> > resultVector = s_.runSpan(d_);
	ASSERT_EQ(resultVector.size(),1);
}

