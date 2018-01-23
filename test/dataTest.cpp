#include "gtest/gtest.h"
#include "../SPAN/Data.cpp"
#include "../core/SPANInputGenerator.cpp"                                                                                                                                                                                                                                      
#include "../core/ExpressionReader.cpp"
#include "config.h"

class DATATest : public ::testing::Test{
	protected:
	DATATest(){
	}

	void virtual SetUp(){
	
	}

	public:
	Data d_;
};

TEST_F(DATATest, ReadInput){
          d_.read_file(TEST_DATA_PATH("Integrated_featureMatrix_ENSG00000107581_5000.tab"),true,false,'g',1,false);
	SUCCEED();
}

TEST_F(DATATest, checkRowCount){
          d_.read_file(TEST_DATA_PATH("Integrated_featureMatrix_ENSG00000107581_5000.tab"),true,false,'g',1,false);
	ASSERT_EQ(d_.getRowCount(),60);
}

TEST_F(DATATest, checkColCount){
          d_.read_file(TEST_DATA_PATH("Integrated_featureMatrix_ENSG00000107581_5000.tab"),true,false,'g',1,false);
	ASSERT_EQ(d_.getColCount(),5001);
}

TEST_F(DATATest, getElement){
          d_.read_file(TEST_DATA_PATH("Integrated_featureMatrix_ENSG00000107581_5000.tab"),true,false,'g',1,false);
	ASSERT_EQ(d_.getElement(0,0),1);
	ASSERT_EQ(d_.getElement(1,0),0);
	ASSERT_EQ(d_.getElement(11,0),56.2);
	ASSERT_EQ(d_.getElement(17,0),2.7);
	ASSERT_EQ(d_.getElement(31,5000),33.6);
	ASSERT_EQ(d_.getElement(56,5000),1.8);
}

TEST_F(DATATest,setDat){
	ExpressionReader exp_ = ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample.txt"));
	exp_.loadExpressionData("ENSG00000184990");
	std::tuple<std::string, unsigned int, unsigned int> coordinates_=std::make_tuple("chr3",10,40);
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
          std::vector<std::string> sampleNames = sig.getSampleNames();
          std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
          std::vector<double> dataVector={0,0,0,0,0,0,0,0,0,0,10,10,10,10,10,20,20,20,20,20,10,10,10,10,10,0,0,0,0,0};
	d_.setData(sig.getInputMatrix(),true,false,'g',1,false);
	ASSERT_EQ(d_.getRowCount(),4);
          ASSERT_EQ(d_.getColCount(),30);
          ASSERT_EQ(d_.getRow(std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11")-sampleNames.begin()),dataVector);
}

