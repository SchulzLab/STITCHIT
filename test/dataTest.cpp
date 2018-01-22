#include "gtest/gtest.h"
#include "../SPAN/Data.cpp"
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
