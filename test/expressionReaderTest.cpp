#include "gtest/gtest.h"
#include "../core/ExpressionReader.cpp"
#include "config.h"

class ExpressionReaderTest : public ::testing::Test{
	protected:
	ExpressionReaderTest()
		:expR_(ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample.txt"))),
		 expR2_(ExpressionReader(TEST_DATA_PATH("Expression_Data_Phony_Sample.txt"))),
		 expR3_(ExpressionReader(TEST_DATA_PATH("Expression_Data_Does_Not_Exist.txt"))),
		 expR4_(ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample_Double.txt")))
	{
	}

	void virtual SetUp(){
	
	}
	public:
	ExpressionReader expR_;
	ExpressionReader expR2_;
	ExpressionReader expR3_;
	ExpressionReader expR4_;

};

TEST_F(ExpressionReaderTest, Constructor){
	SUCCEED();
}

TEST_F(ExpressionReaderTest, getfileName){
	ASSERT_EQ(TEST_DATA_PATH("Expression_Data_Sample.txt"),expR_.getFilename());
}

TEST_F(ExpressionReaderTest, loadGenomeSizeFile){
	expR_.loadExpressionData("ENSG00000184990");
	SUCCEED();
}

TEST_F(ExpressionReaderTest, ExpressionFileDoesNotExist){
	ASSERT_THROW(expR3_.loadExpressionData("ENSG00000184990"),std::invalid_argument);
}

TEST_F(ExpressionReaderTest, ExpressionFileIsPhony){
	ASSERT_THROW(expR2_.loadExpressionData("ENSG00000184990"),std::invalid_argument);
}

TEST_F(ExpressionReaderTest, IsEmpty){
	ASSERT_TRUE(expR_.getExpressionMap().empty());
}

TEST_F(ExpressionReaderTest, IsNotEmpty){
	expR_.loadExpressionData("ENSG00000184990");
	ASSERT_FALSE(expR_.getExpressionMap().empty());
}

TEST_F(ExpressionReaderTest, SIVA1){
	expR_.loadExpressionData("ENSG00000184990");
	std::map<std::string, double> expressionData = expR_.getExpressionMap();
	ASSERT_EQ(expressionData["B_S004BT12"],0);
	ASSERT_EQ(expressionData["B_S00CWT11"],1);
	ASSERT_EQ(expressionData["B_S00DFM11"],1);
	ASSERT_EQ(expressionData["B_S00Y4Y11"],1);
}

TEST_F(ExpressionReaderTest, SIVA2){
	expR4_.loadExpressionData("ENSG00000184990");
	std::map<std::string, double> expressionData = expR4_.getExpressionMap();
	ASSERT_EQ(expressionData["B_S004BT12"],0.3);
	ASSERT_EQ(expressionData["B_S00CWT11"],110.);
	ASSERT_EQ(expressionData["B_S00DFM11"],103.1);
	ASSERT_EQ(expressionData["B_S00Y4Y11"],612.001);
}

TEST_F(ExpressionReaderTest, GeneDoesNotExist){
	ASSERT_THROW(expR_.loadExpressionData("ENSG11111111111"),std::invalid_argument);
}
