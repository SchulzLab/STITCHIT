#include "gtest/gtest.h"
#include "../core/SPANInputGenerator.cpp"
#include "../core/ExpressionReader.cpp"
#include "config.h"

class SPANInputGeneratorTest : public ::testing::Test{
	protected:
	SPANInputGeneratorTest()
		 :exp_(ExpressionReader(TEST_DATA_PATH("Expression_Data_Sample.txt")))
	{
		exp_.loadExpressionData("ENSG00000184990");
		coordinates_=std::make_tuple("chr3",10,40,"+");
		coordinatesOOR_=std::make_tuple("chr3",1000,1030,"-");
	}

	void virtual SetUp(){

	}
	public:
		ExpressionReader exp_;
		std::tuple<std::string, unsigned int, unsigned int,std::string> coordinates_;
		std::tuple<std::string, unsigned int, unsigned int,std::string> coordinatesOOR_;

};

TEST_F(SPANInputGeneratorTest, Constructor){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	SUCCEED();
}

TEST_F(SPANInputGeneratorTest, VariablesAreEmpty){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	ASSERT_EQ(sig.getSampleNames().size(),0);
	ASSERT_EQ(sig.getInputMatrix().size(),0);
}


TEST_F(SPANInputGeneratorTest, AllFilesRead){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	ASSERT_TRUE(std::find(sampleNames.begin(),sampleNames.end(),"B_S004BT12") != sampleNames.end());
	ASSERT_TRUE(std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11") != sampleNames.end());
	ASSERT_TRUE(std::find(sampleNames.begin(),sampleNames.end(),"B_S00DFM11") != sampleNames.end());
	ASSERT_TRUE(std::find(sampleNames.begin(),sampleNames.end(),"B_S00Y4Y11") != sampleNames.end());
	ASSERT_FALSE(std::find(sampleNames.begin(),sampleNames.end(),"B_S1111112") != sampleNames.end());
}

TEST_F(SPANInputGeneratorTest, NoPhonyFilesIncluded){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	ASSERT_FALSE(std::find(sampleNames.begin(),sampleNames.end(),"B_S1111112") != sampleNames.end());
	ASSERT_EQ(sampleNames.size(),4);
}

TEST_F(SPANInputGeneratorTest, NoBigWigButWigFiles){
          SPANInputGenerator sig (TEST_DATA_PATH("WigFiles/"),exp_.getExpressionMap());
          ASSERT_THROW(sig.generateSPANInput(coordinates_),std::invalid_argument);
}

TEST_F(SPANInputGeneratorTest, PhonyPath){
          SPANInputGenerator sig (TEST_DATA_PATH("PhonyPath/"),exp_.getExpressionMap());
          ASSERT_THROW(sig.generateSPANInput(coordinates_),std::invalid_argument);
}

TEST_F(SPANInputGeneratorTest, ExpressionMatrixNotComplet){
	ExpressionReader temp_= ExpressionReader(TEST_DATA_PATH("Expression_Data_One_Sample_Missing.txt"));
          SPANInputGenerator sig (TEST_DATA_PATH("WigFiles/"),exp_.getExpressionMap());
          ASSERT_THROW(sig.generateSPANInput(coordinates_),std::invalid_argument);
}

TEST_F(SPANInputGeneratorTest, CorrectRowCountOfInputMatrix){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	ASSERT_EQ(inputMatrix.size(),4);
}

TEST_F(SPANInputGeneratorTest, CorrectColumnCountsOfEachRowInInputMatrix){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	ASSERT_EQ(inputMatrix[0].size(),31);
	ASSERT_EQ(inputMatrix[1].size(),31);
	ASSERT_EQ(inputMatrix[2].size(),31);
	ASSERT_EQ(inputMatrix[3].size(),31);
}


TEST_F(SPANInputGeneratorTest, CorrectOutputSample_B_S00CWT11){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	std::vector<double> dataVector={0,0,0,0,0,0,0,0,0,0,10,10,10,10,10,20,20,20,20,20,10,10,10,10,10,0,0,0,0,0,1};
	ASSERT_EQ(inputMatrix[std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11")-sampleNames.begin()],dataVector);
}

TEST_F(SPANInputGeneratorTest, CorrectOutputSample_B_S004BT12){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	std::vector<double> dataVector={20,20,20,20,20,30,30,30,30,30,20,20,20,20,20,0,0,0,0,0,10,10,10,10,10,5,5,5,5,5,0};
	ASSERT_EQ(inputMatrix[std::find(sampleNames.begin(),sampleNames.end(),"B_S004BT12")-sampleNames.begin()],dataVector);
}

TEST_F(SPANInputGeneratorTest, CorrectOutputSample_B_S00DFM11){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	std::vector<double> dataVector={10,10,10,10,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	ASSERT_EQ(inputMatrix[std::find(sampleNames.begin(),sampleNames.end(),"B_S00DFM11")-sampleNames.begin()],dataVector);
}

TEST_F(SPANInputGeneratorTest, CorrectOutputSample_B_S00Y4Y11){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinates_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	std::vector<double> dataVector={20,20,20,20,20,30,30,30,30,30,20,20,20,20,20,0,0,0,0,0,10,10,10,10,10,5,5,5,5,5,1};
	ASSERT_EQ(inputMatrix[std::find(sampleNames.begin(),sampleNames.end(),"B_S00Y4Y11")-sampleNames.begin()],dataVector);
}

TEST_F(SPANInputGeneratorTest, CorrectRowCountOfInputMatrix_OOR){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinatesOOR_);
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	ASSERT_EQ(inputMatrix.size(),4);
}

TEST_F(SPANInputGeneratorTest, CorrectColumnCountsOfEachRowInInputMatrix_OOR){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinatesOOR_);
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	ASSERT_EQ(inputMatrix[0].size(),31);
	ASSERT_EQ(inputMatrix[1].size(),31);
	ASSERT_EQ(inputMatrix[2].size(),31);
	ASSERT_EQ(inputMatrix[3].size(),31);
}

TEST_F(SPANInputGeneratorTest, CorrectOutputSample_B_S00CWT11_OOR){
	SPANInputGenerator sig (TEST_DATA_PATH("BigWigFiles/"),exp_.getExpressionMap());
	sig.generateSPANInput(coordinatesOOR_);
	std::vector<std::string> sampleNames = sig.getSampleNames();
	std::vector<std::vector<double>> inputMatrix =sig.getInputMatrix();
	std::vector<double> dataVector={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	ASSERT_EQ(inputMatrix[std::find(sampleNames.begin(),sampleNames.end(),"B_S00CWT11")-sampleNames.begin()],dataVector);
}
