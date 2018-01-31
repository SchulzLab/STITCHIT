#include <omp.h>
#include "gtest/gtest.h"
#include "../core/BinSelection.cpp"
#include "../core/ExpressionReader.cpp"
#include "../core/CorComp.cpp"

#include "config.h"

class CORCOMPTest : public ::testing::Test{
	protected:
	CORCOMPTest(){
	}

	void virtual SetUp(){
		v1 = {1,2,3,4,5};
		v2 = {10,11,12,13,14};
		v3 = {1,2,3,4,5,6};
		v4 = {60,50,40,30,20,10};		
		v5 = {1,1,1,1,1,1};
		v6 = {1,2,4,3,5};
		v7={};
		v8={};
		v9={0.01,0.011,0.009,0.012,0.3,0.005};
		v10={0.012,0.0092,0.00134,0.0118,0.298,0.004};
	}
	public:
	std::vector<double> v1,v2,v3,v4,v5,v6,v7,v8,v9,v10;
};


TEST_F(CORCOMPTest, constructor){
	CorComp cc (v1,v2);
	ASSERT_EQ(cc.getObservationCount(),5);
}

TEST_F(CORCOMPTest, constructor2){
	ASSERT_THROW(CorComp cc(v1,v3), std::runtime_error);
}

TEST_F(CORCOMPTest, MeanTest1){
	CorComp cc (v1,v2);
	ASSERT_EQ(cc.computeMean(v1),3);
	ASSERT_EQ(cc.computeMean(v2),12);	
}

TEST_F(CORCOMPTest, MeanTest2){
	CorComp cc (v3, v4);
	ASSERT_EQ(cc.computeMean(v3),3.5);
	ASSERT_EQ(cc.computeMean(v4),35);	
}

TEST_F(CORCOMPTest, MeanTest3){
	CorComp cc (v3, v5);
	ASSERT_EQ(cc.computeMean(v3),3.5);
	ASSERT_EQ(cc.computeMean(v5),1);	
}

TEST_F(CORCOMPTest, VarTest1){
	CorComp cc (v1,v2);
	ASSERT_DOUBLE_EQ(cc.computeVariance(v1),2.5);
	ASSERT_DOUBLE_EQ(cc.computeVariance(v2),2.5);	
}

TEST_F(CORCOMPTest, VarTest2){
	CorComp cc (v3, v4);
	ASSERT_DOUBLE_EQ(cc.computeVariance(v3),3.5);
	ASSERT_DOUBLE_EQ(cc.computeVariance(v4),350);	
}

TEST_F(CORCOMPTest, VarTest3){
	CorComp cc (v3, v5);
	ASSERT_DOUBLE_EQ(cc.computeVariance(v3),3.5);
	ASSERT_DOUBLE_EQ(cc.computeVariance(v5),0);	
}

TEST_F(CORCOMPTest, PearsonTest1){
	CorComp cc (v1,v2);
	ASSERT_DOUBLE_EQ(cc.computePearsonCorrelation(),1.0);
}

TEST_F(CORCOMPTest, PearsonTest2){
	CorComp cc (v3,v4);
	ASSERT_DOUBLE_EQ(cc.computePearsonCorrelation(),-1);
}

TEST_F(CORCOMPTest, PearsonTest3){
	CorComp cc(v6,v1);
	ASSERT_DOUBLE_EQ(cc.computePearsonCorrelation(),0.9);
}

TEST_F(CORCOMPTest, PearsonTest4){
	CorComp cc(v3,v5);
	ASSERT_DOUBLE_EQ(cc.computePearsonCorrelation(),0.0);
}


TEST_F(CORCOMPTest, PearsonTest5){
	CorComp cc(v9,v10);
	ASSERT_NEAR(cc.computePearsonCorrelation(),0.99,0.01);
}

TEST_F(CORCOMPTest, RankTest1){
	CorComp cc(v1,v2);
	cc.convertToRanks();
	std::vector<double> r1 = cc.getRanksV1();
	std::vector<double> r2 = cc.getRanksV2();
	ASSERT_DOUBLE_EQ(r1[0],5);
	ASSERT_DOUBLE_EQ(r1[1],4);
	ASSERT_DOUBLE_EQ(r1[2],3);
	ASSERT_DOUBLE_EQ(r1[3],2);
	ASSERT_DOUBLE_EQ(r1[4],1);
	ASSERT_DOUBLE_EQ(r2[0],5);
	ASSERT_DOUBLE_EQ(r2[1],4);
	ASSERT_DOUBLE_EQ(r2[2],3);
	ASSERT_DOUBLE_EQ(r2[3],2);
	ASSERT_DOUBLE_EQ(r2[4],1);
}

TEST_F(CORCOMPTest, RankTest2){
	CorComp cc(v5,v4);
	cc.convertToRanks();
	std::vector<double> r1 = cc.getRanksV1();
	std::vector<double> r2 = cc.getRanksV2();
	ASSERT_DOUBLE_EQ(r1[0],3.5);
	ASSERT_DOUBLE_EQ(r1[1],3.5);
	ASSERT_DOUBLE_EQ(r1[2],3.5);
	ASSERT_DOUBLE_EQ(r1[3],3.5);
	ASSERT_DOUBLE_EQ(r1[4],3.5);
	ASSERT_DOUBLE_EQ(r1[5],3.5);
	ASSERT_DOUBLE_EQ(r2[0],1);
	ASSERT_DOUBLE_EQ(r2[1],2);
	ASSERT_DOUBLE_EQ(r2[2],3);
	ASSERT_DOUBLE_EQ(r2[3],4);
	ASSERT_DOUBLE_EQ(r2[4],5);
	ASSERT_DOUBLE_EQ(r2[5],6);
}

TEST_F(CORCOMPTest, RankTest3){
	CorComp cc(v1,v6);
	cc.convertToRanks();
	std::vector<double> r1 = cc.getRanksV1();
	std::vector<double> r2 = cc.getRanksV2();
	ASSERT_DOUBLE_EQ(r1[0],5);
	ASSERT_DOUBLE_EQ(r1[1],4);
	ASSERT_DOUBLE_EQ(r1[2],3);
	ASSERT_DOUBLE_EQ(r1[3],2);
	ASSERT_DOUBLE_EQ(r1[4],1);
	ASSERT_DOUBLE_EQ(r2[0],5);
	ASSERT_DOUBLE_EQ(r2[1],4);
	ASSERT_DOUBLE_EQ(r2[2],2);
	ASSERT_DOUBLE_EQ(r2[3],3);
	ASSERT_DOUBLE_EQ(r2[4],1);
}

TEST_F(CORCOMPTest, SpearmanTest1){
	CorComp cc (v1,v2);
	ASSERT_DOUBLE_EQ(cc.computeSpearmanCorrelation(),1.0);
}

TEST_F(CORCOMPTest, SpearmanTest2){
	CorComp cc (v3,v4);
	ASSERT_DOUBLE_EQ(cc.computeSpearmanCorrelation(),-1);
}

TEST_F(CORCOMPTest, SpearmanTest3){
	CorComp cc(v6,v1);
	ASSERT_DOUBLE_EQ(cc.computeSpearmanCorrelation(),0.9);
}

TEST_F(CORCOMPTest, SpearmanTest4){
	CorComp cc(v3,v5);
	ASSERT_DOUBLE_EQ(cc.computeSpearmanCorrelation(),0.0);
}

TEST_F(CORCOMPTest, SpearmanTest5){
	CorComp cc(v9,v10);
	ASSERT_NEAR(cc.computeSpearmanCorrelation(),0.7714,0.01);
}

TEST_F(CORCOMPTest, FisherZTest1){
	CorComp cc(v9,v10);
	double cor = cc.computePearsonCorrelation();
	std::pair<double,double> fisherZ = cc.getFisherZ(cor);
	ASSERT_NEAR(fisherZ.first,4.298838,0.01);
	ASSERT_NEAR(fisherZ.second,0.5773503,0.01);
}

TEST_F(CORCOMPTest, FisherZTest2){
	CorComp cc(v9,v10);
	double cor= cc.computeSpearmanCorrelation();
	std::pair<double,double> fisherZ = cc.getFisherZ(cor);
	ASSERT_NEAR(fisherZ.first,1.023776,0.01);
	ASSERT_NEAR(fisherZ.second,0.5773503,0.01);
}

TEST_F(CORCOMPTest, FisherZTest3){
	CorComp cc(v6,v1);
	double cor = cc.computePearsonCorrelation();
	std::pair<double,double> fisherZ = cc.getFisherZ(cor);
	ASSERT_NEAR(fisherZ.first,1.47222,0.01);
	ASSERT_NEAR(fisherZ.second,0.707107,0.01);
}

TEST_F(CORCOMPTest, cdf1){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(0),0.5,0.01);
}

TEST_F(CORCOMPTest, cdf2){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(1),0.8413447,0.01);
}

TEST_F(CORCOMPTest, cdf3){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(1.1),0.8643,0.01);
}

TEST_F(CORCOMPTest, cdf4){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(1.5),0.93319,0.01);
}

TEST_F(CORCOMPTest, cdf5){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(2),0.97724,0.01);
}

TEST_F(CORCOMPTest, cdf6){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(2.5),0.9937903,0.01);
}

TEST_F(CORCOMPTest, cdf7){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(3.0),0.9986,0.01);
}

TEST_F(CORCOMPTest, cdf8){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(4.0),0.9999683,0.01);
}

TEST_F(CORCOMPTest, cdf9){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(4.38),0.9999941,0.01);
}

TEST_F(CORCOMPTest, cdf10){
	CorComp cc (v1,v2);
	ASSERT_NEAR(cc.cdf(-0.5),0.3085,0.01);
}

TEST_F(CORCOMPTest, pValueTest1){
	CorComp cc (v6,v1);
	double cor = cc.computePearsonCorrelation();
	ASSERT_NEAR(cc.getPvalue(cor),0.01866973,0.01);
}

TEST_F(CORCOMPTest, pValueTest2){
	CorComp cc (v9,v10);
	double cor = cc.computePearsonCorrelation();
	ASSERT_NEAR(cc.getPvalue(cor),0.0,0.01);
}

TEST_F(CORCOMPTest, pValueTest3){
	CorComp cc (v9,v10);
	double cor = cc.computeSpearmanCorrelation();
	ASSERT_NEAR(cc.getPvalue(cor),0.03815232,0.01);
}
