import os
import sys
import random

def getSampleList(path):
	return(os.listdir(path))

def generateTrainingVec(IndexVec,nOfSamples):
	sampleBinVec=[0]*nOfSamples
	for index in IndexVec:
		sampleBinVec[index]=1
	return(sampleBinVec)

def toString(vector):
	tmp=""
	for value in vector:
		tmp+=str(value)
	return(tmp)

def generateTestFromTrainingVec(vec):
	return(map(lambda x : abs(x-1),vec))

def main():
	bwPath=sys.argv[1]
	annoFile=sys.argv[2]
	disExp=sys.argv[3]
	orgiExp=sys.argv[4]
	outPut=sys.argv[5]
	sizeFile=sys.argv[6]
	geneID=sys.argv[7]
	nFolds=int(sys.argv[8])
	ROutput=sys.argv[9]
	stitchitCommand="./build/core/STITCH -b "+bwPath+" -a "+annoFile+ " -d "+disExp+" -o "+orgiExp+" -s "+sizeFile+" -f "+outPut+" -g "+geneID+ " -w 25000 -c 8 -p 0.05 -z 10 -m Spearman"
	modelFitCommand="Rscript ModelFitting_NestedSTITCHIT.R --outDir="+ROutput+" --dataDir="+outPut+" --response=Expression"

	sampleList=getSampleList(bwPath)
	nOfSamples=len(sampleList)
	nOfTraining=int(0.8*nOfSamples)
	nOfTest=int(0.2*nOfSamples)

	#Cleaning up output folders
	os.system("rm "+outPut+"/*")
	os.system("rm "+ROutput+"/*")

	for fold in range(1,nFolds+1):
		print("Monte Carlo Cross Validation fold: ",fold)
		trainingIndex = random.sample(range(0,nOfSamples),nOfTraining)
		trainingVector = generateTrainingVec(trainingIndex,nOfSamples)
		testVector = generateTestFromTrainingVec(trainingVector)
		print(stitchitCommand+" -y "+toString(trainingVector))

		#Exceuting STITCHIT
		os.system(stitchitCommand+" -y "+toString(trainingVector))

		#Performing model fitting
		print(modelFitCommand)
		os.system(modelFitCommand)
		regionFile=ROutput+"/Selected_Regions_ElasticNet_Segmentation_"+geneID+"_Spearman_10.bed"
		regionFile2=ROutput+"/Selected_Regions_ElasticNet_Segmentation_"+geneID+"_Spearman_10Tab.bed"

		reformat="cut -f 2 "+regionFile+ " | sed 's/\./\t/g' > "+regionFile2
		print(reformat)
		os.system(reformat)

		#Generating Test data
		testDataGenCommand="./build/core/TEST_DATA_GENERATOR -k "+regionFile2+" -b "+bwPath+" -g "+geneID+ " -o "+orgiExp+" -f "+outPut+" -y "+toString(testVector)
		print(testDataGenCommand)
		os.system(testDataGenCommand)

		#Assing performance on hold out data
		testExecuteCommand="Rscript ModelTest_NestedSTITCHIT_ElaNet.R --outDir="+ROutput+" --testdata="+outPut+" --response=Expression --targetGeneID="+geneID+" >> "+sys.argv[10]+"_Nested_STICHIT_ElaNet.txt"

		print(testExecuteCommand)
		os.system(testExecuteCommand)
main()
