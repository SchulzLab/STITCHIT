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
	stitchitCommand="./build/core/STITCH -b "+bwPath+" -a "+annoFile+ " -d "+disExp+" -o "+orgiExp+" -s "+sizeFile+" -f "+outPut+" -g "+geneID+ " -w 25000 -c 8 -p 0.05 -z 10 -m Spearman"

	sampleList=getSampleList(bwPath)
	nOfSamples=len(sampleList)
	nOfTraining=int(0.8*nOfSamples)
	nOfTest=int(0.2*nOfSamples)

	for fold in range(1,nFolds+1):
		print("Monte Carlo Cross Validation fold: ",fold)
		trainingIndex = random.sample(range(0,nOfSamples),nOfTraining)
		trainingVector = generateTrainingVec(trainingIndex,nOfSamples)
		testVector = generateTestFromTrainingVec(trainingVector)
		#Exceuting STITCHIT
		print(stitchitCommand+" -y "+toString(trainingVector))
		#Performing model fitting
		#ToDo
		regionFile="test.bed"
		#ToDo
		#Generating Test data
		testDataGenCommand="./build/core/TEST_DATA_GENERATOR -k "+regionFile+" -b "+bwPath+" -g "+geneID+ " -o "+orgiExp+" -f "+outPut+" -y "+toString(testVector)
		print(testDataGenCommand)
		#Assing performance on hold out data
		#ToDo

main()
