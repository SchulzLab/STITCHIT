import os
import sys
import subprocess

def main():
	integratedFilePath=sys.argv[1]
	integratedFiles=os.listdir(integratedFilePath)
	for intFile in integratedFiles:
		print("Computing bins in file: "+intFile)
		geneID=str(intFile.split(".")[0])
		command="./span.run -i "+integratedFilePath+intFile+" -c -p 10 -o ../"+geneID+".span.txt"
		print(command)
		subprocess.call(command,shell=True)
main()
