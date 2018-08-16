import sys

def main():
	mappingDict={}
	mappingfile=open(sys.argv[2],"r")
	for l in mappingfile:
		s=l.split()
		mappingDict[s[1]]=s[0]
	mappingfile.close()

	infile=open(sys.argv[1],"r")
	header=infile.readline()
	for l in infile:
		s=l.split(";")
		temp=s[0]+"\t"+s[3]+"\t"+s[4]
		for i in xrange(0,1+(len(s)-10)/2):
			geneName=s[2*i+9].replace("connected_gene=","")
			if (geneName in mappingDict):
				print(temp+"\t"+mappingDict[geneName]+"\t"+geneName)
	infile.close()

main()
