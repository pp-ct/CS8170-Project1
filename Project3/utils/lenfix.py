def fix_length(fileName):
	#Output file is the 
	outFile = open(fileName[:-4]+"_fix" + fileName[-4:], "w")
	inFile = open(fileName, "r")
	while(inFile.readline()):
		outFile.write(inFile.readline()[:64] + "\n")
		
fix_length("complex.1.pdb")