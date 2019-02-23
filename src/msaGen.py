import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

inputDir = "../data/target/"
blastDir = "../data/blast/"
outputDir = "../data/msa/"

def blastTarget():
	"""
	This reads a fasta file containing a protein sequence ../data/target/ and queries NCBI to retrieve similar proteins. The output is written to ../data/blast/.
	"""
	for file in os.listdir(inputDir):
		if file.endswith("fasta"):
			sequence = SeqIO.read(inputDir+file, format = "fasta")
			#Query the sequence on the blast website using the pdb dictionary
			resultHandle = NCBIWWW.qblast("blastp", "pdbaa", sequence.seq)
			#Save the BLAST output to the output directory
			with open(blastDir + file.split(".")[0] + "msa.xml", "a") as outHandle:
				outHandle.write(resultHandle.read())
			resultHandle.close()

def parseMsa(eValueThreshold = 3.04, print_match=False, print_subject=False):
	"""
	This parses the blast .xml files from *../data/blast/* and generates a fasta file in ../data/msa/. 
	eValueThreshold: This sets the maximum eValue allowed for sequences returned.
	print_match: This writes match data to the fasta file.
	print_subject: This writes the input sequence to the fast file.
	"""
	for file in os.listdir(blastDir):
		if file.endswith("xml"):
			blast = open(blastDir + file)
			#parse the XML data to retrieve data from each xml file
			#Open a file handle to use to write all alignments into
			fastaHandle = open(outputDir + "msa.fasta", "a")
			blastRecords = NCBIXML.parse(blast)
			for record in blastRecords:
				#Get the alignment object in the record
				for alignment in record.alignments:
					for hsp in alignment.hsps:
						if hsp.expect < eValueThreshold:
							#Write the title
							#print(hsp)
							fastaHandle.write(">"+alignment.accession[0:4]+":"+alignment.accession[5]+"|PDBID|CHAIN|SEQUENCE"+ "\n")
							#Write the length
							#fastaHandle.write("length: " + str(alignment.length) + "\n")
							#Write the e-value
							#fastaHandle.write("e value: " + str(hsp.expect)+ "\n")
							#Write the input input, match and the output
							fastaHandle.write(hsp.query+ "\n")
							if(print_match==True):
								fastaHandle.write(hsp.match+ "\n")
							if(print_subject==True):
								fastaHandle.write(hsp.sbjct+ "\n")
							#fastaHandle.write("\n")
			fastaHandle.close()
			
blastTarget()
parseMsa()