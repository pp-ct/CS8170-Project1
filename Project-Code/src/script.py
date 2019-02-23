import os
from os.path import dirname, abspath
import sys
import textwrap
DIRECTORY = dirname(dirname(abspath(__file__)))
sys.path.append(DIRECTORY)
from lib import library 
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


TARGET_DIR = "../data/target/" #place target fasta files here i.e. one fasta file per target protein
OUTPUT_DIR = "../data/msa/"


def parse_result(msa_xml_file):
	result=open(msa_xml_file,"r")
	output+= OUTPUT_DIR+"msa.fasta";
	records= NCBIXML.parse(result)
	item=next(records)
	for alignment in item.alignments:
		for hsp in alignment.hsps:
			library.write_to_file(output, ">"+alignment.accession[0:4]+":"+alignment.accession[5]+"|PDBID|CHAIN|SEQUENCE","a")
			library.write_to_file(output, textwrap.fill(hsp.query, 60), "a")
		library.write_to_file(output+"msa.fasta", "", "a")

def main():
	for filename in os.listdir(TARGET_DIR):
	    if filename.endswith(".fasta") : 
	    	 target_file = TARGET_DIR+filename
	    	 msa_file = OUTPUT_DIR+"msa.xml";
	    	 record = SeqIO.read(target_file, format="fasta")
	    	 result_handle = NCBIWWW.qblast("blastp", "pdbaa", record.seq)
	    	 library.create_dir(OUTPUT_DIR)
	    	 library.write_stream(msa_file, result_handle)
	    	 parse_result(msa_file)
	    else:
	        continue


main()
