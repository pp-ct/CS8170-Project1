import os
import sys
import textwrap
sys.path.append('../')
from lib import library 
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


TARGET_DIR = "../data/target/"
OUTPUT_DIR = "../data/msa/"

def write_output(filename,handle, data_stream):
	with open(filename, "w+") as handle:
	    handle.write(data_stream)
	    data_stream.close()

def parse_result(msa_xml_file):
	result=open(msa_xml_file,"r")
	records= NCBIXML.parse(result)
	item=next(records)
	for alignment in item.alignments:
		for hsp in alignment.hsps:
			library.write_to_file(OUTPUT_DIR+"msa.fasta", ">"+alignment.accession[0:4]+":"+alignment.accession[5]+"|PDBID|CHAIN|SEQUENCE","a")
			library.write_to_file(OUTPUT_DIR+"msa.fasta", textwrap.fill(hsp.query, 60), "a")
		library.write_to_file(OUTPUT_DIR+"msa.fasta", "", "a")

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




# # result=open("output.xml","r")
# # records= NCBIXML.parse(result)
# # item=next(records)
# # for alignment in item.alignments:
# # 	for hsp in alignment.hsps:
# # 		print(hsp.match)
#         # print(hsp.sbjct)
# with open("output.xml", "w+") as save_to:
#     save_to.write(result_handle.read())
#     result_handle.close()


# result=open("output.xml","r")
# records= NCBIXML.parse(result)
# item=next(records)
# for alignment in item.alignments:
# 	print(alignment.accession[0:4])
# 	# for hsp in alignment.hsps:
# 		# print(hsp)
# 		# print(hsp.match)
  #       print(hsp.sbjct)

print("done")
# https://www.biostars.org/p/98693/
# env_nr
# nr
# pataa
# pdbaa
# refseq_protein
# swissprot