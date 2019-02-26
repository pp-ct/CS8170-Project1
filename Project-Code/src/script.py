import os
from os.path import dirname, abspath
import sys
import glob
import re
import textwrap
import numpy as np

DIRECTORY = dirname(dirname(abspath(__file__)))
sys.path.append(DIRECTORY)
from lib import library
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Alignment import Alignment

DATA_DIR = "../data"
TARGET_DIR = "../data/target/"  # place target fasta files here i.e. one fasta file per target protein
OUTPUT_DIR = "../data/msa/"
TEMPLATE_PDB_DIR = "../data/template_pdbs"
TEMPLATE_FASTA_DIR = "../data/template_fasta"


def parse_result(msa_xml_file, file_no):
    print('*' * 50)
    print("Parsing Alignment XML\n\n")

    result = open(msa_xml_file, "r")
    output = OUTPUT_DIR + "msa" + str(file_no) + ".fasta"
    records = NCBIXML.parse(result)
    item = next(records)
    for alignment in item.alignments:
        for hsp in alignment.hsps:
            library.write_to_file(output, ">" + alignment.accession[0:4] + ":" + alignment.accession[
                5] + "|PDBID|CHAIN|SEQUENCE", "a")
            library.write_to_file(output, hsp.query, "a")
            library.write_to_file(output, hsp.match, "a")
            library.write_to_file(output, hsp.sbjct, "a")
        library.write_to_file(output, "", "a")
    print(str(file_no) + " msa fasta file has been generated")


def gen_alignment_list(msa_xml_file):
    print('*' * 50)
    print("Generating Alignment List\n\n")

    result = open(msa_xml_file, "r")
    records = NCBIXML.parse(result)
    alignment_obj_list = []
    item = next(records)
    for alignment in item.alignments:
        for hsp in alignment.hsps:
            hit_id = alignment.accession[0:4]
            chain_id = alignment.accession[5]
            query_range = (hsp.query_start, hsp.query_end)
            hit_range = (hsp.sbjct_start, hsp.sbjct_end)
            query_seq = hsp.query
            hit_seq = hsp.sbjct
            midline = hsp.match

            # dont want native
            if hit_range == query_range:
                break

            alignment_obj = Alignment(hit_id, chain_id, query_range, hit_range,
                                      query_seq, hit_seq, midline)
            alignment_obj_list.append(alignment_obj)

    return alignment_obj_list


def download_pdbs_for_alignments(alignment_list):
    print('*' * 50)
    print("Downloading PDBs for Alignments\n\n")
    if not os.path.exists(TEMPLATE_PDB_DIR):
        os.mkdir(TEMPLATE_PDB_DIR)
    if not os.path.exists(TEMPLATE_FASTA_DIR):
        os.mkdir(TEMPLATE_FASTA_DIR)

    # some alignments will not actually have PDBs, so we need to discard them
    # (this can also happen if request times out)
    bad_alignments = []

    # download PDBs, and keep track of bad alignments
    print("Getting PDBs from database\n")
    for alignment in alignment_list:
        if not library.download_pdb(TEMPLATE_PDB_DIR, alignment.hit_id, alignment.chain_id):
            bad_alignments.append(alignment)

    for alignment in bad_alignments:
        alignment_list.remove(alignment)

    # need full fasta for renumbering PDBs
    print("Getting FASTAs from database\n")
    for alignment in alignment_list:
        library.get_fasta_for_id(TEMPLATE_FASTA_DIR, alignment.hit_id, alignment.chain_id)

    # now reindex PDBs
    print("Reindexing PDBs\n")
    for alignment in alignment_list:
        hit_chain_id = "{}_{}".format(alignment.hit_id, alignment.chain_id)
        if not os.path.exists(os.path.join(TEMPLATE_PDB_DIR, hit_chain_id + ".reindex.pdb")):
            os.system("python ../lib/zhang_python_scripts/reindex_pdb.py {} {} {} -clean=True".format(
                os.path.join(TEMPLATE_FASTA_DIR, hit_chain_id + ".fasta"),
                os.path.join(TEMPLATE_PDB_DIR, hit_chain_id + ".pdb"),
                os.path.join(TEMPLATE_PDB_DIR, hit_chain_id + ".reindex.pdb")
            ))
        alignment.pdb_path = os.path.join(TEMPLATE_PDB_DIR, hit_chain_id + ".reindex.pdb")


def build_target_distance_pdfs(length, alignment_list):
    print('*' * 50)
    print("Building Target Distance PDF Matrix\n\n")

    # this LxLx2 array will contain mean and SD for each pair of residues
    target_distance_pdfs = np.ndarray((length, length, 2))

    # this is the LxL matrix of lists (sort of)
    # each list contains all of the distances for the aligned pair
    for i in range(length):
        for j in range(length):
            ij_distance_list = []
            for alignment in alignment_list:
                alignment_ij_distance = alignment.get_distance_for_query_residues(i, j)
                if alignment_ij_distance:
                    ij_distance_list.append(alignment_ij_distance)

            if ij_distance_list:
                target_distance_pdfs[i][j][0] = np.mean(ij_distance_list)
                target_distance_pdfs[i][j][1] = np.std(ij_distance_list)
            else:
                target_distance_pdfs[i][j][0] = np.nan
                target_distance_pdfs[i][j][1] = np.nan
            # print if you want to see it
            # print(i, j, target_distance_pdfs[i][j])

    return target_distance_pdfs


def main():
    library.create_dir(TARGET_DIR)
    for filename in os.listdir(TARGET_DIR):
        file_no = 0
        if filename.endswith(".fasta"):
            file_no += 1
            target_file = TARGET_DIR + filename
            msa_file = OUTPUT_DIR + "msa.xml"
            record = SeqIO.read(target_file, format="fasta")
            # result_handle = NCBIWWW.qblast("blastp", "pdbaa", record.seq)
            # library.create_dir(OUTPUT_DIR)
            # library.write_stream(msa_file, result_handle)

            alignment_list = gen_alignment_list(msa_file)
            download_pdbs_for_alignments(alignment_list)
            target_distance_pdfs = build_target_distance_pdfs(len(record.seq), alignment_list)

            print(target_distance_pdfs[:,:,0].shape)
            library.output_distance_matrix(DATA_DIR, target_distance_pdfs[:,:,0])
        else:
            continue


main()
