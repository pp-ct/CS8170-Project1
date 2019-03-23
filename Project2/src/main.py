import sys
import os
from os.path import dirname, abspath
from Bio import SeqIO
from tqdm import tqdm
import numpy as np

DIRECTORY = dirname(dirname(abspath(__file__)))
sys.path.append(DIRECTORY)

from lib import library
from lib import rrmaps
from lib import gradient_descent
from glob import glob

DATA_DIR = "../data/"
TARGET_DIR = "../data/target/"
OUTPUT_DIR = "../data/target/{}/output"
for file in glob(TARGET_DIR + "*"):
    target_name = os.path.split(file)[1]
    output_dir = OUTPUT_DIR.format(target_name)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    target_file = os.path.join(TARGET_DIR, target_name, "target.fasta")
    record = SeqIO.read(target_file, format="fasta")

    bound_file = os.path.join(TARGET_DIR, target_name, "bound.txt")
    distance_matrix = rrmaps.distance_matrix(bound_file, target_name)

    # if not os.path.exists(OUTPUT_DIR):
    #     os.mkdir(OUTPUT_DIR)
    # # output_dir = os.path.join(OUTPUT_DIR, record.name)
    # if not os.path.exists(output_dir):
    #     os.mkdir(output_dir)

    residue_matrix = gradient_descent.initialize_residue_matrix(len(record.seq))
    new_residue_matrix = residue_matrix.copy()
    r_update_previous = np.zeros(new_residue_matrix.shape)

    iterations = 10000
    output_interval = 100
    a = 0.01
    b = 0.0
    length = len(record.seq)
    for i in tqdm(range(0, iterations + 1)):
        # update residue matrix for residues up to given depth
        folding_depth = min(length - 1, int((2 * i / iterations) * length))
        # folding_depth = length - 1
        new_residue_matrix, r_update_previous = gradient_descent.update_r(a, b, [distance_matrix], new_residue_matrix,
                                                                          r_update_previous,
                                                                          folding_depth)
        new_residue_matrix = gradient_descent.resolve_clashes(new_residue_matrix, folding_depth)

        if i % output_interval == 0:
            library.write_pdb(output_dir, new_residue_matrix, record.seq, "structure-{}".format(i))

    # library.write_pdb(output_dir, new_residue_matrix, record.seq, "final_structure")
    os.system("python ../lib/zhang_python_scripts/reindex_pdb.py {} {} {} -clean=True".format(
        target_file,
        os.path.join(output_dir, "structure-{}.pdb".format(iterations)),
        os.path.join(output_dir, "final_structure.pdb")
    ))








