import os
from Bio import PDB
import numpy as np
import urllib

def write_to_file(file_name, data, mode):
	with open(file_name, mode) as file:  
		file.write(data+"\n") 

def write_stream(filename,data_stream):
	with open(filename, "w+") as handle:
	    handle.write(data_stream.read())
	    data_stream.close()

def create_dir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def load_ca_distnaces_from_pdb(path, length, chain_id=None):
    parser = PDB.PDBParser()
    chain = parser.get_structure(id='temp', file=path)[0]
    if chain_id is not None:
        chain = parser.get_structure(id='temp', file=path)[0][chain_id]
    distance_matrix = np.zeros((length, length))

    for r1, residue1 in enumerate(chain):
        for r2, residue2 in enumerate(chain):
            if residue1 != residue2:
                # compute distance between CA atoms
                try:
                    distance = residue1['CA'] - residue2['CA']
                    distance_matrix[r1][r2] = distance
                except KeyError:
                    continue

    return distance_matrix

def download_pdb(id, save_folder):
    if not os.path.isfile(os.path.join(save_folder, id + '.pdb')):
        try:
            urllib.request.urlretrieve("http://files.rcsb.org/view/" + id.lower() + ".pdb", os.path.join(save_folder, id + '.pdb'))
        except:
            print('PDB {} not found.'.format(id))