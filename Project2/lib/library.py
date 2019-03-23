import os
from Bio import PDB
import numpy as np
import urllib
import matplotlib.pyplot as plt


aa_dict = {
    'A':'ALA',
    'C':'CYS',
    'D':'ASP',
    'E':'GLU',
    'F':'PHE',
    'G':'GLY',
    'H':'HIS',
    'I':'ILE',
    'K':'LYS',
    'L':'LEU',
    'M':'MET',
    'N':'ASN',
    'P':'PRO',
    'Q':'GLN',
    'R':'ARG',
    'S':'SER',
    'T':'THR',
    'V':'VAL',
    'W':'TRP',
    'Y':'TYR',
    'ALA':'A',
    'CYS':'C',
    'ASP':'D',
    'GLU':'E',
    'PHE':'F',
    'GLY':'G',
    'HIS':'H',
    'ILE':'I',
    'LYS':'K',
    'LEU':'L',
    'MET':'M',
    'ASN':'N',
    'PRO':'P',
    'GLN':'Q',
    'ARG':'R',
    'SER':'S',
    'THR':'T',
    'VAL':'V',
    'TRP':'W',
    'TYR':'Y'
}


def write_to_file(file_name, data, mode):
    with open(file_name, mode) as file:
        file.write(data + "\n")


def write_stream(filename, data_stream):
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


def download_pdb(save_folder, hit_id, chain_id):
    file_name = "{}_{}.pdb".format(hit_id, chain_id)
    if not os.path.isfile(os.path.join(save_folder, file_name)):
        try:
            urllib.request.urlretrieve("http://files.rcsb.org/view/" + hit_id.lower() + ".pdb",
                                       os.path.join(save_folder, file_name))
        except:
            print("PDB {} not found.".format(hit_id))
            return False

    return True


def get_fasta_for_id(save_folder, hit_id, chain_id):
    file_name = "{}_{}.fasta".format(hit_id, chain_id)
    if not os.path.isfile(os.path.join(save_folder, file_name)):
        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId={}&chainId={}"\
            .format(hit_id, chain_id)

        dest = os.path.join(save_folder, file_name)
        urllib.request.urlretrieve(url, dest)


def output_distance_matrix(save_folder, distance_matrix, prefix=None):
    map = plt.imshow(distance_matrix, cmap='gray')
    plt.colorbar(map)
    plt.savefig(os.path.join(save_folder, "{}distance_matrix.png".format((prefix + "_") if prefix else "")))
    plt.close()


def write_pdb(path, residue_matrix, seq, name):
    "ATOM      2  CA  MET A   1      -2.112  80.521  14.723  1.00 27.66           C  "
    with open(os.path.join(path, "{}.pdb".format(name)), 'w+') as pdb_file:
        for i, residue in enumerate(residue_matrix):
            # pdb_file.write(
            #     "ATOM\t{}\tCA\t{}\tA\t{}\t{}\t{}\t{}\t1\t1\tC\n".format(
            #         i,
            #         aa_dict[seq[i]],
            #         i,
            #         np.round(residue_matrix[i][0], 4),
            #         np.round(residue_matrix[i][1], 4),
            #         np.round(residue_matrix[i][2], 4)
            #     )
            # )
            pdb_file.write("{}{} {} {} {}{}    {}{}{}{}{}\n".format(
                "ATOM".ljust(6),
                str(i).rjust(5),
                ("CA" if seq[i] == 'G' else "CB").center(4),
                aa_dict[seq[i]].ljust(3),
                "A".rjust(1),
                str(i).rjust(4),
                str('%8.3f' % (float(residue_matrix[i][0]))).rjust(8),
                str('%8.3f' % (float(residue_matrix[i][1]))).rjust(8),
                str('%8.3f' % (float(residue_matrix[i][2]))).rjust(8),
                str('%6.2f' % (float(1))).rjust(6),
                str('%6.2f' % (float(1))).rjust(6),
                "C".rjust(12)
            ))

        pdb_file.write("TER")

    # parser = PDB.PDBParser()
    # chain = parser.get_structure(id='temp', file="{}.pdb".format(name))
    # io = PDB.PDBIO()
    # io.set_structure(chain)
    # io.save("{}.pdb".format(name))