import numpy as np
from Bio import PDB
import os
import pickle

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
    'TYR':'Y',

    'X':'MSE',
    'MSE':'X'
}


class Alignment:

    def __init__(self, hit_id, chain_id, query_range, hit_range, query_seq, hit_seq, midline):
        self.hit_id = hit_id
        self.chain_id = chain_id
        self.query_range = query_range
        self.hit_range = hit_range

        self.query_seq = query_seq
        self.hit_seq = hit_seq
        self.midline = midline

        self.query_span = query_range[1] - query_range[0]
        self.hit_span = hit_range[1] - hit_range[0]
        self.query_hit_dict, self.hit_query_dict = self.generate_query_hit_dict()
        self.pdb_path = ""

        self.ca_hit_distance_matrix = None
        self.no_hit_distance_matrix = None

    def generate_query_hit_dict(self):
        query_hit_dict = {}
        hit_query_dict = {}
        q_index, h_index = self.query_range[0] - 1, self.hit_range[0] - 1
        for q, h, m in zip(self.query_seq, self.hit_seq, self.midline):
            # if gap in query sequence
            if q == '-':
                h_index += 1
            # if gap in hit sequence
            elif h == '-':
                q_index += 1
            # if residues align
            elif m != ' ':
                query_hit_dict[q_index] = h_index
                hit_query_dict[h_index] = q_index
                q_index += 1
                h_index += 1
            # if residues do not align
            else:
                q_index += 1
                h_index += 1

        return query_hit_dict, hit_query_dict

    def generate_hit_distance_matrix(self, type='CA'):
        hit_distance_matrix = np.zeros((self.hit_span, self.hit_span))

        r1_type = 'CA'
        r2_type = 'CA'
        if type == 'NO':
            r1_type = 'N'
            r2_type = 'O'

        parser = PDB.PDBParser()
        chains = parser.get_structure(id='temp', file=self.pdb_path)[0]
        chain = chains[self.chain_id] if self.chain_id in chains else chains['A']

        for residue1 in chain.get_residues():
            r1 = residue1.id[1]
            if self.hit_range[0] < r1 < self.hit_range[1]:
                for residue2 in chain.get_residues():
                    r2 = residue2.id[1]
                    if self.hit_range[0] < r2 < self.hit_range[1] and r1_type in residue1 and r2_type in residue2:
                        distance = residue1[r1_type] - residue2[r2_type]
                        hit_distance_matrix[r1 - self.hit_range[0]][r2 - self.hit_range[0]] = distance

        return hit_distance_matrix


    def build_hit_matrices(self, distance_matrix_dir):
        if self.ca_hit_distance_matrix is None:
            pickle_file = os.path.join(distance_matrix_dir, "{}_{}-CA.pickle".format(self.hit_id, self.chain_id))
            if os.path.exists(pickle_file):
                with open(pickle_file, 'rb') as distance_matrix_file:
                    self.ca_hit_distance_matrix = pickle.load(distance_matrix_file)
            else:
                self.ca_hit_distance_matrix = self.generate_hit_distance_matrix(type='CA')
                os.system("touch {}".format(pickle_file))
                with open(pickle_file, 'wb') as distance_matrix_file:
                    pickle.dump(self.ca_hit_distance_matrix, distance_matrix_file)
        if self.no_hit_distance_matrix is None:
            pickle_file = os.path.join(distance_matrix_dir, "{}_{}-NO.pickle".format(self.hit_id, self.chain_id))
            if os.path.exists(pickle_file):
                with open(pickle_file, 'rb') as distance_matrix_file:
                    self.no_hit_distance_matrix = pickle.load(distance_matrix_file)
            else:
                self.no_hit_distance_matrix = self.generate_hit_distance_matrix(type='NO')
                os.system("touch {}".format(pickle_file))
                with open(pickle_file, 'wb') as distance_matrix_file:
                    pickle.dump(self.no_hit_distance_matrix, distance_matrix_file)


    def get_distance_for_query_residues(self,  r1, r2, type='CA'):
        if self.ca_hit_distance_matrix is None:
            self.ca_hit_distance_matrix = self.generate_hit_distance_matrix(type='CA')
        if self.no_hit_distance_matrix is None:
            self.no_hit_distance_matrix = self.generate_hit_distance_matrix(type='NO')

        hit_distance_matrix = self.ca_hit_distance_matrix
        if type == 'NO':
            hit_distance_matrix = self.no_hit_distance_matrix

        if r1 in self.query_hit_dict and r2 in self.query_hit_dict:
            indexed_r1 = self.query_hit_dict[r1] - self.hit_range[0]
            indexed_r2 = self.query_hit_dict[r2] - self.hit_range[0]
            return hit_distance_matrix[indexed_r1][indexed_r2]
        else:
            return False