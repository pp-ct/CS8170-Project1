import numpy as np



def distance(r1, r2):
    return np.sqrt(np.sum(np.power(r1 - r2, 2)))


def d_distance(r1, r2):
    inv_distance = 1 / distance(r1,r2)
    return (r1 - r2) * inv_distance


def normal(d_ij, mu_ij, sd_ij):
    return np.exp(-0.5 * np.power((d_ij - mu_ij) / sd_ij, 2)) / (sd_ij * np.sqrt(2 * np.pi))


def d_normal(d_ij, mu_ij, sd_ij):
    return np.exp(- 0.5 * np.power((d_ij - mu_ij) / sd_ij, 2)) * \
           (mu_ij - d_ij) / \
           (np.power(sd_ij, 3) * np.sqrt(2 * np.pi))


def d_normal_d_r(distance_matrix, residue_matrix, r_index):

    d_normal_d_r = 0
    for i in range(len(distance_matrix)):
        if i != r_index \
                and not np.isnan(distance_matrix[r_index][i][0]) \
                and not np.isnan(distance_matrix[r_index][i][1]):
            r1 = residue_matrix[r_index]
            r2 = residue_matrix[i]
            d_ij = distance(r1, r2)
            mu_ij = distance_matrix[r_index][i][0]
            sd_ij = distance_matrix[r_index][i][1]
            d_normal_val = d_normal(d_ij, mu_ij, sd_ij)
            d_distance_val = d_distance(r1, r2)
            d_normal_d_r += - np.log(d_normal_val * d_distance_val)

    return d_normal_d_r


def initialize_residue_matrix(length):
    residue_matrix = np.zeros((length, 3))
    for i in range(length):
        residue_matrix[i][0] = i * 3.800001

    return residue_matrix