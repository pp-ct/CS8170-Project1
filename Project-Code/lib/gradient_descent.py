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
            # TODO: why are these abs?
            d_normal_val = np.abs(d_normal(d_ij, mu_ij, sd_ij))
            d_distance_val = np.abs(d_distance(r1, r2))
            d_normal_d_r += - np.log(d_normal_val * d_distance_val + 0.0001)

    return d_normal_d_r


def update_r(a, distance_matrix, residue_matrix, folding_depth):
    new_residue_matrix = np.copy(residue_matrix)

    # only fold up to given depth
    for i in range(folding_depth):
        new_residue_matrix = np.copy(residue_matrix)

        r_update = d_normal_d_r(distance_matrix[:folding_depth,:folding_depth], residue_matrix[:folding_depth], i)
        gradient_projected_position = - a * r_update + residue_matrix[i]
        relative_residue_position = residue_matrix[i - 1] if i > 0 else residue_matrix[i + 1]
        relative_residue_vector = gradient_projected_position - relative_residue_position
        relative_residue_vector = 3.8 * relative_residue_vector / np.linalg.norm(relative_residue_vector)

        new_residue_matrix[i] = relative_residue_vector + relative_residue_position

        residue_matrix = new_residue_matrix

    # "Pull" residues that have not begun to be optimized yet
    last_residue = new_residue_matrix[folding_depth - 1]
    next_residue_vector = last_residue - new_residue_matrix[folding_depth]
    update_later_residue_vector = next_residue_vector - 3.8 * next_residue_vector / np.linalg.norm(next_residue_vector)
    new_residue_matrix[folding_depth:] += update_later_residue_vector

    return residue_matrix






def initialize_residue_matrix(length):
    residue_matrix = np.ones((length, 3))
    for i in range(length):
        residue_matrix[i][0] = i * 3.8
        residue_matrix[i][1] = i * np.random.normal(loc=0.0001, scale=0.001)
        residue_matrix[i][2] = i * np.random.normal(loc=0.0001, scale=0.001)

    return residue_matrix