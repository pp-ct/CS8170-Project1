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
            d_normal_val = d_normal(d_ij, mu_ij, sd_ij)
            d_distance_val = d_distance(r1, r2)

            d_normal_d_r -= d_normal_val * d_distance_val + 0.0001

    return d_normal_d_r


def detect_clashes(residue_matrix):
    clashes = []
    for i in range(len(residue_matrix)):
        for j in range(len(residue_matrix)):
            if i > j and distance(residue_matrix[i], residue_matrix[j]) < 3.4:
                clashes.append([i, j, distance(residue_matrix[i], residue_matrix[j])])

    return sorted(clashes, key=lambda clash: clash[2], )


def resolve_clashes(residue_matrix, folding_depth):
    new_residue_matrix = np.copy(residue_matrix)

    clashes = detect_clashes(residue_matrix)
    while clashes:
        clash = clashes[0]
        # print("resolving: ",clash)

        chain_center = np.mean(new_residue_matrix[:folding_depth], axis=0)
        r1_index = clash[0] if distance(clash[0], chain_center) > distance(clash[1], chain_center) else clash[1]
        r2_index = clash[1] if distance(clash[0], chain_center) > distance(clash[1], chain_center) else clash[0]

        r1 = new_residue_matrix[r1_index]
        r2 = new_residue_matrix[r2_index]
        clash_distance = clash[2]

        # for i in range(1, len(new_residue_matrix)):
        #     print(distance(new_residue_matrix[i], new_residue_matrix[i - 1]))
        # print("\n\n")

        resolution_unit_vector = (r1 - r2) / np.linalg.norm(r1 - r2)
        resolution_vector = 3.4 * resolution_unit_vector
        new_residue_matrix[r1_index] = r2 + resolution_vector

        for i in range(r1_index - 1, -1, -1):
            adjusted_position = mediate_residue_placement(new_residue_matrix[i + 1],
                                                          new_residue_matrix[i])
            new_residue_matrix[i] = adjusted_position
        for i in range(r1_index + 1, len(new_residue_matrix), 1):
            adjusted_position = mediate_residue_placement(new_residue_matrix[i - 1],
                                                          new_residue_matrix[i])
            new_residue_matrix[i] = adjusted_position

        # for i in range(1, len(new_residue_matrix)):
        #     print(distance(new_residue_matrix[i], new_residue_matrix[i - 1]))
        clashes = detect_clashes(new_residue_matrix)
    return new_residue_matrix


def mediate_residue_placement(relative_position, projected_position):
    relative_residue_vector = projected_position - relative_position
    relative_residue_vector = (3.8 * np.random.normal(loc=1.0, scale=0.01)) \
                              * relative_residue_vector / np.linalg.norm(relative_residue_vector)
    new_residue_position = relative_position + relative_residue_vector
    return new_residue_position


def update_r(a, b, distance_matrices, residue_matrix, r_update_previous, folding_depth):
    new_residue_matrix = np.copy(residue_matrix)
    r_update_current = np.zeros(residue_matrix.shape)

    # only fold up to given depth
    for i in range(folding_depth):
        new_residue_matrix = np.copy(residue_matrix)

        r_update = a * np.sum([d_normal_d_r(d_mat[:folding_depth,:folding_depth], residue_matrix[:folding_depth], i)
            for d_mat in distance_matrices], axis=0) + b * r_update_previous[i]
        # r_update = d_normal_d_r(distance_matrix[:folding_depth,:folding_depth], residue_matrix[:folding_depth], i) \
        #            + b * r_update_previous[i]
        r_update_current[i] = r_update

        gradient_projected_position = -1 * r_update + residue_matrix[i]
        relative_residue_position = residue_matrix[i - 1] if i > 0 else residue_matrix[i + 1]
        # new_residue_matrix[i] = mediate_residue_placement(residue_matrix[i], relative_residue_position,
        #                                                   gradient_projected_position)
        new_residue_matrix[i] = mediate_residue_placement(relative_residue_position,
                                                          gradient_projected_position)

        residue_matrix = new_residue_matrix

    # "Pull" residues that have not begun to be optimized yet
    last_residue = new_residue_matrix[folding_depth - 1]
    next_residue_vector = last_residue - new_residue_matrix[folding_depth]
    update_later_residue_vector = next_residue_vector - 3.8 * next_residue_vector / np.linalg.norm(next_residue_vector)
    new_residue_matrix[folding_depth:] += update_later_residue_vector

    return residue_matrix, r_update_current


def ranked_update_r(a, distance_matrices, residue_matrix, folding_depth):
    new_residue_matrix = np.copy(residue_matrix)

    for _ in range(folding_depth):
        r_update_list = []
        for i in range(folding_depth):
            r_update = a * np.sum([d_normal_d_r(d_mat[:folding_depth,:folding_depth], residue_matrix[:folding_depth], i)
                               for d_mat in distance_matrices], axis=0)
            r_update_list.append([i, r_update, np.linalg.norm(r_update)])
        r_update_list = sorted(r_update_list, key=lambda r_u: r_u[2], reverse=True)
        update_residue_index = r_update_list[0][0]
        new_residue_matrix[update_residue_index] += -1 * r_update_list[0][1]

        for i in range(update_residue_index - 1, -1, -1):
            adjusted_position = mediate_residue_placement(new_residue_matrix[i + 1],
                                                          new_residue_matrix[i])
            new_residue_matrix[i] = adjusted_position
        for i in range(update_residue_index + 1, folding_depth, 1):
            adjusted_position = mediate_residue_placement(new_residue_matrix[i - 1],
                                                          new_residue_matrix[i])
            new_residue_matrix[i] = adjusted_position

    # "Pull" residues that have not begun to be optimized yet
    last_residue = new_residue_matrix[folding_depth - 1]
    next_residue_vector = last_residue - new_residue_matrix[folding_depth]
    update_later_residue_vector = next_residue_vector - 3.8 * next_residue_vector / np.linalg.norm(
        next_residue_vector)
    new_residue_matrix[folding_depth:] += update_later_residue_vector

    return new_residue_matrix


def initialize_residue_matrix(length):
    residue_matrix = np.ones((length, 3))
    for i in range(length):
        residue_matrix[i][0] = i * 3.8
        residue_matrix[i][1] = i * np.random.normal(loc=0, scale=0.001)
        residue_matrix[i][2] = i * np.random.normal(loc=0, scale=0.001)

    return residue_matrix