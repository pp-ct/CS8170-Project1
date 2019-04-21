import os
import glob
import numpy as np
import subprocess
import re
import tqdm


def evaluate_irmsd(target_name, method, pred_num):
	if not os.path.exists("structures/{}/{}".format(target_name, method)):
		return 0

	target_name_trimmed = re.split('-', target_name)[0]
	profit_script = sorted(list(glob.glob("structures/{}/*.izone".format(target_name))))[0]
	bound_pdb = "structures/{}/{}_bound.pdb".format(target_name, target_name_trimmed)
	pred_pdb = "structures/{}/{}/{}.pdb".format(target_name, method, pred_num)	

	args = ["./ProFitV3.1/profit", "-f", profit_script, bound_pdb, pred_pdb]
	profit_output = str.split(subprocess.check_output(args).decode('utf-8'))

	score = str(profit_output[profit_output.index("RMS:") + 1])
	return score
	

profit_path = "ProFitV3.1/profit"
target_files = sorted(list(glob.glob("structures/*")))

target_names = []
zdock_results = np.array(np.zeros((len(target_files), 10)), dtype=str)
cluspro_results = np.array(np.zeros((len(target_files), 10)), dtype=str)

for i, file in tqdm.tqdm(list(enumerate(target_files))):
	target_names.append(os.path.split(file)[1])

	for j in range(1, 11):
		zdock_results[i][j - 1] =  evaluate_irmsd(target_names[-1], "zdock", j)
		cluspro_results[i][j - 1] =  evaluate_irmsd(target_names[-1], "cluspro", j)


target_names = np.array(target_names).reshape((len(target_names), 1))
zdock_results = np.concatenate((target_names, zdock_results), axis=1)
cluspro_results = np.concatenate((target_names, cluspro_results), axis=1)

np.savetxt("zdock_results.csv", zdock_results, fmt="%s", delimiter=',')
np.savetxt("cluspro_results.csv", cluspro_results, fmt="%s", delimiter=',')
