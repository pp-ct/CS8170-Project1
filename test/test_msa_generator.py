import sys
import os
sys.path.append("..\src")

from script import main

main()

for file in os.listdir("."):
	if file.endswith(".msa"):
		print("MSA generated for protein")