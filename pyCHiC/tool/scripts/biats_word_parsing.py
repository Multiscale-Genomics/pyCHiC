
"""


"""
import re

with open("baits_cancer_loci.txt", "r") as file_in:
    starting = []
    ending = []
    chromo = []
    counter = 0
    set_counter = 0
    for line in file_in:
        counter += 1
        line_hdl = line.rstrip()
        if line_hdl.startswith("rs"):
            set_counter = counter
            chomosome = re.search("\d+", prev_line)
            chromo.append(chomosome.group(0))

        if set_counter+1 == counter:
            starting.append(line_hdl)

        elif set_counter+2 == counter:
            ending.append(line_hdl)
        prev_line = line_hdl

import pandas as pd

starting = [int("".join(i.split(","))) for i in starting[1:]]
ending = [int("".join(i.split(","))) for i in ending[1:]]


df = pd.DataFrame({1:chromo, 2:starting, 3:ending})

df.to_csv("dataframe_baits.baitmap", sep="\t", header=False, index=None)
