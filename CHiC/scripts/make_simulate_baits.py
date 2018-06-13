"""
Takes the genome in a fasta file and select 22.000 sequences,
that are going to represents the capture sequence for the baits.
"""

import random


select_lines = set(random.sample(range(0, 1625477), 10000))


with open("/home/pablo/MuG/C-HiC/tests/data/test_makeBaitmap/chr21_hg19.fa", "r") as file_in:
    counter_line = 0
    flag = "close"
    for line in file_in:
        if line[0] in ("N", ">"):
            continue
        counter_line += 1
        line = line.rstrip()
        if flag == "open":
            for char in line:
                if len(bait) < 120:
                    bait.append(char.upper())
                else:
                    flag = "close"
                    with open("/home/pablo/MuG/C-HiC/tests/data/test_makeBaitmap/baits.fa", "a") as file_out:
                        file_out.write(">probe"+str(counter_line)+"\n")
                        file_out.write("".join(bait)+"\n")
                    break

        else:
            if counter_line in select_lines:
                bait = []
                for char in line:
                    if char == "N":
                        break
                    bait.append(char.upper())
                flag = "open"


                #if len(bait) == 120:
                #   print(bait)
                    #with open("/Users/pacera/developing/baits.txt", "a") as file_out:
                    #   print("".join(bait), file = file_out)