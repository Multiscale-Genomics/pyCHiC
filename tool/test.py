"""
with open("ref_test.txt", "r") as file:
    genome_dict = {}
    sequence = ""
    counter = 0
    for line in file:
        print(line)
        counter += 1
        line = line.rstrip()
        if line[0] is ">":
            if int(line[4:]) not in genome_dict:
                genome_dict[int(line[4:])] = []
                chromo = int(line[4:])
                if sequence is "":
                    continue
                    else:
                     genome_dict[chromo] = sequence
                            sequence = ""
                            continue
                sequence += line
            #the last chromosome was not added
            genome_dict[chromo] = sequence
return genome_dict
"""

with open("ref_test.txt", "r") as file:
    genome_dict = {}
    sequence = ""
    for line in file:
        line = line.rstrip()
        if line[0] is ">":
            if sequence is "":
                genome_dict[int(line[4:])] = []
                chromo = int(line[4:])
                continue
            else:
                genome_dict[chromo] = sequence
                chromo = int(line[4:])
                sequence = ""
                continue

        sequence += line

    #Ad last chromosome
    genome_dict[chromo] = sequence



print(genome_dict)