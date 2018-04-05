

with open("./baits.sam", "r") as file:
	for line in file:
		line = line.rstrip().split("\t")
		if line[0][0] != "@":
			print(line)
			break
