

"""
Trim the name of the fastq files
"""

with open("SRR3535023_1.fastq", "r") as file_in:
	with open("SRR3535023_1_trim.fastq", "w") as file_out:
		counter = 0
		for line in file_in:
			counter +=1
			line_hdl = line.rstrip()
			if line[:7] == "@SRR3535":
				line_hdl = line_hdl.split("\t")
				file_out.write("{}\n".format(line_hdl[0]))
			else:
				file_out.write("{}\n".format(line_hdl))

			if counter == 16:
				break

