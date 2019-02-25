



with open("/home/pacera/MuG/CHi-C/tests/data/test_run_chicago/test.rmap", "r") as rmap_file:
	with open("/home/pacera/MuG/CHi-C/tests/data/test_run_chicago/test.baitmap", "r") as baitmap_file:
		with open("/home/pacera/MuG/CHi-C/tests/data/test_run_chicago/test2.baitmap", "w") as out_baitmap:
			baitmaps = []
			for line in baitmap_file:
				line_hdl = line.rstrip().split("\t")
				baitmaps.append(line_hdl[1]+" "+line_hdl[2])

			for line in rmap_file:
				line_hdl = line.rstrip().split("\t")
				if line_hdl[1]+" "+line_hdl[2] in baitmaps:
					out_baitmap.write("{}\t{}\t{}\t{}\n".format("21",
															  line_hdl[1],
															  line_hdl[2],
															  line_hdl[3]))
