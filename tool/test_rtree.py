
from rtree import index
"""
idx = index.Rtree("rtreetest")

with open("/Users/pacera/developing/C-HiC/tests/data/test_makeRmap/restriction_enzyme_test2.txt", "r") as out:
    counter = 0
    for line in out:
    	counter += 1
    	line = line.rstrip().split("\t")
    	#print(counter, (int(line[1]), int(line[0][3:]), int(line[2]), int(line[0][3:])))

    	idx.insert(counter, (int(line[1]), int(line[0][3:]), int(line[2]), int(line[0][3:])))

idx.close()

"""
idx = index.Rtree("/Users/pacera/developing/C-HiC/tests/data/test_makeRmap/rtree_file")

print(list(idx.intersection((50640197, 22, 50640456, 22))))