import os
import shutil
import tarfile
import shlex
import subprocess



idx_out = "/Users/pacera/test_pipeline/mg-process-test1/tool/output.tar.gz"

"""

idx_out_pregz = idx_out.replace('.tar.gz', '.tar')
        #idx_out_pregz = path/file.tar


index_dir = idx_out.replace('.tar.gz', '')
#index_dir = path/file

try:
	os.mkdir(index_dir)
except:
	pass

idx_split = index_dir.split("/")
#idx_split = ["path" , "file"]

index_folder = idx_split[-1]

tar = tarfile.open(idx_out_pregz, "w")
tar.add(index_dir, arcname=index_folder) #arcname??
tar.close()


command_line = 'pigz ' + idx_out_pregz
args = shlex.split(command_line)
print(args)
process = subprocess.Popen(args)
process.wait()
"""

print("BS - idx_out", idx_out, idx_out.replace('.tar.gz', ''))
idx_out_pregz = idx_out.replace('.tar.gz', '.tar')

index_dir = idx_out.replace('.tar.gz', '')

try:
	os.mkdir(index_dir)
except:
	pass


idx_split = index_dir.split("/")


index_folder = idx_split[-1]

tar = tarfile.open(idx_out_pregz, "w")
tar.add(index_dir, arcname=index_folder)
tar.close()

command_line = 'pigz ' + idx_out_pregz
args = shlex.split(command_line)
print(args)
process = subprocess.Popen(args)
process.wait()

shutil.rmtree(index_dir)

