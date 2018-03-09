import os
import shutil
import tarfile


idx_out = "/Users/pacera/test_pipeline/test_tars/output.tar.gz"


idx_out_pregz = idx_out.replace('.tar.gz', '.tar')
        #idx_out_pregz = path/file.tar


index_dir = idx_out.replace('.tar.gz', '')
#index_dir = path/file

os.mkdir(index_dir)

idx_split = index_dir.split("/")
#idx_split = ["path" , "file"]
"""
shutil.move(amb_loc, index_dir)
shutil.move(ann_loc, index_dir)
shutil.move(bwt_loc, index_dir)
shutil.move(pac_loc, index_dir)
shutil.move(sa_loc, index_dir)
"""
index_folder = idx_split[-1]

tar = tarfile.open(idx_out_pregz, "w")
tar.add(index_dir, arcname=index_folder) #arcname??
tar.close()
