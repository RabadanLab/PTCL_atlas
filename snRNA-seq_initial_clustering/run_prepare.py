import scanpy as sc
import numpy as np
import pandas as pd
import scanorama
import subprocess
import sys

CMD_PREFIX = "aws s3 cp s3://aw3291-projects/PTCL_Palomero/Transfer/CellRangerOuts_ForNotebook/cellbender_outs/"
#FOLDER = "sn01-Bank02-AITL"
FOLDER = sys.argv[1]
#H5_FILE = "sn01-Bank02-AITL_cellbender_output_filtered.h5"
H5_FILE = sys.argv[2]
BATCH = FOLDER[:4]

cmd_download_1 = CMD_PREFIX + FOLDER + "/outs/" + H5_FILE + " ./"
cmd_download_2 = CMD_PREFIX + FOLDER + "/outs/kept_cells_DF.txt ./"

process1 = subprocess.Popen(cmd_download_1.split(), stdout=subprocess.PIPE)
output1, error1 = process1.communicate()

process2 = subprocess.Popen(cmd_download_2.split(), stdout=subprocess.PIPE)
output2, error2 = process2.communicate()

adata=sc.read_10x_h5(H5_FILE)
adata.var_names_make_unique()

kept_cells_DF=pd.read_csv('kept_cells_DF.txt',header=None)[0]
adata=adata[kept_cells_DF,:]

adata.obs.index=[BATCH+'_'+item for item in adata.obs.index]

adata.write(f'./raw_cellbender_DF/{BATCH}.h5ad')
