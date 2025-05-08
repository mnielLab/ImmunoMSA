### IMPORTS AND STATIC STUFF ###
from pathlib import Path
import subprocess
import pdb
import sys

WORKDIR = Path(__file__).parent.resolve()
uniprot30_db = WORKDIR / "uniref30_db"
UNIREF_DATA = "uniref30_2302"
mmseqs_colabfold_uniprot30_dlink = f"https://wwwuser.gwdg.de/~compbiol/colabfold/{UNIREF_DATA}.tar.gz"

### MAIN ###
# download mmseqs uniref 30
run = False
if run: subprocess.run(["wget", "-P", str(uniprot30_db), mmseqs_colabfold_uniprot30_dlink])
# extract gunzipped tar archive 
run = False
if run:
    uniprot30_db_tar = str(uniprot30_db / f"{UNIREF_DATA}.tar.gz")
    subprocess.run(["tar", "xzf", uniprot30_db_tar, "-C", str(uniprot30_db)] )
    Path(uniprot30_db_tar).unlink()

# # create mmseqs2 unired30 database
run = True
if run:
    uniprot30_db_tsv =  str(uniprot30_db / UNIREF_DATA)
    uniprot30_db_out = str(uniprot30_db / f"{UNIREF_DATA}_db")
    subprocess.run(["mmseqs", "tsv2exprofiledb", uniprot30_db_tsv, uniprot30_db_out])
          
# create mmseqs2 indexing 
run = True
if run:
    mmseqs_index_folder = uniprot30_db / "index_folder"
    if not mmseqs_index_folder.is_dir(): mmseqs_index_folder.mkdir(parents=True)
    mmseqs_index_folder = str(mmseqs_index_folder)
    subprocess.run(["mmseqs", "createindex", uniprot30_db_out, mmseqs_index_folder])