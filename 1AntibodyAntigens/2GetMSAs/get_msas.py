### IMPORTS AND STATIC PATHS
from pathlib import Path
import pdb
import pickle
import argparse
import subprocess
import sys
SRC = Path(__file__).parents[2].resolve() / "immunomsa"
sys.path.append( str(SRC) )
from general_functions import load_pickle_file

WORKDIR = Path(__file__).parent.resolve()
FASTADIR = WORKDIR / "fasta_files"
MSADIR = WORKDIR / "msas"
TMP_DIR = WORKDIR / "tmp"

UNIPROT_DB = Path(WORKDIR / "../../data/uniref30_db/uniref30_2302_db").resolve()
UNIPROT_DB_IDX = Path(WORKDIR / "../../data/uniref30_db/index_folder").resolve()
ANTIBODY_DB = Path(WORKDIR / "../1GetAntiBodyDatabase/antibodies/mmseqs_antibody_db/antibody_db").resolve()
ANTIBODY_DB_IDX = Path(WORKDIR / "../1GetAntiBodyDatabase/antibodies/mmseqs_antibody_db/index_folder").resolve()

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Run MMseqs MSAs")
parser.add_argument("-msa_run", required=True, action="store", dest="msa_run", help="ab_uniref, ag_uniref or ab_abdb.")
args = parser.parse_args()
msa_run = args.msa_run

### FUNCTIONS ###

def null_readlines(f, bufsize):
    buf = ""
    data = True
    while data:
        data = f.read(bufsize)
        buf += data
        lines = buf.split('\x00')
        buf = lines.pop()
        for line in lines:
            yield line
    yield buf

### MAIN ###

## Creating antigen + antibody seuqences (NOTE: This is case specific, use your own abag sequence data) ##

run = False
if run:
    cry_abags = Path("/home/projects2/joacl/crystal_data/antibody_antigen_complexes") 
    cry_abags = cry_abags.glob("*")
    cry_abag_annos = [abag / "anno.pickle" for abag in cry_abags]
    ag_seqs, ab_seqs = [], []
    for i, cry_abag_anno in enumerate(cry_abag_annos):

        # read data from annotation file
        pdbid = cry_abag_anno.parent.name.split("_")[0]
        anno = load_pickle_file(cry_abag_anno)
        abagtypes = anno["abagtypes"]
        chainids = anno["chainids"]
        sequences = anno["sequences"]
        
        # get indexes for antigen + antibody
        N = len(abagtypes)
        ag_idxs = [j for j in range(N) if abagtypes[j] == "Ag"]
        ab_idxs = [j for j in range(N) if abagtypes[j] in ("L", "H")]
        
        # collect sequences data
        ag_seqs.extend([(f">{pdbid}_ag_{chainids[idx]}\n{sequences[idx]}") for idx in ag_idxs])
        ab_seqs.extend([(f">{pdbid}_ab_{chainids[idx]}\n{sequences[idx]}") for idx in ab_idxs])
        print(f"Collected sequences from {i+1} antibody-antigen complexes: {cry_abag_anno}")
        
    # antigen write fasta file
    if not FASTADIR.is_dir(): FASTADIR.mkdir(parents=True)
    ag_seqs = "\n".join(ag_seqs)
    with open(FASTADIR / "antigens.fasta", "w") as outfile: outfile.write(ag_seqs)
    # write antibody fasta file
    ab_seqs = "\n".join(ab_seqs)
    with open(FASTADIR / "antibodies.fasta", "w") as outfile: outfile.write(ab_seqs)



def split_mmseqs_a3m_file(a3m_file, outdir):

    # get individual .a3m (they are split by null characters)
    with open(a3m_file, 'r') as f:
        a3ms = [item for item in null_readlines(f, 524288)]
    a3ms = [a3m for a3m in a3ms if a3m] # remove empty trail a3m

    # write individual .a3m files
    if not outdir.is_dir(): outdir.mkdir(parents=True)
    for a3m in a3ms:    
        a3m_filename = a3m.split("\n")[0][1:]
        outfile_path = outdir / f"{a3m_filename}.a3m"
        with open(outfile_path, "w") as outfile: outfile.write(a3m)


## Creating MSA inputs with MMseqs2 ##
        
def create_a3m_files(query_file, against_db, against_db_idx, queryout, combined_msas_outfile, msa_dir, cleanup_queryout=True):

    # create query output directory
    if not queryout.is_dir(): queryout.mkdir(parents=True)
    query_db = str(queryout / "query_db")
    query_file = str(query_file)

    # create query database
    subprocess.run(["mmseqs", "createdb", query_file, query_db])
    # run search 
    result_db = str(queryout / "result_db")
    subprocess.run(["mmseqs", "search", query_db, against_db, result_db, against_db_idx])
    # create msa file (.a3m format)
    subprocess.run(["mmseqs", "result2msa", query_db, against_db, result_db, combined_msas_outfile, "--msa-format-mode", "6"])
    # split into invididual .a3m files 
    split_mmseqs_a3m_file(msa_outfile, msa_dir)

    if cleanup_queryout:
        dirs_to_rm = []
        all_paths = queryout.glob("**/*")

        # remove all files 
        for p in all_paths:
            
            if p.is_file(): p.unlink()
            elif p.is_dir(): dirs_to_rm.append(p)

        # remove empty directories
        for d in dirs_to_rm: d.rmdir()
        
        # remove query outdir
        queryout.rmdir()
    
# create antibody .a3m files using uniref30 database 

if msa_run == "ab_uniref":
    against_db = str(UNIPROT_DB)
    against_db_idx = str(UNIPROT_DB_IDX)
    queryout = TMP_DIR / "antibodies_mmseqs_search_uniref"
    msa_outfile = queryout / "antibodies.a3m"
    outdir = MSADIR / "antibodies" / "uniref30"
    create_a3m_files(FASTADIR / "antibodies.fasta", against_db, against_db_idx, queryout, msa_outfile, outdir, cleanup_queryout=True)

# create antigen .a3m files using uniref30 database
elif msa_run == "ag_uniref":
    against_db = str(UNIPROT_DB)
    against_db_idx = str(UNIPROT_DB_IDX)
    queryout = TMP_DIR / "antigen_mmseqs_search_uniref"
    msa_outfile = queryout / "antigens.a3m"
    outdir = MSADIR / "antigens" / "uniref30"
    create_a3m_files(FASTADIR / "antigens.fasta", against_db, against_db_idx, queryout, msa_outfile, outdir, cleanup_queryout=True)

# create antibody .a3m files using antibody database 
elif msa_run == "ab_abdb":
    against_db = str(ANTIBODY_DB)
    against_db_idx = str(ANTIBODY_DB_IDX)
    queryout = TMP_DIR / "antibody_mmseqs_search_abdb"
    msa_outfile = queryout / "antibodies.a3m"
    outdir = MSADIR / "antibodies" / "antibody_db"
    create_a3m_files(FASTADIR / "antibodies.fasta", against_db, against_db_idx, queryout, msa_outfile, outdir, cleanup_queryout=True)


else:
    print(f"Invalid msa_run: {msa_run}. Please specify one of the following ab_uniref, ag_uniref or ab_abdb.")