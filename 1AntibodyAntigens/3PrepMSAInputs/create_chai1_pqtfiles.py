### IMPORTS AND STATIC STUFF ###
from pathlib import Path
import shutil
import pandas as pd
import sys
import pdb
from chai_lab.data.parsing.msas.aligned_pqt import merge_a3m_in_directory, expected_basename
SRC = Path(__file__).parents[2].resolve() / "immunomsa"
sys.path.append( str(SRC) )
from general_functions import read_accs_and_sequences_from_fasta

WORKDIR = Path(__file__).parents[0].resolve()
TMPDIR = WORKDIR / "tmp"

# msa directories 
MSADIR = Path(WORKDIR / "../2GetMSAs/msas").resolve()
AB_MSADIR = MSADIR / "antibodies"
AB_MSADIR_UNIREF = AB_MSADIR / "uniref30"
AB_MSADIR_ABDB = AB_MSADIR  / "antibody_db"
AG_MSADIR = MSADIR / "antigens"
AG_MSADIR_UNIREF = AG_MSADIR / "uniref30"

# pqt directories
PQTDIR = WORKDIR / "chai1_pqt"
AB_UNIREF_AG_UNIREF = PQTDIR / "abuniref30_aguniref30"
AB_ABDB_AG_UNIREF = PQTDIR / "ababdb_aguniref30"

# paired pqt directories
AB_UNIREF_AG_UNIREF_abpaired = PQTDIR / "abuniref30_aguniref30_abpaired"
AB_UNIREF_AG_UNIREF_agpaired = PQTDIR / "abuniref30_aguniref30_agpaired"
AB_UNIREF_AG_UNIREF_abpaired_agpaired = PQTDIR / "abuniref30_aguniref30_abpaired_agpaired"
AB_UNIREF_AG_UNIREF_allpaired = PQTDIR / "abuniref30_aguniref30_allpaired"

AB_ABDB_AG_UNIREF_abpaired = PQTDIR / "ababdb_aguniref30_abpaired"
AB_ABDB_AG_UNIREF_agpaired = PQTDIR / "ababdb_aguniref30_agpaired"
AB_ABDB_AG_UNIREF_abpaired_agpaired = PQTDIR / "ababdb_aguniref30_abpaired_agpaired"
AB_ABDB_AG_UNIREF_allpaired = PQTDIR / "ababdb_aguniref30_allpaired"


### FUNCTIONS ###
def create_pqtfiles_from_a3mfiles(a3m_files, outdir, tmpdir, database = "uniref90"):

    N = len(a3m_files)
    # copy everything,  
    for i in range(N):
        a3m_file = a3m_files[i]
            
        # copy .a3m file 
        fromfile = a3m_file
        a3m_tmpdir = tmpdir / a3m_file.stem
        tofile = a3m_tmpdir / f"hits_{database}.a3m"
        if not tofile.parent.is_dir(): tofile.parent.mkdir(parents=True)
        shutil.copy(fromfile, tofile)
        
        # check hash 
        _, seqs = read_accs_and_sequences_from_fasta(a3m_file)
        query_seq = seqs[0]
        checkfile = outdir / f"{expected_basename(query_seq)}"
        if checkfile.is_file():
            print(f"Already found hash for query seqeunce {a3m_file.stem}. This just means that a .pqt for this sequence is already in {outdir}. It will be overwritten {outdir}")
            aligned_pqt = pd.read_parquet(checkfile)
            aligned_pqt.head()

        # convert to .pqt
        merge_a3m_in_directory(a3m_tmpdir, output_directory=outdir)

        # remove temporary stuff
        tofile.unlink()
        a3m_tmpdir.rmdir()
        print(f"Created pqt files for {i+1}/{N}")


def create_paired_pqt_files(pqt_files, outdir, ab_pairkey=None, ag_pairkey=None, abag_pairkey=None):

    if not outdir.is_dir(): outdir.mkdir(parents=True)
    N = len(pqt_files)
    for i in range(N):
        pqt_file = pqt_files[i]
        df = pd.read_parquet(pqt_file)

        query = df.loc[0]
        query_name = query.comment

        if abag_pairkey is not None: df["pairing_key"] = "AbAg"
        if ab_pairkey is not None and ab_pairkey in query_name: df["pairing_key"] = "Ab"
        if ag_pairkey is not None and ag_pairkey in query_name: df["pairing_key"] = "Ag"
            
        outfile_path = outdir / pqt_file.name
        df.to_parquet(outfile_path)

        print(f"Created {i}/{N} paired .pqt files: {outfile_path}")


### MAIN ###

# create directories 
if not TMPDIR.is_dir(): TMPDIR.mkdir(parents=True) 
if not AB_UNIREF_AG_UNIREF.is_dir(): AB_UNIREF_AG_UNIREF.mkdir(parents=True) 
if not AB_ABDB_AG_UNIREF.is_dir(): AB_ABDB_AG_UNIREF .mkdir(parents=True) 

## UNIREF30 ##
run = False
if run:
    create_pqtfiles_from_a3mfiles(list(AB_MSADIR_UNIREF.glob("*.a3m")), AB_UNIREF_AG_UNIREF, TMPDIR, database = "uniref90")
    create_pqtfiles_from_a3mfiles(list(AG_MSADIR_UNIREF.glob("*.a3m")), AB_UNIREF_AG_UNIREF, TMPDIR, database = "uniref90")

##  AbDb (Ab) UNIREF30 (Ag) ##
run = False
if run:
    create_pqtfiles_from_a3mfiles(list(AB_MSADIR_ABDB.glob("*.a3m")), AB_ABDB_AG_UNIREF, TMPDIR, database = "uniref90")
    create_pqtfiles_from_a3mfiles(list(AG_MSADIR_UNIREF.glob("*.a3m")), AB_ABDB_AG_UNIREF, TMPDIR, database = "uniref90")

## UNIREF30 ##
run = False
pqt_files = list(AB_UNIREF_AG_UNIREF.glob("*.pqt"))
if run:    
    create_paired_pqt_files(pqt_files, AB_UNIREF_AG_UNIREF_abpaired, ab_pairkey="_ab_")
    create_paired_pqt_files(pqt_files, AB_UNIREF_AG_UNIREF_agpaired, ag_pairkey="_ag_")
    create_paired_pqt_files(pqt_files, AB_UNIREF_AG_UNIREF_abpaired_agpaired, ab_pairkey="_ab_", ag_pairkey="_ag_")
    create_paired_pqt_files(pqt_files, AB_UNIREF_AG_UNIREF_allpaired, abag_pairkey="yes")

## AbDb (Ab) UNIREF30 (Ag) Paired PQT ##
run = False
pqt_files = list(AB_ABDB_AG_UNIREF.glob("*.pqt"))
if run:
    create_paired_pqt_files(pqt_files, AB_ABDB_AG_UNIREF_abpaired, ab_pairkey="_ab_")
    create_paired_pqt_files(pqt_files, AB_ABDB_AG_UNIREF_agpaired, ag_pairkey="_ag_")
    create_paired_pqt_files(pqt_files, AB_ABDB_AG_UNIREF_abpaired_agpaired, ab_pairkey="_ab_", ag_pairkey="_ag_")
    create_paired_pqt_files(pqt_files, AB_ABDB_AG_UNIREF_allpaired, abag_pairkey="yes")