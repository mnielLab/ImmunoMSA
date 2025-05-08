### IMPORTS AND STATIC STUFF ###
import subprocess
from pathlib import Path
import json
from joblib import delayed, Parallel
import gzip
import pandas as pd
import pdb
import sys
SRC = Path(__file__).parents[2].resolve() / "immunomsa"
sys.path.append( str(SRC) )
from pdb_and_cif_utilities import read_pdb_structure, write_pdb_chain_seq
from general_functions import split_list, read_accs_and_sequences_from_fasta

WORKDIR = Path(__file__).parent.resolve()
AB_DATADIR = WORKDIR / "antibodies"
SABDAB_DATADIR = AB_DATADIR / "SABDAB"
SABDAB_PDBS = SABDAB_DATADIR / "all_structures" / "raw"
OSA_HTTP_LINKS_FILE = AB_DATADIR / "oas_antibody_https.txt" # from https://opig.stats.ox.ac.uk/webapps/oas/oas_paired/ 
TMP = WORKDIR / "tmp"
num_workers= 10

#NOTE
# This scripts will extract download all antibody sequences form OSA database.
# Currently (02/05/2025) no automated way to download the the sequences SABDAB (download script no longer works) 
# Instead, download all structures + summary file from SABDAB and store them in SABDAB_DATADIR
# This script will then all antibody sequences.
# Finally, it combines all unique antibody sequnces, and creates an antibody mmseqs database

### FUNCTIONS ###

def get_osa_data(osa_http_links_file, heavy_sequences_outpath, light_sequences_outpath):
    osa_http_links = []
    infile = open(osa_http_links_file, "r" )
    for line in infile: osa_http_links.append(line.strip())
    infile.close()

    # collect all paired antibody sequences 
    datas = Parallel(n_jobs = num_workers)(delayed(download_and_parse_osa_csv)(osa_http_link, TMP) for osa_http_link in osa_http_links)
    heavy_seqs, light_seqs, = [], []
    for data in datas:
        heavy_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[2], data[0])] ) 
        light_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[3], data[1])] ) 
        
    heavy_seqs = "\n".join(heavy_seqs)
    light_seqs = "\n".join(light_seqs)

    # write sequences
    with open(heavy_sequences_outpath, "w") as outfile: outfile.write(heavy_seqs)
    with open(light_sequences_outpath, "w") as outfile: outfile.write(light_seqs)

def download_and_parse_osa_csv(link, tmpdir):

    heavy_seqs, light_seqs, light_headers, heavy_headers = [], [], [], []
    disease, species = "", ""

    # download + read file
    download_filename = link.split("/")[-1]
    outfile_path = tmpdir / download_filename
    subprocess.run(["wget", "-P", tmpdir, link])
    
    infile = gzip.open(outfile_path, "r") 
    file_lines = [l.decode("utf-8") for l in  infile.readlines()]  
    infile.close()

    # get header data
    header = file_lines[0]
    cleaned_header = header.strip().strip('"')
    json_str = cleaned_header.replace('""', '"')
    header_dict = json.loads(json_str)
    disease = header_dict["Disease"]
    species = header_dict["Species"]

    df = pd.read_csv(outfile_path, compression='gzip', skiprows=1)

    heavy_seqs = df["v_sequence_alignment_aa_heavy"].tolist()
    light_seqs = df["v_sequence_alignment_aa_light"].tolist()
    
    N = len(heavy_seqs)
 
    # get all the filenames
    # create accession headers 
    study_name = download_filename.split("_")[0]
    for i in range(N-2):
        heavy_headers.append(f">{study_name}_HEAVY_{i}_disease_{disease}_species_{species}")
        light_headers.append(f">{study_name}_LIGHT_{i}_disease_{disease}_species_{species}")
        
    # remove downloaded file
    outfile_path.unlink()

    return heavy_seqs, light_seqs, heavy_headers, light_headers, disease, species

def extract_light_heavy_chains_from_sabdab_pdbs(sabdab_summaryfile, pdb_filedir, heavy_sequences_outpath, light_sequences_outpath):

    df = pd.read_csv(sabdab_summaryfile, sep='\t')
    pdb_paths = []
    for pdb_acc in df["pdb"]:
        
        pdb_path = pdb_filedir / f"{pdb_acc}.pdb"
        if not pdb_path.is_file(): pdb.set_trace()
        pdb_paths.append(pdb_path)

    heavy_ids = df["Hchain"]
    light_ids = df["Lchain"]
    model_nums = df["model"]
    data = list( zip(pdb_paths, heavy_ids, light_ids, model_nums) )

    # split data into data lists + parrellize
    data_lists = split_list(data, num_workers)
    datas = Parallel(n_jobs = num_workers)(delayed(get_light_heavy_chains_from_sabdab_pdb_wrapper)(data_list) for data_list in data_lists)

    heavy_seqs, light_seqs = [], []
    for data in datas:
        heavy_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[2], data[0])] ) 
        light_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[3], data[1])] ) 
        
    
    heavy_seqs = "\n".join(heavy_seqs)
    light_seqs = "\n".join(light_seqs)

    # write sequences
    with open(heavy_sequences_outpath, "w") as outfile: outfile.write(heavy_seqs)
    with open(light_sequences_outpath, "w") as outfile: outfile.write(light_seqs)


def get_light_heavy_chains_from_sabdab_pdb_wrapper(data):
  
    # extract sequences
    heavy_seqs, light_seqs, heavy_headers, light_headers = [], [], [], []
    N = len(data)
    for i in range(N):
        d = data[i]
        pdb_path, heavy_id, light_id, model_num = d
        heavy_seq, light_seq, heavy_header, light_header = get_light_heavy_chains_from_sabdab_pdb(pdb_path, heavy_id, light_id, model_num)
       
        # append stuff
        if heavy_seq:
            heavy_seqs.append(heavy_seq)
            heavy_headers.append(heavy_header)
        
        if light_seq:
            light_seqs.append(light_seq)
            light_headers.append(light_header)

        print(f"Collect antibody sequence from {i}/{N}")


    return heavy_seqs, light_seqs, heavy_headers, light_headers


def get_light_heavy_chains_from_sabdab_pdb(pdb_file, heavy_id, light_id, model_num):

    # get load light + heavy chain    
    pdb_structure = read_pdb_structure(pdb_file, pdb_id="foo", modelnr=model_num, return_all_models = False)
    chainids = [c.get_id() for c in pdb_structure.get_chains()]

    # get heavy chain data
    heavy_seq, heavy_header = "", "" 
    if heavy_id in chainids:
        heavy = pdb_structure[heavy_id]
        heavy_seq = write_pdb_chain_seq( heavy )
        heavy_header = f">{pdb_file.stem}_HEAVY_{heavy_id}"

    else: print(f"Heavy chain ID {heavy_id} not found in {pdb_file.stem}")
 
    # get light chain data
    light_seq, light_header = "", ""
    if light_id in chainids:
        light = pdb_structure[light_id]
        light_seq = write_pdb_chain_seq( light )
        light_header = f">{pdb_file.stem}_LIGHT_{light_id}"

    else: print(f"Light chain ID {light_id} not found in {pdb_file.stem}")

   
    return heavy_seq, light_seq, heavy_header, light_header
    
### MAIN ###

## get OSA data
run = False
if run: get_osa_data(OSA_HTTP_LINKS_FILE, AB_DATADIR / "OSA_heavy_sequences.fasta", AB_DATADIR / "OSA_light_sequences.fasta")

## get SABDAB sequences ##
run = False
if run: extract_light_heavy_chains_from_sabdab_pdbs(SABDAB_DATADIR / "sabdab_summary_all.tsv", SABDAB_PDBS, AB_DATADIR / "SABDAB_heavy_sequences.fasta", AB_DATADIR / "SABDAB_light_sequences.fasta" )
    
# concatanate all antibody sequences 
run = False
if run:

    # 10x10 oas database antibody sequences
    with open( AB_DATADIR / "OSA_heavy_sequences.fasta") as infile: osa_heavy = infile.read()
    with open( AB_DATADIR / "OSA_light_sequences.fasta") as infile: osa_light = infile.read()

    # strutural antibody sequences 
    with open( AB_DATADIR / "SABDAB_heavy_sequences.fasta") as infile: sabdab_heavy = infile.read()
    with open( AB_DATADIR / "SABDAB_light_sequences.fasta") as infile: sabdab_light = infile.read()

    # write all antibody sequences to fasta
    all_antibody_sequences = "\n".join( [osa_heavy, osa_light, sabdab_heavy, sabdab_light] )
    with open(AB_DATADIR / "antibody_sequences.fasta", "w") as outfile: outfile.write(all_antibody_sequences)

    # read accession IDs and sequences
    accs, seqs = read_accs_and_sequences_from_fasta(AB_DATADIR / "antibody_sequences.fasta")

    # pair and sort by sequence to group duplicates
    sorted_acc_and_seqs = sorted(zip(accs, seqs), key=lambda x: x[1])

    # use a dict to collect the first accession per unique sequence
    seen = {}
    for acc, seq in sorted_acc_and_seqs:
        if seq not in seen:
            seen[seq] = acc

    uniq_acc_and_seqs = [f">{acc}\n{seq}" for seq, acc in seen.items()]
    print(f"Number of unique antibody sequences: {len(uniq_acc_and_seqs)}")
    
    # write unique antibody sequence ffasta file
    unique_antibody_sequences = "\n".join(uniq_acc_and_seqs)
    with open(AB_DATADIR / "antibody_sequences_uniq.fasta", "w") as outfile: outfile.write(unique_antibody_sequences)

# this step takes a very long time
run = False
if run:
    antibody_sequences_fastafile = str(AB_DATADIR / "antibody_sequences_uniq.fasta")
    cluster_results = str(AB_DATADIR / "ClusterRes")
    subprocess.run(["mmseqs", "easy-cluster", antibody_sequences_fastafile, cluster_results, TMP, "--min-seq-id", "0.99", "-c", "0.8", "--cov-mode", "1"])

# create mmseqs database
run = False
if run:

    # create mmseqs database directories
    mmmseqs_ab_database_dir = AB_DATADIR / "mmseqs_antibody_db"
    mmmseqs_ab_database_preidx = mmmseqs_ab_database_dir / "index_folder"
    if not mmmseqs_ab_database_preidx.is_dir(): mmmseqs_ab_database_preidx.mkdir(parents=True)

    # create mmseqs database
    antibody_sequences_fastafile = str(AB_DATADIR / "ClusterRes_rep_seq.fasta")
    subprocess.run(["mmseqs", "createdb", antibody_sequences_fastafile, str(mmmseqs_ab_database_dir / "antibody_db")])
    
    # create mmseqs preprocessed index directory
    if not mmmseqs_ab_database_preidx.is_dir(): mmmseqs_ab_database_preidx.mkdir(parents=True)
    # create mmseqs preprocessed index
    subprocess.run(["mmseqs", "createindex", str(mmmseqs_ab_database_dir / "antibody_db"), str(mmmseqs_ab_database_preidx)])