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
from pdb_and_cif_utilities import read_pdb_structure, write_pdb_chain_seq, AA3to1_DICT
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



def write_fasta(filename, accs, seqs):
    with open(filename, "w") as f:
        for acc, seq in zip(accs, seqs):
            f.write(f">{acc}\n{seq}\n")

def split_into_fasta_files(accs, seqs, output_dir, num_files=10):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    total = len(accs)
    chunk_size = (total + num_files - 1) // num_files  # Ceiling division

    for i in range(num_files):
        start = i * chunk_size
        end = min(start + chunk_size, total)
        chunk_accs = accs[start:end]
        chunk_seqs = seqs[start:end]
        output_path = output_dir / f"split_{i+1}.fasta"
        write_fasta(output_path, chunk_accs, chunk_seqs)

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

def extract_heavy_light_chains_from_sabdab_pdbs(sabdab_summaryfile, pdb_filedir, outdir):

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
    datas = Parallel(n_jobs = num_workers)(delayed(get_heavy_light_chains_from_sabdab_pdb_wrapper)(data_list) for data_list in data_lists)

    # paired heavy and light chains ()
    paired_heavy_seqs, paired_light_seqs = [], []
    single_heavy_seqs, single_light_seqs = [], []
    for data in datas:
        paired_heavy_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[0], data[1])] ) 
        paired_light_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[2], data[3])] ) 
        single_heavy_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[4], data[5])] ) 
        single_light_seqs.extend( [f"{acc}\n{seq}" for acc, seq in zip(data[6], data[7])] ) 

    # write paired sequence files
    heavy_seqs = "\n".join(paired_heavy_seqs)
    light_seqs = "\n".join(paired_light_seqs)
    heavy_sequences_outpath = outdir / "SABDAB_paired_heavy_sequences.fasta"
    light_sequences_outpath = outdir / "SABDAB_paired_light_sequences.fasta"
    with open(heavy_sequences_outpath, "w") as outfile: outfile.write(heavy_seqs)
    with open(light_sequences_outpath, "w") as outfile: outfile.write(light_seqs)
    
    # write single heavy or light chains sequence files
    heavy_seqs = "\n".join(single_heavy_seqs)
    light_seqs = "\n".join(single_light_seqs)
    heavy_sequences_outpath = outdir / "SABDAB_single_heavy_sequences.fasta"
    light_sequences_outpath = outdir / "SABDAB_single_light_sequences.fasta"
    with open(heavy_sequences_outpath, "w") as outfile: outfile.write(heavy_seqs)
    with open(light_sequences_outpath, "w") as outfile: outfile.write(light_seqs)


def get_heavy_light_chains_from_sabdab_pdb_wrapper(data):
  
    # extract sequences
    paired_heavy_headers, paired_heavy_seqs, paired_light_headers, paired_light_seqs = [], [], [], []
    single_heavy_headers, single_heavy_seqs, single_light_headers, single_light_seqs = [], [], [], []
   
    
    
    N = len(data)
    for i in range(N):
        d = data[i]
        pdb_path, heavy_id, light_id, model_num = d
        heavy_seq, light_seq, heavy_header, light_header = get_heavy_light_chains_from_sabdab_pdb(pdb_path, heavy_id, light_id, model_num)

        # paired antibody sequences (both heavy and light chain were found)
        if heavy_seq and light_seq:
            paired_heavy_headers.append(heavy_header)
            paired_heavy_seqs.append(heavy_seq)
            paired_light_headers.append(light_header)
            paired_light_seqs.append(light_seq)
            
        # only heavy chain was found
        if heavy_seq and not light_seq:
            single_heavy_headers.append(heavy_header)
            single_heavy_seqs.append(heavy_seq)
           
        # only light chain was found
        if light_seq and not heavy_seq:
            single_light_headers.append(light_header)
            single_light_seqs.append(light_seq)
            
  
        print(f"Collect antibody sequence from {i+1}/{N}")

    return (paired_heavy_headers, paired_heavy_seqs, paired_light_headers, paired_light_seqs,
            single_heavy_headers, single_heavy_seqs, single_light_headers, single_light_seqs)


def get_heavy_light_chains_from_sabdab_pdb(pdb_file, heavy_id, light_id, model_num):

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
if run: extract_heavy_light_chains_from_sabdab_pdbs(SABDAB_DATADIR / "sabdab_summary_all.tsv", SABDAB_PDBS, AB_DATADIR)


# concatanate all antibody sequences 
run = False
if run:
    
    # 10x10 oas database + sabdab structural paired antibody sequences 
    osa_heavy_accs, osa_heavy_seqs = read_accs_and_sequences_from_fasta(AB_DATADIR / "OSA_heavy_sequences.fasta")
    osa_light_accs, osa_light_seqs = read_accs_and_sequences_from_fasta(AB_DATADIR / "OSA_light_sequences.fasta")
    
    sabdab_heavy_accs, sabdab_heavy_seqs = read_accs_and_sequences_from_fasta(AB_DATADIR / "SABDAB_paired_heavy_sequences.fasta")
    sabdab_light_accs, sabdab_light_seqs = read_accs_and_sequences_from_fasta(AB_DATADIR / "SABDAB_paired_light_sequences.fasta")
    
    # concatenate sequences
    N1, N2 = len(osa_heavy_accs), len(sabdab_heavy_accs)
    heavy_lights = ["".join([osa_heavy_seqs[i], osa_light_seqs[i]]) for i in range(N1)] + ["".join([sabdab_heavy_seqs[i], sabdab_light_seqs[i]]) for i in range(N2)]    
    heavy_accs = osa_heavy_accs + sabdab_heavy_accs
    heavy_seqs = osa_heavy_seqs + sabdab_heavy_seqs
    light_accs = osa_light_accs + sabdab_light_accs
    light_seqs = osa_light_seqs + sabdab_light_seqs
    sorted_acc_and_seqs = sorted(zip(heavy_accs, heavy_seqs, light_accs, light_seqs, heavy_lights), key=lambda x: x[-1] ) 
    
    print(f"Paired antibody sequences: {len(sorted_acc_and_seqs)}")

    unique_heavy, unique_light, unique_lightheavy = [], [], [] 
    seen = set()
    N = len(sorted_acc_and_seqs)
    for i in range(N):
        heavy_acc, heavy_seq, light_acc, light_seq, heavy_light = sorted_acc_and_seqs[i]
        if heavy_light not in seen:
            unique_heavy.append( f">{heavy_acc}_{i}\n{heavy_seq}" )
            unique_light.append( f">{light_acc}_{i}\n{light_seq}" )
            heavy_light_acc = "split".join([heavy_acc, light_acc])
            unique_lightheavy.append(f">heavylight_{i}\n{heavy_light}")
            seen.update( {heavy_light} )
            

    with open(AB_DATADIR / "paired_heavy_sequences.fasta", "w") as outfile: outfile.write("\n".join(unique_heavy) )
    with open(AB_DATADIR / "paired_light_sequences.fasta", "w") as outfile: outfile.write("\n".join(unique_light) ) 
    with open(AB_DATADIR / "heavylight_sequences.fasta", "w") as outfile: outfile.write("\n".join(unique_lightheavy) ) 
    

    print(f"Paired antibody sequences after removing duplicates (antibodies with identical heavy + light chain): {len(seen)}")

def run_mmmseqs_linclust(antibody_sequences_fastafile, linclust_out, cluster_outfile, tmp, seqid="0.99", coverage="0.8", covmode="1"):

    # create database
    if not linclust_out.is_dir(): linclust_out.mkdir(parents=True)
    input_db_path = linclust_out  / "input_db"
    subprocess.run(["mmseqs", "createdb", str(antibody_sequences_fastafile), str(input_db_path)])

    # run mmmseqs linclust
    result_db_path = linclust_out  / "result_db"
    subprocess.run(["mmseqs", "linclust", str(input_db_path), str(result_db_path), str(tmp), "--min-seq-id", seqid, "-c", coverage, "--cov-mode", covmode] )

    # extract the clusters
    clus_db_path =  linclust_out  / "clustered_db"
    subprocess.run(["mmseqs", "createseqfiledb", str(input_db_path), str(result_db_path), str(clus_db_path)])
    subprocess.run(["mmseqs", "result2flat", str(input_db_path), str(input_db_path), str(clus_db_path), str(cluster_outfile)])

    
def extract_mmseq_cluster_representatives(mmseqs_clusterfile):


    infile = open(mmseqs_clusterfile)
    on_acc_line = False
    cluster_representatives = []

    for line in infile:

        if line.startswith(">"):
            acc = line
            
            # two consecutive acc lines = cluster rep
            if on_acc_line:
                cluster_representatives.append(acc.strip())
                on_acc_line = False

            else: on_acc_line = True


    return cluster_representatives


def write_paired_clustered_heavy_light(all_paired_heavy, all_paired_light, mmseqs_clusterfile, heavy_outfile, light_outfile, antibody_outfile):
    
    # read all heavy + light paired chains
    heavy_accs, heavy_seqs = read_accs_and_sequences_from_fasta(all_paired_heavy)
    light_accs, light_seqs = read_accs_and_sequences_from_fasta(all_paired_light)
    N = len(heavy_accs)

    cluster_representatives = extract_mmseq_cluster_representatives(mmseqs_clusterfile)
    identifiers = set( [acc.split("_")[1] for acc in cluster_representatives] )
    paired_heavy, paired_light, paired_antibody = [], [], []

    for i in range(N):
        heavy_acc, heavy_seq = heavy_accs[i], heavy_seqs[i]
        light_acc, light_seq = light_accs[i], light_seqs[i] 
        
        identifier = heavy_acc.split("_")[-1]
 
        if identifier in identifiers:
            paired_antibody.extend( [f">{heavy_acc}_{i}\n{heavy_seq}", f">{light_acc}_{i}\n{light_seq}"] )
            paired_heavy.append( f">{heavy_acc}_{i}\n{heavy_seq}" )
            paired_light.append( f">{light_acc}_{i}\n{light_seq}" )

    # write clustered paired heavy + light pairs
    with open(heavy_outfile, "w") as outfile: outfile.write("\n".join(paired_heavy) )
    with open(light_outfile, "w") as outfile: outfile.write("\n".join(paired_light) )

    with open(antibody_outfile, "w") as outfile: outfile.write("\n".join(paired_antibody) )
    

def create_mmseqs_database(mmmseqs_ab_database_dir, antibody_sequences_fastafile):
    
    # create directories
    mmmseqs_ab_database_preidx = mmmseqs_ab_database_dir / "index_folder"
    if not mmmseqs_ab_database_preidx.is_dir(): mmmseqs_ab_database_preidx.mkdir(parents=True)

    # create mmseqs database
    subprocess.run(["mmseqs", "createdb", antibody_sequences_fastafile, str(mmmseqs_ab_database_dir / "antibody_db")])    
    # create mmseqs preprocessed index directory
    if not mmmseqs_ab_database_preidx.is_dir(): mmmseqs_ab_database_preidx.mkdir(parents=True)
    # create mmseqs preprocessed index
    subprocess.run(["mmseqs", "createindex", str(mmmseqs_ab_database_dir / "antibody_db"), str(mmmseqs_ab_database_preidx)])

run = False
if run:
    
    # cluster conctatenated heavy+light paired chains
    antibody_sequences_fastafile = str(AB_DATADIR / "heavylight_sequences.fasta")
    run_mmmseqs_linclust(antibody_sequences_fastafile, TMP / "linclust99", TMP / "linclust99" / "clustered.fasta", TMP, seqid="0.99")
    run_mmmseqs_linclust(antibody_sequences_fastafile, TMP / "linclust95", TMP / "linclust95" / "clustered.fasta", TMP, seqid="0.95")
    run_mmmseqs_linclust(antibody_sequences_fastafile, TMP / "linclust90", TMP / "linclust90" / "clustered.fasta", TMP, seqid="0.9")
    
    # write antibody cluster represenatives 
    write_paired_clustered_heavy_light(AB_DATADIR / "paired_heavy_sequences.fasta", AB_DATADIR / "paired_light_sequences.fasta", TMP / "linclust99" / "clustered.fasta",
                                       AB_DATADIR / "paired_heavy_sequences99.fasta", AB_DATADIR / "paired_light_sequences99.fasta", AB_DATADIR / "paired_antibody_sequences99.fasta")
    write_paired_clustered_heavy_light(AB_DATADIR / "paired_heavy_sequences.fasta", AB_DATADIR / "paired_light_sequences.fasta", TMP / "linclust95" / "clustered.fasta",
                                       AB_DATADIR / "paired_heavy_sequences95.fasta", AB_DATADIR / "paired_light_sequences95.fasta", AB_DATADIR / "paired_antibody_sequences95.fasta")
    write_paired_clustered_heavy_light(AB_DATADIR / "paired_heavy_sequences.fasta", AB_DATADIR / "paired_light_sequences.fasta", TMP / "linclust90" / "clustered.fasta",
                                       AB_DATADIR / "paired_heavy_sequences90.fasta", AB_DATADIR / "paired_light_sequences90.fasta", AB_DATADIR / "paired_antibody_sequences90.fasta")

    # clean up temporary files
    for f in Path(TMP / "linclust99").glob("*"): f.unlink()
    for f in Path(TMP / "linclust95").glob("*"): f.unlink()
    for f in Path(TMP / "linclust90").glob("*"): f.unlink()
    
    # remove temporary directories
    Path(TMP / "linclust99").rmdir()
    Path(TMP / "linclust95").rmdir()
    Path(TMP / "linclust90").rmdir()    

# create mmseqs database
run = False
if run:
    create_mmseqs_database(AB_DATADIR / "mmseqs_antibody_db99", AB_DATADIR / "paired_antibody_sequences99.fasta")
    create_mmseqs_database(AB_DATADIR / "mmseqs_antibody_db95", AB_DATADIR / "paired_antibody_sequences95.fasta")
    create_mmseqs_database(AB_DATADIR / "mmseqs_antibody_db90", AB_DATADIR / "paired_antibody_sequences90.fasta")
