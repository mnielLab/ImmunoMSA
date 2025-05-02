### IMPORTS  ###
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Selection

### FUNCTIONS ###

AA3to1_DICT = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

def is_pdb_file(file_path):
    """
    File extension needs to be pdb/PDB, as well as readable by PDB parser.
    """

    # check file ending first
    if Path(file_path).suffix != '.pdb': return False
    # check if can be parsed as pdb, if no exception is raised, likely valid pdb file
    try:
        parser = PDBParser(QUIET=True)  
        structure = parser.get_structure('pdb_structure', file_path)
        return True  
    except Exception: return False
  
def is_cif_file(file_path):
    """
    File extension needs to be cif/CIF, as well as readable by CIF parser.
    """

    # check file ending first
    if Path(file_path).suffix != '.cif': return False
    # check if can be parsed as cif, if no exception is raised, likely valid cif file
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('cif_structure', file_path)
        return True

    except Exception: return False


def read_pdb_structure(pdb_file, pdb_id="foo", modelnr=0, return_all_models = False):
        """
        pdb_id: PDB acession, string
        pdb_file: path to 

        """
        #reading model 0 by default
    
        assert isinstance(modelnr, int), f"Model number needs to be a valid integer, it was {modelnr}"
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_file)

        #return all models
        if return_all_models:
            models = list()
            for m in structure: models.append(m)
            return models

        #return only desired model
        else: return structure[modelnr]

def read_cif_structure(cif_file, pdb_id="foo", modelnr=0, return_all_models=False):
    """
    Reads a CIF file and returns a specific model or all models in the structure.
    
    cif_file: path to the CIF file.
    pdb_id: PDB accession, string (default is 'foo').
    modelnr: The model number to return (default is 0).
    return_all_models: If True, returns all models in the structure.
    
    Returns the specified model or all models if return_all_models is True.
    """
    assert isinstance(modelnr, int), f"Model number needs to be a valid integer, it was {modelnr}"
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id, cif_file)

    if return_all_models:
        models = list(structure.get_models())
        return models
    else:
        return structure[modelnr]
    

def cif_to_pdb(cif_file, pdb_file, verbose=True):
    """
    Converts a CIF file to a PDB file using Biopython.

    Parameters:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to the output PDB file.

    Returns:
        None
    """
    try:
        # Parse the CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", cif_file)

        # Write the structure to a PDB file
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)

        if verbose: print(f"Successfully converted {cif_file} to {pdb_file}.")
    except Exception as e:
        print(f"An error occurred: {e}")


def get_and_remove_heteroatoms(chain):
    """
   Heteroatoms in the form of water and other solvents need to be removed from the chain.
   Inputs: chain id

    """
    residues = Selection.unfold_entities(chain, "R")
    heteroatom_residue_ids = list()
    for residue in residues:
        residue_id = residue.get_full_id()

        #residue is a heteroatom
        if residue_id[3][0] != " ":
            heteroatom_residue_ids.append(residue_id[3])
    #remove all heteroatoms
    [chain.detach_child(ids) for ids in heteroatom_residue_ids]

def write_pdb_res_to_seq(residues):
    """
    residues: Bio PDB residues
    """
    AA_seq = str()
    get_and_remove_heteroatoms(residues)

    for residue in residues:
        try:
            aa = AA3to1_DICT[residue.get_resname()]
        #when residue is something nonstandard
        except KeyError:
            print(f"Non-standard amino acid detected, {residue.get_resname()}")
            aa = "X"

        AA_seq += aa

    return AA_seq

def write_pdb_chain_seq(chain):
    """
    residues: Bio PDB residues
    """
    AA_seq = str()
    # remove hetero atoms
    get_and_remove_heteroatoms(chain)
    residues = list( chain.get_residues() )

    for residue in residues:
        try:
            aa = AA3to1_DICT[residue.get_resname()]
        #when residue is something nonstandard
        except KeyError:
            print(f"Non-standard amino acid detected, {residue.get_resname()}")
            aa = "X"

        AA_seq += aa

    return AA_seq
