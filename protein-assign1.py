import os
from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer, PDBList, PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Levenshtein import distance as levenshtein_distance
import pandas as pd

def lev_sequence_distance(seq1, seq2):
    return levenshtein_distance(seq1, seq2)

def per_id_sequence_distance(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    
    best_alignment = alignments[0]
    
    matches = sum(1 for a, b in zip(best_alignment.seqA, best_alignment.seqB) if a == b)
    
    pid = matches / max(len(seq1), len(seq2)) * 100
    
    return pid

def get_sequence_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    ppb = PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(structure):
        sequence += pp.get_sequence()
    return str(sequence)

def structure_distance(pdb_id1, pdb_id2):
    pdbl = PDBList()
    pdb_file1 = pdbl.retrieve_pdb_file(pdb_id1, file_format='pdb')
    pdb_file2 = pdbl.retrieve_pdb_file(pdb_id2, file_format='pdb')
    
    parser = PDBParser(QUIET=True)
    
    structure1 = parser.get_structure('structure1', pdb_file1)
    structure2 = parser.get_structure('structure2', pdb_file2)

    atoms1 = [atom for atom in structure1.get_atoms() if atom.get_id() == 'CA']
    atoms2 = [atom for atom in structure2.get_atoms() if atom.get_id() == 'CA']
    
    min_length = min(len(atoms1), len(atoms2))
    atoms1 = atoms1[:min_length]
    atoms2 = atoms2[:min_length]
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    
    super_imposer.apply(structure2.get_atoms())
    rmsd = super_imposer.rms

    return rmsd


# Base PDB ID and lists to compare
base_pdb_id = "1BTP"  # Trypsin PDB ID selected
sequence_search_pdb_ids = ["1BTW", "1BTX", "1BTY", "1BTZ", "1F0T"]
structure_search_pdb_ids = ["1O3L", "1O30", "1O2Q", "1G3D", "1XUG"]
random_pdb_ids = ["2C2W", "4CQJ", "1Z7L", "1EZA", "2EZB"]

sequence_results_dict = {}
structure_results_dict = {}

pdbl = PDBList()
base_pdb_file = pdbl.retrieve_pdb_file(base_pdb_id, file_format='pdb')
base_sequence = get_sequence_from_pdb(base_pdb_file)

def calculate_distances(base_pdb_id, pdb_ids_list, list_name):
    pdb_ids = []
    lev_dists = []
    percent_identities = []
    structure_dists = []
    
    for pdb_id in pdb_ids_list:
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
        sequence = get_sequence_from_pdb(pdb_file)
        
        lev_dist = lev_sequence_distance(base_sequence, sequence)
        percent_identity = per_id_sequence_distance(base_sequence, sequence)
        
        struct_dist = structure_distance(base_pdb_id, pdb_id)
        
        pdb_ids.append(pdb_id)
        lev_dists.append(lev_dist)
        percent_identities.append(percent_identity)
        structure_dists.append(struct_dist)
        
    
    results_df = pd.DataFrame({
        "PDB ID": pdb_ids,
        "Levenshtein Distance": lev_dists,
        "Percentage Identity": percent_identities,
        "Structure Distance (RMSD)": structure_dists
    })
    
    try:
        from IPython.display import display
        display(results_df)
    except ImportError:
        print(f"{list_name} Comparison Results:\n", results_df)

# Calculate distances for all lists
calculate_distances(base_pdb_id, sequence_search_pdb_ids, "Sequence Search PDB IDs")
calculate_distances(base_pdb_id, structure_search_pdb_ids, "Structure Search PDB IDs")
calculate_distances(base_pdb_id, random_pdb_ids, "Random PDB IDs")


