from Bio.PDB import PDBParser, PDBList
import numpy as np
import matplotlib.pyplot as plt
import os

def calculate_intraresidue_distances(structure):
    model = structure[0]  
    chain = model['A']  
    ca_atoms = []
    for residue in chain:
        if 'CA' in residue:
            ca_atoms.append(residue['CA'])
    
    num_residues = len(ca_atoms)
    
    distance_matrix = np.zeros((num_residues, num_residues))
    
    for i in range(num_residues):
        for j in range(num_residues):
            distance_matrix[i, j] = ca_atoms[i] - ca_atoms[j]
    
    return distance_matrix

def fetch_pdb_file(pdb_id):
    pdbl = PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
    return pdb_file_path


pdb_id = '1BTP' 

pdb_file = fetch_pdb_file(pdb_id)

parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdb_file)

distance_map = calculate_intraresidue_distances(structure)

plt.imshow(distance_map, cmap='viridis')
plt.colorbar(label='Distance (Å)')
plt.title(f'Intraresidue Distance Map (CA Atoms) for PDB ID: {pdb_id}')
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.show()

