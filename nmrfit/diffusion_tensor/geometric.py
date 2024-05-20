from Bio.PDB import PDBParser
import numpy as np
import pandas as pd


## Functions for inertia tensor calculations

def get_coordinates_from_pdb(pdb_file, isotopes):
    masses = {'12C': 12.0, '13C': 13.0033548378, '14N': 14.0030740048, '15N': 15.0001088984}
    # Load PDB file and extract coordinates
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", pdb_file)

    atoms = []

    # Collect atom coordinates and calculate total mass
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    number = atom.get_full_id()[3][1]
                    aa = residue.resname
                    position = atom.get_coord()
                    mass = atom.mass
                    atom_str = str(atom).replace("<Atom ","")
                    atom_str = str(atom_str.replace(">",""))
                    if '15N' in isotopes:
                        if atom.element == 'N':
                            mass = masses['15N']
                    if '13C' in isotopes:
                        if atom.element == 'C':
                            mass = masses['13C']
                    atoms.append((number, aa, atom_str, position[0], position[1], position[2], mass))
    return pd.DataFrame(atoms, columns=["num", "AA", "atom", "x", "y", "z", "mass"])

def calculate_center_of_mass(coordinates):
    total_mass = coordinates['mass'].sum()
    
    center_x = (coordinates['x'] * coordinates['mass']).sum() / total_mass
    center_y = (coordinates['y'] * coordinates['mass']).sum() / total_mass
    center_z = (coordinates['z'] * coordinates['mass']).sum() / total_mass
    
    return center_x, center_y, center_z

def calculate_relative_coordinates(coordinates):
    com = calculate_center_of_mass(coordinates)
    coordinates[["x0", "y0", "z0"]] = coordinates[["x", "y", "z"]] - com
    # Normalized coordinates in case
    # norms = np.linalg.norm(coordinates[["x0", "y0", "z0"]].values, axis=1)
    # coordinates[["xn", "yn", "zn"]] = coordinates[['x0', 'y0', 'z0']].div(norms, axis=0)
    return coordinates

def calculate_inertia_tensor(normalized_coord): # Eigenvalues and Eigenvectors are given as z,y,x while the Inertia Tensor is Ixx, Iyy, Izz
    data = normalized_coord
    Ixx = sum(data["mass"]*((data["y0"]**2)+data["z0"]**2))
    Iyy = sum(data["mass"]*((data["x0"]**2)+data["z0"]**2))
    Izz = sum(data["mass"]*((data["x0"]**2)+data["y0"]**2))
    Ixy = -sum(data["mass"]*data["x0"]*data["y0"])
    Iyz = -sum(data["mass"]*data["y0"]*data["z0"])
    Ixz = -sum(data["mass"]*data["x0"]*data["z0"])
    
    inertia_tensor = np.array([[Ixx, Ixy, Ixz],
                               [Ixy, Iyy, Iyz],
                               [Ixz, Iyz, Izz]])
    
    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    order = np.argsort(eigenvalues)
    
    eve3, eve2, eve1 = eigenvectors[:, order].transpose()
    eva3, eva2, eva1 = eigenvalues[order]
    
    eigenvectors = np.array([eve3, eve2, eve1])
    eigenvalues = np.array([eva3, eva2, eva1])
    
    return inertia_tensor, eigenvectors, eigenvalues

def calculate_inertia_tensor_from_pdb(pdb_file, isotopes):
    data = get_coordinates_from_pdb(pdb_file, isotopes)
    data = calculate_relative_coordinates(data)
    inertia_tensor, eigenvectors, eigenvalues = calculate_inertia_tensor(data)
    return data, inertia_tensor, eigenvectors, eigenvalues

### Functions for NH_bond_vector calculation
def calculate_NH_bond_vectors(coordinates):
    nh_bond_vectors = []
    # nh_bond_vectors0 = []
    nh_direction_cos = []
    unique_residues = coordinates['num'].unique()

    for residue in unique_residues:
        residue_data = coordinates[coordinates['num'] == residue]
        
        n_coords = residue_data[residue_data['atom'] == 'N'][['x', 'y', 'z']].values
        h_coords = residue_data[residue_data['atom'] == 'HN'][['x', 'y', 'z']].values

        # n_coords0 = residue_data[residue_data['atom'] == 'N'][['x0', 'y0', 'z0']].values
        # h_coords0 = residue_data[residue_data['atom'] == 'HN'][['x0', 'y0', 'z0']].values

        if len(n_coords) == 1 and len(h_coords) == 1:
            nh_vector = h_coords - n_coords
            nh_unit_vector = nh_vector / np.linalg.norm(nh_vector)
            nh_bond_vectors.append(nh_unit_vector.squeeze())
            
            nh = nh_vector.squeeze()
            nn = nh[0]**2 + nh[1]**2 + nh[2]**2
            nh_direction_cos.append((nh_unit_vector/(np.sqrt(nn))).squeeze())
            
        # if len(n_coords0) == 1 and len(h_coords0) == 1:
        #     nh_vector0 = h_coords0 - n_coords0
        #     # nh_unit_vector0 = nh_vector0 / np.linalg.norm(nh_vector0)
        #     nh_bond_vectors0.append(nh_vector0.squeeze())

    df = pd.DataFrame(nh_bond_vectors, columns=["NH_x", "NH_y", "NH_z"])
    df2 = pd.DataFrame(nh_direction_cos, columns=["NH_l", "NH_m", "NH_n"])
    # df2 = pd.DataFrame(nh_bond_vectors0, columns=["NH_x0", "NH_y0", "NH_z0"])
    data = pd.concat([df, df2], axis=1)
    data.insert(0, "num", coordinates['num'].unique())
    return data

def get_NH_bond_vectors(pdb, isotopes, offset, relaxation_data):
    coordinates, intens, eve, eva = calculate_inertia_tensor_from_pdb(pdb, isotopes=isotopes)
    coordinates["num"] = coordinates["num"] + offset
    coordinates = coordinates[coordinates['num'].isin(relaxation_data['num'])]
    NH_bond_vectors = calculate_NH_bond_vectors(coordinates)    
    # return NH_bond_vectors, coordinates, intens, eve, eva
    return NH_bond_vectors

# Functions for NH_bond_vector_rotation
def rotate_NH_bond_vector(NH_bond_vectors, angles): #Woessner
      xyz = NH_bond_vectors[["NH_x","NH_y","NH_z"]].to_numpy()
      theta, phi, psi = angles
      xd = []
      a11 = np.cos(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.sin(psi)
      a12 = -np.sin(psi)*np.cos(psi) - np.cos(theta)*np.sin(phi)*np.cos(psi)
      a13 = np.sin(theta)*np.sin(phi)
      a21 = np.cos(psi) *np.sin(phi) + np.cos(theta)*np.cos(phi)*np.sin(psi)
      a22 = -np.sin(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.cos(psi)
      a23 = -np.sin(theta)*np.cos(phi)
      a31 = np.sin(theta)*np.sin(psi)
      a32 = np.sin(theta)*np.cos(psi)
      a33 = np.cos(theta)   
      
      R = np.array([
          [a11, a12, a13],
          [a21, a22, a23],
          [a31, a32, a33]
      ])
      
      for i in range(len(xyz)):
            xd.append(np.dot(np.dot(R, xyz[i]), R.T))
      return np.array(xd).ravel()

def rotate_NH_bond_vector_individual(NH_bond_vectors, angles): #Woessner
      xyz = NH_bond_vectors[["NH_x","NH_y","NH_z"]].to_numpy()
      theta, phi, psi = angles
      a11 = np.cos(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.sin(psi)
      a12 = -np.sin(psi)*np.cos(psi) - np.cos(theta)*np.sin(phi)*np.cos(psi)
      a13 = np.sin(theta)*np.sin(phi)
      a21 = np.cos(psi) *np.sin(phi) + np.cos(theta)*np.cos(phi)*np.sin(psi)
      a22 = -np.sin(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.cos(psi)
      a23 = -np.sin(theta)*np.cos(phi)
      a31 = np.sin(theta)*np.sin(psi)
      a32 = np.sin(theta)*np.cos(psi)
      a33 = np.cos(theta)   
      
      R = np.array([
          [a11, a12, a13],
          [a21, a22, a23],
          [a31, a32, a33]
      ])
      
      x, y, z = np.dot(np.dot(R, xyz), R.T)
      return x, y, z
  
def rotate_NH_bond_vector_individual_unique(NH_bond_vectors, angles): #Woessner
    xyz = np.array(NH_bond_vectors.to_numpy()).ravel()
    theta, phi, psi = angles
    a11 = np.cos(psi)*np.cos(psi) - np.cos(theta)*np.sin(phi)*np.sin(psi)
    a12 = -np.sin(psi)*np.cos(psi) - np.cos(theta)*np.sin(phi)*np.cos(psi)
    a13 = np.sin(theta)*np.sin(phi)
    a21 = np.cos(psi) *np.sin(phi) + np.cos(theta)*np.cos(phi)*np.sin(psi)
    a22 = -np.sin(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.cos(psi)
    a23 = -np.sin(theta)*np.cos(phi)
    a31 = np.sin(theta)*np.sin(psi)
    a32 = np.sin(theta)*np.cos(psi)
    a33 = np.cos(theta)   
    
    R = np.array([
        [a11, a12, a13],
        [a21, a22, a23],
        [a31, a32, a33]
    ])
    
    x, y, z = np.dot(np.dot(R, xyz), R.T)
    return x, y, z