import rdkit.Chem as Chem
from xyz2mole import xyz2mole
import mod
import json

bond_label = {
    Chem.rdchem.BondType.SINGLE: "-",
    Chem.rdchem.BondType.DOUBLE: "=",
    Chem.rdchem.BondType.TRIPLE: "#",
    Chem.rdchem.BondType.HYDROGEN: "-",
}


def convert_xyz_to_mol(pathString):
    # Read the .xyz file
    raw_mol = Chem.MolFromXYZFile(pathString)
    # Convert to RDKit mol object
    mol = xyz2mol(raw_mol)
    # Print the SMILES representation
    print(Chem.MolToSmiles(mol))
    return mol


# convert to MØD graph:
def convert_mol_to_gml(mol):
    gml = """graph [\n"""
    for vertex in mol.getAtoms():
        gml += f"""node [ id {vertex.GetIdx()} label {vertex.GetSymbol()}]\n"""
        for edge in vertex.GetBonds():
            source = edge.GetBeginAtomIdx()
            target = edge.GetEndAtomIdx()
            if source > target:
                gml += f"""edge [ source {source} target {target} label {bond_label[edge.GetBondType()]}]\n"""
    gml += """]"""
    return gml


def convert_rule_to_gml(fileName):
    with open(fileName, 'r') as file:
        data = json.load(file)



