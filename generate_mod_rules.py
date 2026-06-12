#!/usr/bin/env python

import os, sys
import json
from mgsm.qnet_molecule import qnet_molecule
from rdmc import RDKitMol
from rdmc.fix import fix_mol
from rdkit.Chem.rdmolfiles import MolFromXYZFile, MolToSmiles
import os, sys
import glob



def convert_xyz(xyz_file, charge=0, unpaired=0):

    print(xyz_file, charge, unpaired, type(charge), type(unpaired))
    if xyz_file[-4:] == ".xyz":
       with open(xyz_file, "r") as f:
           xyz = f.read()
    try:
      try:
        mol = RDKitMol.FromXYZ(xyz, sanitize=False, backend="jensen", charge=int(charge), allow_charged_fragments=False, force_rdmc=False)
      except:
        mol = RDKitMol.FromXYZ(xyz, sanitize=False, backend="jensen", charge=int(charge), allow_charged_fragments=False, force_rdmc=True)
      mol = fix_mol(mol, fix_spin_multiplicity=True, mult=int(unpaired)+1, sanitize=True)
      smi = mol.ToSmiles(removeHs=False, removeAtomMap=True)
    except:
      try:
        try:
          mol = RDKitMol.FromXYZ(xyz, sanitize=False, backend="jensen", charge=int(charge), allow_charged_fragments=True, force_rdmc=False)
        except:
          mol = RDKitMol.FromXYZ(xyz, sanitize=False, backend="jensen", charge=int(charge), allow_charged_fragments=True, force_rdmc=True)
        mol = fix_mol(mol, fix_spin_multiplicity=True, mult=int(unpaired)+1, sanitize=True)
        smi = mol.ToSmiles(removeHs=False, removeAtomMap=True)
      except:
        try:
          smi = mol.ToSmiles(removeHs=False, removeAtomMap=True)
          mol = RDKitMol.FromXYZ(xyz, sanitize=False)
          mol = fix_mol(mol, fix_spin_multiplicity=True, mult=int(unpaired)+1, sanitize=True)
          smi = mol.ToSmiles(removeHs=False, removeAtomMap=True)
        except:
          print("Problem with ToSmiles conversion:", xyz_file, charge, unpaired)
          sys.exit()

    sdf = mol.ToSDFFile(xyz_file.replace(".xyz", "_new.sdf"))

    return smi


def read_SDF_with_bond_order(file):

    """ Read in the connectivity from an SDF file and return the connectivity as a list """

    conn = []
    inp = open(file, "r")
    inpdata = inp.readlines()
    start = int(inpdata[3].split()[0]) + 4
    end = start + int(inpdata[3].split()[1])
    for i in range(start, end):
        val1 = int(inpdata[i].split()[0])
        val2 = int(inpdata[i].split()[1])
        bo = int(inpdata[i].split()[2])
        if val1 > val2 and (val2, val1) not in conn:
            conn.append((val2, val1, bo))
        elif val2 > val1 and (val1, val2) not in conn:
            conn.append((val1, val2, bo))

    radicals = []
    charges = []

    for line in inpdata:
        if "M  RAD" in line:
            nrads = line.split()[2]
            if nrads == "1":
                radicals = [line.split()[3]]
            elif nrads == "2":
                radicals = [line.split()[3], line.split()[5]]
        if "M  CHG" in line:
            ncharges = line.split()[2]
            if ncharges == "1":
                charges = [[line.split()[3], line.split()[4]]]
            elif ncharges == "2":
                charges = [[line.split()[3], line.split()[4]], [line.split()[5], line.split()[6]]]

    return conn, radicals, charges

def read_SDF_without_bond_order(file):

    """ Read in the connectivity from an SDF file and return the connectivity as a list """

    conn = []
    inp = open(file, "r")
    inpdata = inp.readlines()
    start = int(inpdata[3].split()[0]) + 4
    end = start + int(inpdata[3].split()[1])
    for i in range(start, end):
        val1 = int(inpdata[i].split()[0])
        val2 = int(inpdata[i].split()[1])
        if val1 > val2 and (val2, val1) not in conn:
            conn.append((val2, val1))
        elif val2 > val1 and (val1, val2) not in conn:
            conn.append((val1, val2))

    return conn

def read_SDF_atoms(file):
    """ Read in the atom XYZ info from an SDF file and return  """

    inp = open(file, "r")
    inpdata = inp.readlines()
    start = 4 
    end = int(inpdata[3].split()[0]) + 4

    num = 1
    label = {}

    for i in range(start, end):
        label[str(num)] = inpdata[i].split()[3]
        num += 1
         
    return label

def compare_connectivity(bondlist1, bondlist2):

    """ Compare two connectivity lists and return the difference: bondlist2-bondlist1 """

    diff1 = set(bondlist2) - set(bondlist1)

    diff = []
    for item1 in diff1:
        if len(item1)==3:
            diff.append(str(item1[0])+"."+str(item1[1])+"."+str(item1[2]))
        else:
            diff.append(str(item1[0])+"."+str(item1[1]))


    return sorted(diff)

def convert_bond_order(number):

    if int(number) == 1:
        return "-"
    if int(number) == 2:
        return "="
    if int(number) == 3:
        return "#"
    else:
        return None

if  __name__=="__main__":

    if sys.argv[-1] == "full":
        full_context = True
        for file in glob.glob("full_context_rules/qnet_rule*.gml"):
            os.remove(file)
    else:
        full_context = False
        for file in glob.glob("no_context_rules/qnet_rule*.gml"):
            os.remove(file)

    with open("relations.json", "r") as json_file:
        relations = json.load(json_file)
    with open("molecules.json", "r") as json_file:
        mol_json = json.load(json_file)

    if len(sys.argv) == 3:
        keys = [sys.argv[-2]]
    else:
        keys = relations.keys()
    
    mols = mol_json.keys()

    for mol in mols:
        print(mol)
        mol_xyz =  mol + "/" + mol + ".xyz"

        with open("molecules/"+mol+"/charge.txt", "r") as inp:
            tot_charge = int(inp.readlines()[0])
            print(tot_charge)
        with open("molecules/"+mol+"/unpaired.txt", "r") as inp:
            tot_unpaired = int(inp.readlines()[0].strip("u"))
            print(tot_unpaired)

        convert_xyz(mol_xyz, charge=tot_charge, unpaired=tot_unpaired)

        mol_conn = read_SDF_without_bond_order(mol+"/"+mol+"_new.sdf")

        mol_conn_BO, mol_radicals, mol_charges = read_SDF_with_bond_order(mol+"/"+mol+"_new.sdf")
        mol_charge_atoms = []
        for atom in mol_charges:
            mol_charge_atoms.append(atom[0])

        atom_labels = read_SDF_atoms(mol+"/"+mol+"_new.sdf")

        mol_gml = "graph [\n"

        mol_nodes = []
        for item in mol_conn_BO:
            mol_gml += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bond_order(item.split(".")[2]))+"\" ]\n"
            for it in item.split(".")[:2]:
                if it not in mol_nodes:
                    if it in mol_radicals:
                        mol_gml += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+".\" ]\n"
                    if it in mol_charge_atoms:
                        if mol_charges[it][1] == "1":
                            sign = "+"
                        elif mol_charges[it][1] == "2":
                            sign = "+2"
                        elif mol_charges[it][1] == "-1":
                            sign = "-"
                        elif mol_charges[it][1] == "-2":
                            sign = "-2"
                        mol_gml += "\t\t\tnode [ id "+str(it[0])+" label \""+str(atom_labels[it[0]])+sign+"\" ]\n"
                    mol_nodes.append(it)
        mol_gml += "]"
        
        mod_file_name = f"mod_molecules/{mol}.gml"

        with open(mod_file_name, "w") as out:
            out.write(mol_gml)




    for key in keys:
        print(key)
        educt_xyz = key+"/educt.xyz"
        product_xyz = key+"/product.xyz"
        # read charge and unpaired of complex/TS!
        if not os.path.isdir("ts/"+key.replace("/PATH","_PATH")):
            continue
        if os.path.isfile("ts/"+key.replace("/PATH","_PATH")+"/STOP"):
            continue
        with open("ts/"+key.replace("/PATH","_PATH")+"/charge.txt", "r") as inp:
            tot_charge = int(inp.readlines()[0])
            print(tot_charge)
        with open("ts/"+key.replace("/PATH","_PATH")+"/unpaired.txt", "r") as inp:
            tot_unpaired = int(inp.readlines()[0].strip("u"))
            print(tot_unpaired)
        educts = relations[key]["educts"]
        #for ed in educts:
        #    ed_smi = convert_xyz("molecules/"+ed+"/"+ed+".xyz", charge=mol_json[ed]["charge"], unpaired=mol_json[ed]["charge"])
        #    ed_conn, ed_radicals, ed_charges = read_SDF_withBO("molecules/"+ed+"/"+ed+"_new.sdf")
        #    ed_conn_final = compare_connectivity([],ed_conn)
        #    nodes = []
        #    for items in ed_conn_final:
        #        for it in items.split(".")[:2]:
        #            if it not in nodes:
        #                nodes.append(it)
        #    atom_labels = read_SDF_atoms("molecules/"+ed+"/"+ed+"_new.sdf")
        #    context_string = ""
        #    for node in nodes:
        #        context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+"\" ]\n"
        #    for item in ed_conn_final:
        #        context_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bo(item.split(".")[2]))+"\" ]\n"
        #    print(context_string)
        #sys.exit()
        if not os.path.isfile(educt_xyz) or not os.path.isfile(product_xyz):
            print("educt/product.xyz missing for: ", key)
            continue

        convert_xyz(educt_xyz, charge=tot_charge, unpaired=tot_unpaired)
        convert_xyz(product_xyz, charge=tot_charge, unpaired=tot_unpaired)

        ed_conn = read_SDF_without_bond_order(key+"/educt_new.sdf")
        prod_conn = read_SDF_without_bond_order(key+"/product_new.sdf")

        ed_conn_BO, ed_radicals, ed_charges = read_SDF_with_bond_order(key+"/educt_new.sdf") 
        prod_conn_BO, prod_radicals, prod_charges = read_SDF_with_bond_order(key+"/product_new.sdf")

        print("ed_radicals, prod_radicals:", ed_radicals, prod_radicals)

        bo_changes1 = []
        bo_changes2 = []

        atom_labels = read_SDF_atoms(key+"/educt_new.sdf")

        educt_bonds = compare_connectivity(prod_conn_BO, ed_conn_BO)
        product_bonds = compare_connectivity(ed_conn_BO, prod_conn_BO)

        # educt_string = ""
        # product_string = ""

        all_nodes = []
        for item in ed_conn_BO:
            for it in item[:2]:
                if str(it) not in all_nodes:
                    all_nodes.append(str(it))

        nodes = []
        for items in educt_bonds+product_bonds:
            for it in items.split(".")[:2]:
                if it not in nodes:
                    nodes.append(it)

        ed_conn = compare_connectivity([],ed_conn_BO)
        ed_conn_final = []
        for item in ed_conn:
            if item in educt_bonds or item in product_bonds:
                continue
            else:
                ed_conn_final.append(item)

        prod_charge_ids = []
        for prod_charge in prod_charges:
            prod_charge_ids.append(prod_charge[0])

        left_string = ""
        for_right_string = ""
        for item in educt_bonds:
            left_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bond_order(item.split(".")[2]))+"\" ]\n"
            # educt_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bond_order(item.split(".")[2]))+"\" ]\n"
        if ed_radicals != []:
            for_right_string =""
            for it in ed_radicals:
                if (prod_radicals == [] or it not in prod_radicals) and it not in prod_charge_ids:
                    left_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+".\" ]\n"
                    for_right_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+"\" ]\n"
                    if it in all_nodes:
                        all_nodes.remove(it)
        if ed_charges != []:
            for_right_string =""
            for it in ed_charges:
                if (prod_charges == [] or it not in prod_charges) and it[0] not in prod_radicals:
                    if it[1] == "1":
                        sign = "+"
                    elif it[1] == "2":
                        sign = "+2"
                    elif it[1] == "-1":
                        sign = "-"
                    elif it[1] == "-2":
                        sign = "-2"
                    left_string += "\t\t\tnode [ id "+str(it[0])+" label \""+str(atom_labels[it[0]])+sign+"\" ]\n"
                    for_right_string += "\t\t\tnode [ id "+str(it[0])+" label \""+str(atom_labels[it[0]])+"\" ]\n"
                    if it[0] in all_nodes:
                        all_nodes.remove(it[0])

        right_string = ""
        for_left_string = ""

        ed_charge_ids = []
        for ed_charge in ed_charges:
            ed_charge_ids.append(ed_charge[0])

        for item in product_bonds:
            right_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bond_order(item.split(".")[2]))+"\" ]\n"
        if prod_radicals != []:
            for it in prod_radicals:
                if ed_radicals == [] or it not in ed_radicals:
                    right_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+".\" ]\n"
                    for_left_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+"\" ]\n"
                    if it in all_nodes:
                        all_nodes.remove(it)
        if prod_charges != []:
          for it in prod_charges:
              if (ed_charges == [] or it not in ed_charges):# and it[0] not in ed_radicals:
                  if it[1] == "1":
                      sign = "+"
                  elif it[1] == "2":
                      sign = "+2"
                  elif it[1] == "-1":
                      sign = "-"
                  elif it[1] == "-2":
                      sign = "-2"
                  right_string += "\t\t\tnode [ id "+str(it[0])+" label \""+str(atom_labels[it[0]])+sign+"\" ]\n"
                  for_left_string += "\t\t\tnode [ id "+str(it[0])+" label \""+str(atom_labels[it[0]])+"\" ]\n"
                  if it[0] in all_nodes:
                    all_nodes.remove(it[0])

        right_string += for_right_string
        left_string += for_left_string

        # Create context for the full context rules
        full_context_string = ""
        for node in all_nodes:
            if node in ed_radicals and node in prod_radicals:
                full_context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+".\" ]\n"
            else:
                full_context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+"\" ]\n"
        for item in ed_conn_final:
            full_context_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bond_order(item.split(".")[2]))+"\" ]\n"

        # Create context for the no context rules (nodes that take part in reaction)
        no_context_string = ""
        for node in nodes:
            no_context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+"\" ]\n" 

        rule_name = key.replace("/PATH","_PATH").replace("-", "_")
        file_full_context = f"full_context_rules/qnet_rule_{rule_name}.gml"
        file_no_context = f"no_context_rules/qnet_rule_{rule_name}.gml"

        with open(file_full_context, "a") as out:
            out.write("""rule [
\truleID \" """+str(rule_name)+"""\"
\tleft [\n"""
+left_string+"""
\t]
\tcontext [\n"""
+full_context_string+"""
\t]
\tright [\n"""
+right_string+"""
\t]
]""")

        with open(file_no_context, "a") as out:
            out.write("""rule [
\truleID \" """+str(rule_name)+"""\"
\tleft [\n"""
+left_string+"""
\t]
\tcontext [\n"""
+no_context_string+"""
\t]
\tright [\n"""
+right_string+"""
\t]
]""")
#rule1_B = Rule.fromGMLString("""+str(rule_name)+""", invert=True)
#""")
