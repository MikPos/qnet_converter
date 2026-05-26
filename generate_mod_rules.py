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


def read_SDF_withBO(file):

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

def read_SDF_woBO(file):

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

def convert_bo(number):

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
    else:
        full_context = False

    for file in glob.glob("mod_rules*.py"):
        os.remove(file)

    with open("relations.json", "r") as json_file:
        relations = json.load(json_file)
    with open("molecules.json", "r") as json_file:
        mol_json = json.load(json_file)

    if len(sys.argv) == 3:
        keys = [sys.argv[-2]]
    else:
        keys = relations.keys()

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
        if not os.path.isfile(key.replace("/PATH","_PATH")+"/unpaired.txt"):
            print(key)
            print("No "+key.replace("/PATH","_PATH")+"/unpaired.txt")
            tot_unpaired = 0
        else:
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

        ed_conn = read_SDF_woBO(key+"/educt_new.sdf")
        prod_conn = read_SDF_woBO(key+"/product_new.sdf")

        ed_conn_BO, ed_radicals, ed_charge = read_SDF_withBO(key+"/educt_new.sdf") 
        prod_conn_BO, prod_radicals, prod_charge = read_SDF_withBO(key+"/product_new.sdf")

        print("ed_radicals, prod_radicals:", ed_radicals, prod_radicals)
        #breaks = compare_connectivity(prod_conn, ed_conn)
        #adds = compare_connectivity(ed_conn, prod_conn)
        #print("breaks", breaks)
        #print("adds", adds)

        bo_changes1 = []
        bo_changes2 = []

        atom_labels = read_SDF_atoms(key+"/educt_new.sdf")

        educt_bonds = compare_connectivity(prod_conn_BO, ed_conn_BO)
        #for item in bo_changes:
        #    if ".".join(item.split(".")[:-1]) in breaks + adds:
        #        continue
        #    else:
        #        bo_changes1.append(item)

        product_bonds = compare_connectivity(ed_conn_BO, prod_conn_BO)
        #for item in bo_changes:
        #    if ".".join(item.split(".")[:-1]) in breaks + adds:
        #        continue
        #    else:
        #        bo_changes2.append(item)

        nodes = []
        for items in educt_bonds+product_bonds:
            for it in items.split(".")[:2]:
                if it not in nodes:
                    nodes.append(it)

        all_nodes = []
        for item in ed_conn_BO:
            for it in item[:2]:
                if str(it) not in all_nodes:
                    all_nodes.append(str(it))

        ed_conn = compare_connectivity([],ed_conn_BO)
        ed_conn_final = []
        for item in ed_conn:
            if item in educt_bonds or item in product_bonds:
                continue
            else:
                ed_conn_final.append(item)

        #charge treatment missing!!!

        left_string = ""
        for_right_string = ""
        for item in educt_bonds:
            left_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bo(item.split(".")[2]))+"\" ]\n"
        if ed_radicals != []:
            for_right_string =""
            for it in ed_radicals:
                if prod_radicals == [] or it not in prod_radicals:
                    left_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+".\" ]\n"
                    for_right_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+"\" ]\n"
                    all_nodes.remove(it)

        right_string = ""
        for_left_string = ""
        for item in product_bonds:
            right_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bo(item.split(".")[2]))+"\" ]\n"
        if prod_radicals != []:
            for it in prod_radicals:
                if ed_radicals == [] or it not in ed_radicals:
                    right_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+".\" ]\n"
                    for_left_string += "\t\t\tnode [ id "+str(it)+" label \""+str(atom_labels[it])+"\" ]\n"
                    all_nodes.remove(it)
  
        right_string += for_right_string
        left_string += for_left_string

        context_string = ""
        if full_context:
            for node in all_nodes:
                if node in ed_radicals and node in prod_radicals:
                    context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+".\" ]\n"
                else:
                    context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+"\" ]\n"
            for item in ed_conn_final:
                context_string += "\t\t\t edge [ source "+str(item.split(".")[0])+" target "+str(item.split(".")[1])+" label \""+str(convert_bo(item.split(".")[2]))+"\" ]\n"

        else:
            for node in nodes:
                context_string += "\t\t\tnode [ id "+str(node)+" label \""+str(atom_labels[node])+"\" ]\n" 

        rule_name = key.replace("/PATH","_PATH").replace("-", "_")

        mod_rules_file = "mod_rules.py"

        with open(mod_rules_file, "a") as out:
            out.write(str(rule_name)+""" = \"\"\"rule [
\truleID \" """+str(rule_name)+"""\"
\tleft [\n"""
+left_string+"""
\t]
\tcontext [\n"""
+context_string+"""
\t]
\tright [\n"""
+right_string+"""
\t]
]\"\"\"
rule1_F = Rule.fromGMLString("""+str(rule_name)+""")
""")
#rule1_B = Rule.fromGMLString("""+str(rule_name)+""", invert=True)
#""")
