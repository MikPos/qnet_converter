######################################
# Converter Using Connectivity Files #
######################################
# # # # # # # 
# This converter works by using *.xyz.connectivity files to convert q_net reactions to MØD graphs
# It shoudl still use the *.xyz files since they contain the needed energy information.
# Graph structure should come form connectivity only though
#
#
#
import os
import json
import re


######################
# Reaction Constants #
######################
TEMP = "25"
SOLU = "h2o"
FUNC = "sp_m06-2x_qz"




####################
# Helper Functions #
####################

'''
These functions help us split the connectivity files created by QNet. 
These files have a nice structure of all needed information:
- Atom Labels
- Atom IDs
- Bond Type
- Connected Atoms

Example:
1C: 2O(-) 3H(-) 4H(-) 5H(-)
2O: 1C(-) 6O(-)
3H: 1C(-)
4H: 1C(-)
5H: 1C(-)
6O: 2O(-)
7O: 8C(=)
'''

def split_connectivity_file(filename):
    graph_strings = ["", "", "", "", "", ""]
    graph_id_sets = []
    with open(filename, "r") as file:
        for line in file:
            ids = []
            nums = re.sub(r"[^0-9\s]", "", line)
            for id in nums.split(" "):
                # if c.isdigit(): # This only looks at the current digit. We need something that looks at everything except for the character.
                ids.append(id.strip())
            added_to_sets = []
            if not ids:
                continue
            for s in graph_id_sets:
                if ids[0] in s: # If the current atom is in the set.
                    s.update(ids)
                    added_to_sets.append(graph_id_sets[graph_id_sets.index(s)])
            added_to_sets = trim_sets(graph_strings, graph_id_sets, added_to_sets)
            if not added_to_sets:
                new_set = set(ids)
                graph_id_sets.append(new_set)
                added_to_sets.append(new_set)
            index_of_set = graph_id_sets.index(added_to_sets[0])
            graph_strings[index_of_set] += line
    result_strings = []
    for s in graph_strings:
        filtered_string = s.strip()
        if filtered_string:
            result_strings.append(filtered_string)
    return result_strings


def trim_sets(graph_strings, graph_id_sets, added_to_sets):
    if len(added_to_sets) > 1:
        updated_set = set()
        joined_graph = ""
        sets_to_be_removed = []
        graphs_to_be_removed = []
        for added_set in added_to_sets:
            index = added_to_sets.index(added_set)
            sets_to_be_removed.append(graph_id_sets[index])
            graphs_to_be_removed.append(graph_strings[index])
            updated_set.update(graph_id_sets[index])
            joined_graph += graph_strings[index] + "\n"
        for s in sets_to_be_removed:
            graph_id_sets.remove(s)
        for s in graphs_to_be_removed:
            graph_strings.remove(s)
        graph_id_sets.append(updated_set)
        added_to_sets = [updated_set]
        index = graph_id_sets.index(updated_set)
        graph_strings[index] = joined_graph
    return added_to_sets
            
            
#################################################################################################################################
# The folowing functions help convert the connectivity files into the nodes and edges needed for building mod graphs and rules. #
#################################################################################################################################


def split_atom_label(atom):
    id = ""
    label = ""
    for c in atom:
        if c.isdigit():
            id += c
        else:
            label += c
    return id, label

def find_edge_target_and_label(edge_string):
    label = ""
    target = ""
    parts = edge_string.split("(")
    for c in parts[0]:
        if c.isdigit():
            target += c
    label = parts[1][0] # The first character after the opening parenthesis is always the edge label.
    return label, target

def create_nodes_and_edges(graph_string):
    graph_gml = ""
    for line in graph_string.split("\n"):
        if line:
            line_parts = line.split(":")
            node_id, node_label = split_atom_label(line_parts[0])
            graph_gml += f"""node [ id {node_id} label "{node_label}" ]\n"""
            for edge in line_parts[1].split(" "):
                if edge:
                    edge_label, target_id = find_edge_target_and_label(edge)
                    if node_id > target_id:
                        graph_gml += f"""edge [ source {node_id} target {target_id} label "{edge_label}" ]\n"""
    return graph_gml


##################
# GML Generators #
##################

def generate_graph_gml(graph_string):
    gml = """graph [\n"""
    gml += graph_string
    gml += """]"""
    return gml

def generate_rule_left_gml(graph_strings):
    gml = """left [\n"""
    for g_string in graph_strings:
        gml += g_string
    gml += """]"""
    return gml

def generate_rule_right_gml(graph_strings):
    gml = """right [\n"""
    for g_string in graph_strings:
        gml += g_string
    gml += """]"""
    return gml

def generate_rule_gml(list_of_left_graph_string, list_of_right_graph_string):
    rule_gml = """rule [\n"""
    rule_gml += generate_rule_left_gml(list_of_left_graph_string)
    rule_gml += generate_rule_right_gml(list_of_right_graph_string)
    rule_gml += """]"""
    # print(rule_gml)
    return rule_gml


######################
# Full Functionality #
######################
def make_mod_representation():
    reactions = {}
    graph_pair_db = {}
    rule_db = {}
    # These dictionaries should hold both the GML as well as the energy for each educt and product.
    # Energy comes from relation.json.
    # GMl comes from connectivity.
    # Filename (graph name) is the key.

    print("To run this code properly, please make sure that all .xyz files related of educts and products, as well as all connectivity files are available within the same folder.\n")
    print("If you have not curated the files you need, then please do so.\n")
    print("If you do not have the files available yet, please run q_net, q_net analyze and q_spider.\n\n")

    
    place = "files_for_mod"
    os.chdir(place)
    print(os.listdir())
    for filename in os.listdir():
        if os.path.isfile(filename):
            if filename.endswith(".json"):
                reactions = find_energies(filename)
                print("done with relations")
            # These filenmaes should be changed to fit the new connectivity files.
            elif filename.endswith(".xyz.connectivity"):
                read_qnet_graphs(graph_pair_db, rule_db, filename, ".xyz.connectivity")
                print(f"Done with the file: {filename}")
    print("Done, importing files!")
    os.chdir(os.pardir)
    resulting_dg = build_mod_dg(reactions, graph_pair_db, rule_db)
    p = DGPrinter()
    p.withRuleName = True
    graph_p = p.graphPrinter
    # graph_p.collapseHydrogens = False
    # graph_p.simpleCarbons = False
    # p.withInlineGraphs = True
    resulting_dg.print(p)


# Constructs a dictionary of the graphs in the different connectivity files. The graphs are saved as a list and the key is the filename (without extension).
def read_qnet_graphs(graph_dict, rule_dict, filename, file_ending):
    graphs = split_connectivity_file(filename)
    gmlStrings = []
    rule_string = ""
    for i in range(len(graphs)):
        if graphs[i]:
            gmlStrings.append(create_nodes_and_edges(graphs[i]))
    for graph_string in gmlStrings:
        rule_string += graph_string
    ending_len = len(file_ending) * -1
    key = filename[0:ending_len]
    graph_dict[key] = gmlStrings
    rule_dict[key] = rule_string


# Check to make sure: Reaction is "ch3oo-vitc/PATH_0_1"
# Dict like this: {"ch3oo-vitc/PATH_0_1" : [100.37, 98.91, ["ch3oo", "vitc"], ["mol_0"]]}
def find_energies(json_file):
    result_dict = {}
    with open(json_file, 'r') as file:
        reaction_json = json.load(file)
        for reaction, data in reaction_json.items():
            educts = []
            products = []
            ga_solv = None
            gr_solv = None
            for option, value in data.items():
                if option == "educts":
                    for mol in value:
                        educts.append(mol)
                if option == "products":
                    for mol in value:
                        products.append(mol)
            for type, energy in data["G"][TEMP][SOLU][FUNC].items():
                if type == "G_R_solv":
                    gr_solv = energy
                if type == "G_A_solv":
                    ga_solv = energy
            result_dict[reaction] = [ga_solv, gr_solv, educts, products]
    return result_dict


def build_mod_dg(reaction_dict, connectivity_dict, connectivity_rule_dict):
    graph_db = {}
    gml_db = {}
    rule_db = {}
    # Iterate through keys of energy dictionary.
    # We look for the educts and products in the graph DB.
    # If they don't exist we read the dictionary of the connectivity files and create them.
    # Then we make a derivation from the educts to the products of the reaction.
    # # The hyperedge should be labelled by the path function, and the edges should also be labelled with the energies.
    dg = mod.DG(graphDatabase=[])
    b = dg.build()

    for key in reaction_dict.keys():
        e_a = reaction_dict[key][0]
        e_r = reaction_dict[key][1]
        educts_names = reaction_dict[key][2]
        product_names = reaction_dict[key][3]

        educt_product = key.split("/")
        educts = educt_product[0]
        products = educt_product[1]
        print(educt_product)

        # I need to import the remaining molecules...
        if ("mol" in educts):
            break

        educt_gmls, educts_graphs = collect_graphs(connectivity_dict, graph_db, gml_db, educts, educts_names)
        product_gmls, products_graphs = collect_graphs(connectivity_dict, graph_db, gml_db, products, product_names)

        rule_left_gml = connectivity_rule_dict[educts]
        rule_right_gml = connectivity_rule_dict[products]

        new_rule_string = generate_rule_gml(rule_left_gml, rule_right_gml)
        new_rule_name = f"{products.replace('_', '-')}, E_A: {round(e_a, 2)}, E_R: {round(e_r,2)}"
        new_rule = mod.ruleGMLString(new_rule_string, name=new_rule_name)

        with open(f"{new_rule_name}.txt", "w") as f:
            f.writelines(new_rule_string.split("\n"))

        print(f"Generated Rule for: {key}")

        d = mod.Derivation()
        d.left = educts_graphs
        d.right = products_graphs
        d.rule = new_rule
        b.addDerivation(d)
        new_rule.print()
        print(f"Did Derivation for: {key}")
        print("-----------------------------------")

        rule_db[key] = new_rule

    return dg

def collect_graphs(connectivity_dict, graph_db, gml_db, molecules, molecule_names):
    result = []
    gmls = []
    mols = connectivity_dict[molecules]
    for i in range(len(mols)):
        mol_name = molecule_names[i]
        print("Graph name: " + mol_name)
        print(mols[i])
        if mol_name not in graph_db.keys():
            graph_string = mols[i]
            gml_string = generate_graph_gml(graph_string)
            gmls.append(graph_string) # This is the the gml without the surrounding graph [] tag.
            mod_graph = mod.graphGMLString(gml_string, name=mol_name)
            result.append(mod_graph)
            graph_db[mol_name] = mod_graph
            gml_db[mol_name] = graph_string
            # Create graphs
            with open(f"{mol_name}.txt", "w") as f:
                f.writelines(gml_string.split("\n"))
        else:
            result.append(graph_db[mol_name])
            gmls.append(gml_db[mol_name])
    return gmls, result


make_mod_representation()