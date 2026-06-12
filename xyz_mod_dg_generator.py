import mod
import json

######################
# Reaction Constants #
######################
TEMP = "25"
SOLU = "h2o"
FUNC = "sp_m06-2x_qz"


######################
# Full Functionality #
######################
def make_mod_representation():
    reactions = {}
    graph_db = {}

    # print("To run this code properly, please make sure that all .xyz files related of educts and products, as well as all connectivity files are available within the same folder.\n")
    # print("If you have not curated the files you need, then please do so.\n")
    # print("If you do not have the files available yet, please run q_net, q_net analyze and q_spider.\n\n")
    with open("molecules.json", "r") as json_file:
        mol_json = json.load(json_file)
    mol_names = mol_json.keys()
    for mol in mol_names:
        graph = mod.Graph.FromGMLFile("mod_molecules/"+mol+".gml")
        graph_db[mol] = graph

    reactions = find_energies("relations.json")

    print("Done, importing files!")
    resulting_dg = build_mod_dg(reactions, graph_db)
    p = mod.DGPrinter()
    p.withRuleName = True
    graph_p = p.graphPrinter
    # graph_p.collapseHydrogens = False
    # graph_p.simpleCarbons = False
    # p.withInlineGraphs = True
    resulting_dg.print(p)


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


def build_mod_dg(reaction_dict, graph_dict):
    rule_db = {}
    # Iterate through keys of energy dictionary.
    # We look for the educts and products in the graph DB.
    # If they don't exist we read the dictionary of the connectivity files and create them.
    # Then we make a derivation from the educts to the products of the reaction.
    # # The hyperedge should be labelled by the path function, and the edges should also be labelled with the energies.
    dg = mod.DG(graphDatabase=[])
    b = dg.build()

    for key in reaction_dict.keys():
        print(key)
        e_a = reaction_dict[key][0]
        e_r = reaction_dict[key][1]
        educts_names = reaction_dict[key][2]
        products_names = reaction_dict[key][3]

        rule_name = key.split("/")[1]

        new_rule_name = f"{rule_name.replace('_', '-')}, E_A: {round(e_a, 2)}, E_R: {round(e_r,2)}"
        new_rule = mod.fromGMLFile(f"full_context_rules/qnet_rule{key.replace("/", "_")}.gml", name=new_rule_name)

        educts_graphs = []
        for ed in educts_names:
            educts_graphs.append(graph_dict[ed])
        products_graphs = []
        for prod in products_names:
            products_graphs.append(graph_dict[prod])

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

# Running code:
make_mod_representation()