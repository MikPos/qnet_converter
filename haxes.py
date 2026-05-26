import mod

from typing import Union

def rule_from_reaction_smiles(smiles_string: str, name=None, allowAbstract=False, *, invert: bool=False, add: bool=True) -> mod.Rule:
	"""
	Creates a mod rule from a reaction smarts string.
	Starts by splitting the smarts string into its left and right side, and each of those into their respective components.
	Using the external ids -> create the gml lines for a gml rule.
	Return the final GML rule.
	"""
	smiles_left, smiles_right = smiles_string.split(">>")
	smiles_left_components = smiles_left.split(".")
	smiles_right_components = smiles_right.split(".")
	mod_left_graph_strings = [mod.Graph.fromSMILES(s, allowAbstract=allowAbstract, add=False) for s in smiles_left_components]
	mod_right_graph_strings = [mod.Graph.fromSMILES(s, allowAbstract=allowAbstract, add=False) for s in smiles_right_components]
	def create_mod_string(graph):
		ext_from_int = {}
		for iExt in range(graph.minExternalId, graph.maxExternalId + 1):
			v = graph.getVertexFromExternalId(iExt)
			if not v.isNull():
				ext_from_int[v] = iExt
		mod_string = ""
		for v in graph.vertices:
			assert v in ext_from_int
			mod_string += '\t\tnode [ id %d label "%s" ]\n' % (ext_from_int[v], v.stringLabel)
		for e in graph.edges:
			mod_string += '\t\tedge [ source %d target %d label "%s" ]\n' % (ext_from_int[e.source], ext_from_int[e.target], e.stringLabel)
		return mod_string
	mod_rule_string = "rule [\n\tleft [\n"
	for mod_graph in mod_left_graph_strings:
		mod_rule_string += create_mod_string(mod_graph)
	mod_rule_string += "\t]\n\tright [\n"
	for mod_graph in mod_right_graph_strings:
		mod_rule_string += create_mod_string(mod_graph)
	mod_rule_string += "\t]\n]\n"
	return mod.Rule.fromGMLString(mod_rule_string, name=name, invert=invert, add=add)


def reaction_smiles_from_rule(r: mod.Rule) -> str:
	"""
	Uses two subfunctions in order to convert a mod rule to a reaction smiles rule.
	"""
	def replace_ids(smiles_string: str, mod_graph: mod.Graph) -> str:
		"""
		Converts the ids of vertices in a smiles string, to unique ones.
		"""
		import re
		# put an 'o' in front of each ID to distinguish it as original
		### That is not what the line below does...
		smiles_string = re.sub(":([0-9]+)]", ":o\\1]", smiles_string)
		# now replace the original IDs that are the IDs from loaded graphs
		# replace with the original IDs in the rule
		for i in range(mod_graph.minExternalId, mod_graph.maxExternalId + 1):
			v = mod_graph.getVertexFromExternalId(i)
			if not v:
				continue
			smiles_string = smiles_string.replace(f":o{v.id}]", f":{i}]")
		return smiles_string

	def make_side_graphs(side_graphs: Union[mod.Rule.LeftGraph, mod.Rule.RightGraph]) -> str:
		"""
		Creates the reaction smarts string from a single side of a mød rule.
		"""
		smiles_string = "graph [\n"
		for v in side_graphs.vertices:
			smiles_string += f'	node [ id {v.id} label "{v.stringLabel}" ]\n'
		for e in side_graphs.edges:
			smiles_string += f'	edge [ source {e.source.id} target {e.target.id} label "{e.stringLabel}" ]\n'
		smiles_string += "]\n"
		mod_graphs = mod.Graph.fromGMLStringMulti(smiles_string, add=False)
		for mod_graph in mod_graphs:
			smiles_string = mod_graph.smilesWithIds
			smiles_string = replace_ids(smiles_string, mod_graph)
		return ".".join(replace_ids(g.smilesWithIds, g) for g in mod_graphs)
	return f"{make_side_graphs(r.left)}>>{make_side_graphs(r.right)}"
