import mod

# import molecules as smiles strings:
#First we define our input molecules
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
formaldehyde.print()

glycolaldehyde = Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")
glycolaldehyde.print()

input_molecules = [formaldehyde, glycolaldehyde]

#String containing a GML rule representing keto-enol isomerization
ketoEnolGML = """rule [
   ruleID "Keto-enol isomerization" 
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
   ]   
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
   ]   
   right [
      edge [ source 1 target 2 label "=" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]   
]"""

#String containing a GML rule representing aldol Addition 
aldolAddGML = """rule [
   ruleID "Aldol Addition"
   left [
      edge [ source 1 target 2 label "=" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 3 target 4 label "-" ]
      edge [ source 5 target 6 label "=" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
      node [ id 5 label "O" ]
      node [ id 6 label "C" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 5 target 6 label "-" ]
      edge [ source 4 target 5 label "-" ]
      edge [ source 6 target 1 label "-" ]
   ]
]"""

aldolAdd_F = Rule.fromGMLString(aldolAddGML)
aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)
inputRules = [aldolAdd_B, aldolAdd_F, ketoEnol_B, ketoEnol_F]

for rule in inputRules:
    rule.print()

#We can make even more complicated strategies
dg = mod.DG(graphDatabase=input_molecules)
reaction_network = dg.build()

#We define a strategy as follows
strategy = (
   addSubset(input_molecules) 
   >> repeat[4](inputRules)
   )

reaction_network.execute(strategy)
del reaction_network
dg.print()