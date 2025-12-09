from slapmapper.aam import SlapAAM

mapper = SlapAAM(binary=True)

# Reaction SMILES
rxn = "COC(=O)c1ccc(C=C2CC2)cc1>>COC(=O)c1ccc(C2=CCC2)cc1"
mapper.reset()
mapper.map_smiles(rxn)
mapped_rxns = [ r["smiles"] for r in mapper.results ]
print(mapped_rxns)

# Reaction SMILES with partial constraints
rxn = "COC(=O)c1ccc(C=[C:1]2CC2)cc1>>COC(=O)c1ccc(C2=CC[C:1]2)cc1"
mapper.reset()
mapper.map_smiles(rxn)
mapped_rxns = [ r["smiles"] for r in mapper.results ]
print(mapped_rxns)

# 3D-structure
react = "react.xyz"
prod = "prod.xyz"
mapper.reset()
mapper.map_3d(react,prod)
mapped_rxns = [ r["mapping"] for r in mapper.results ]
print(mapped_rxns)

