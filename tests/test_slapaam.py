from slapmapper.aam import SlapAAM

def test_smiles():

    mapper = SlapAAM(binary=True)

    # Reaction SMILES
    rxn = "COC(=O)c1ccc(C=C2CC2)cc1>>COC(=O)c1ccc(C2=CCC2)cc1"
    mapper.reset()
    mapper.map_smiles(rxn)
    mapped_rxns = [ r["smiles"] for r in mapper.results ]

    assert len(mapped_rxns) == 1
    assert mapped_rxns[0] == '[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6][cH:7][c:8]([CH:9]=[C:10]2[CH2:11][CH2:12]2)[cH:13][cH:14]1>>[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6][cH:7][c:8]([C:10]2=[CH:9][CH2:11][CH2:12]2)[cH:13][cH:14]1'

def test_smiles_with_constraints():

    mapper = SlapAAM(binary=True)
    # Reaction SMILES with partial constraints
    rxn = "COC(=O)c1ccc(C=[C:1]2CC2)cc1>>COC(=O)c1ccc(C2=CC[C:1]2)cc1"
    mapper.reset()
    mapper.map_smiles(rxn)
    mapped_rxns = [ r["smiles"] for r in mapper.results ]

    assert len(mapped_rxns) == 1
    assert mapped_rxns[0] == '[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6][cH:7][c:8]([CH:9]=[C:10]2[CH2:11][CH2:12]2)[cH:13][cH:14]1>>[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6][cH:7][c:8]([C:9]2=[CH:11][CH2:12][C:10]2)[cH:13][cH:14]1'

def test_3d():

    mapper = SlapAAM(binary=True)
    # 3D-structure
    react = "tests/react.xyz"
    prod = "tests/prod.xyz"
    mapper.reset()
    mapper.map_3d(react,prod)
    mapped_rxns = [ r["mapping"] for r in mapper.results ]

    assert len(mapped_rxns) == 1
    assert mapped_rxns[0] == '0>>4;1>>5;2>>0;3>>3;4>>1;5>>2;6-8>>6-8'
