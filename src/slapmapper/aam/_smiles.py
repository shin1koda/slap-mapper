try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
from collections import Counter,defaultdict
from ..core import LabeledGraph


def count_elements_by_atomic_num(mol):

    element_count = Counter()
    for atom in mol.GetAtoms():
        element_count[atom.GetAtomicNum()] += 1
    return element_count


def balance_elements(mol1, mol2):

    diff = count_elements_by_atomic_num(mol2)
    diff.subtract(count_elements_by_atomic_num(mol1))

    if not diff:
        return Chem.Mol(mol1), Chem.Mol(mol2)

    rwmol1 = Chem.RWMol(mol1)
    rwmol2 = Chem.RWMol(mol2)

    for atomic_num, count in diff.items():
        if count>0:
            for _ in range(abs(count)):
                rwmol1.AddAtom(Chem.Atom(atomic_num))
        elif count<0:
            for _ in range(abs(count)):
                rwmol2.AddAtom(Chem.Atom(atomic_num))

    return rwmol1.GetMol(), rwmol2.GetMol()


def get_labeled_graph_from_mol(mol):

    labels = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    graph = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()

        others = {}
        for bond in atom.GetBonds():
            other_idx = bond.GetOtherAtomIdx(idx)
            bond_order = int(2*bond.GetBondTypeAsDouble()+0.1)
            others[other_idx] = bond_order
        graph[idx] = others

    return LabeledGraph(graph,labels)


def smiles2lgp(rxn_smiles,add_Hs=True):

    natoms_pair = [[],[]]
    sources_pair = [[],[]]

    ss = rxn_smiles.split('>')
    r = Chem.MolFromSmiles(ss[0])
    p = Chem.MolFromSmiles(ss[2])

    sources_pair[0].append(r)
    sources_pair[1].append(p)
    natoms_pair[0].append(r.GetNumAtoms())
    natoms_pair[1].append(p.GetNumAtoms())

    if add_Hs:
        r = Chem.AddHs(r)
        p = Chem.AddHs(p)

    sources_pair[0].append(r)
    sources_pair[1].append(p)
    natoms_pair[0].append(r.GetNumAtoms())
    natoms_pair[1].append(p.GetNumAtoms())

    r,p = balance_elements(r,p)

    sources_pair[0].append(r)
    sources_pair[1].append(p)
    natoms_pair[0].append(r.GetNumAtoms())
    natoms_pair[1].append(p.GetNumAtoms())


    atomic_nums_pair = [[atom.GetAtomicNum() for atom in mol.GetAtoms()] for mol in [r,p]]

    lgp = []
    for mol in [r,p]:
        lgp.append(get_labeled_graph_from_mol(mol))

    ini_l2i_pair = [defaultdict(list),defaultdict(list)]
    for mol,ini_l2i in zip([r,p],ini_l2i_pair):
        for atom in mol.GetAtoms():
            l = atom.GetAtomMapNum()
            if l>0:
                ini_l2i[l*1000].append(atom.GetIdx())

    for l,idxs0 in ini_l2i_pair[0].items():
        if l in ini_l2i_pair[1]:
            idxs1 = ini_l2i_pair[1][l]

            if len(idxs0)!=len(idxs1):
                print('Unbalanced initial labels detected')

            for lg,idxs in zip(lgp,[idxs0,idxs1]):
                for idx in idxs:
                    lg.labels[idx] = l
                lg.build_label2idxs()

    for lg, atomic_nums, natoms_slices, sources in zip(
        lgp, atomic_nums_pair, natoms_pair, sources_pair
    ):
        lg.set_prop("atomic numbers", atomic_nums)
        lg.set_prop("natoms slices", natoms_slices)
        lg.set_prop("sources", sources)

    return lgp


#elg = extended labeled graph
def smiles2elg(rxn_smiles,add_Hs=True,binarize=True,weight=1000):

    natoms = []
    sources = []

    ss = rxn_smiles.split('>')
    s = ss[0]+"."+ss[2]

    natoms_r = Chem.MolFromSmiles(ss[0]).GetNumAtoms()
    natoms_p = Chem.MolFromSmiles(ss[2]).GetNumAtoms()

    mol = Chem.MolFromSmiles(s)

    sources.append(mol)
    natoms.append(mol.GetNumAtoms())

    if add_Hs:
        mol = Chem.AddHs(mol)

    sources.append(mol)
    natoms.append(mol.GetNumAtoms())

    atomic_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    atom_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]

    elg = get_labeled_graph_from_mol(mol)
    if binarize:
        elg.binarize_graph()

    amn2i_r = defaultdict(list)
    for i in range(natoms_r):
        if atom_map_nums[i]>0 and atomic_nums[i]>1:
            amn2i_r[atom_map_nums[i]].append(i)

    amn2i_p = defaultdict(list)
    for i in range(natoms_r,natoms_p+natoms_r):
        if atom_map_nums[i]>0 and atomic_nums[i]>1:
            amn2i_p[atom_map_nums[i]].append(i)

    amns = set(amn for amn in amn2i_r.keys() if len(amn2i_r[amn])==1) \
         & set(amn for amn in amn2i_p.keys() if len(amn2i_p[amn])==1)

    for amn in amns:
        i = amn2i_r[amn][0]
        j = amn2i_p[amn][0]
        elg.graph[i][j] = weight
        elg.graph[j][i] = weight

    elg.set_prop("atomic numbers", atomic_nums)
    elg.set_prop("natoms slices", natoms)
    elg.set_prop("sources", sources)

    return elg


def get_numbered_rxn_smiles(rxn_smiles,map_nums_pair):
    new_smiles_pair = []
    for smi,map_nums in zip(rxn_smiles.split('>>'),map_nums_pair):
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.RWMol(mol)
        natoms = mol.GetNumAtoms()
        nnums = len(map_nums)
        if nnums<natoms:
            adjusted_map_nums = map_nums + [0]*(natoms-nnums)
        else:
            adjusted_map_nums = map_nums[:natoms]
        for atom,map_num in zip(mol.GetAtoms(),adjusted_map_nums):
            atom.SetAtomMapNum(map_num)
        new_smiles_pair.append(Chem.MolToSmiles(mol,canonical=False))
    return '>>'.join(new_smiles_pair)


def canonicalize_rxn_smiles(rxn_smiles):
    components_cano = []
    ss = rxn_smiles.split('>')
    ss.pop(1)
    for s in ss:
        mol = Chem.MolFromSmiles(s)
        mol_cano = Chem.RWMol(mol)
        mapnums = []
        for atom in mol_cano.GetAtoms():
            mapnums.append(atom.GetAtomMapNum())
            atom.SetAtomMapNum(0)
        mol_cano = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(mol_cano)))
        matches = mol.GetSubstructMatches(mol_cano)
        if matches:
            for atom, idx in zip(mol_cano.GetAtoms(), matches[0]):
                atom.SetAtomMapNum(mapnums[idx])
            s_cano = Chem.MolToSmiles(mol_cano, canonical=False,allHsExplicit=True)
        components_cano.append(s_cano)
    return '>>'.join(components_cano)

