try:
    from ase.data import covalent_radii
    from ase.io import write, read

    HAS_ASE = True
except ImportError:
    HAS_ASE = False
import numpy as np
from ..core import LabeledGraph
from .._idxmapstr import parse_index_mapping_string


def geoms2lgp(
    files_react,
    files_prod,
    ase_read_kwargs={},
    idx_map_str=None,
    base=0,
    bond_scale=1.2,
):

    files_pair = []
    for files in [files_react, files_prod]:
        if isinstance(files, list):
            files_pair.append(files)
        else:
            files_pair.append([files])

    lgp = []
    atomic_nums_pair = []
    natoms_slices_pair = []
    sources_pair = []
    graph_pair = []
    labels_pair = []

    for files in files_pair:
        atomic_nums = []
        natoms_slices = [0]
        sources = []
        graph = {}
        for file in files:
            atoms = read(file, **ase_read_kwargs)
            sources.append(atoms)
            temp_graph, temp_atomic_nums, temp_natoms = atoms2lg_inps(
                atoms, bond_scale, idx_shift=natoms_slices[-1]
            )
            natoms_slices.append(temp_natoms + natoms_slices[-1])
            atomic_nums.extend(temp_atomic_nums)
            graph = {**graph, **temp_graph}

        atomic_nums_pair.append(atomic_nums)
        labels_pair.append(atomic_nums.copy())
        natoms_slices_pair.append(natoms_slices)
        sources_pair.append(sources)
        graph_pair.append(graph)

    if sorted(atomic_nums_pair[0]) != sorted(atomic_nums_pair[1]):
        raise ValueError(f"Unbalanced reactions are not supported in geometric mapping")

    if idx_map_str is not None:
        pairs = parse_index_mapping_string(idx_map_str, base=base)

        for i, pair in enumerate(pairs):
            for idxs, labels in zip(pair, labels_pair):
                for idx in idxs:
                    labels[idx] = i + 1000

    if sorted(labels_pair[0]) != sorted(labels_pair[1]):
        raise ValueError(f"Unbalanced label assignment is detected")

    lgp = [
        LabeledGraph(graph, labels) for graph, labels in zip(graph_pair, labels_pair)
    ]

    for lg, atomic_nums, natoms_slices, sources in zip(
        lgp, atomic_nums_pair, natoms_slices_pair, sources_pair
    ):
        lg.set_prop("atomic numbers", atomic_nums)
        lg.set_prop("natoms slices", natoms_slices)
        lg.set_prop("sources", sources)

    return lgp

def geoms2lg(
    files,
    ase_read_kwargs={},
    bond_scale=1.2,
):

    if not isinstance(files, list):
        files = [files]

    atomic_nums = []
    natoms_slices = [0]
    sources = []
    graph = {}
    for file in files:
        atoms = read(file, **ase_read_kwargs)
        sources.append(atoms)
        temp_graph, temp_atomic_nums, temp_natoms = atoms2lg_inps(
            atoms, bond_scale, idx_shift=natoms_slices[-1]
        )
        natoms_slices.append(temp_natoms + natoms_slices[-1])
        atomic_nums.extend(temp_atomic_nums)
        graph = {**graph, **temp_graph}

    labels = atomic_nums.copy()

    lg = LabeledGraph(graph, labels)

    lg.set_prop("atomic numbers", atomic_nums)
    lg.set_prop("natoms slices", natoms_slices)
    lg.set_prop("sources", sources)

    return lg

def atoms2lg_inps(atoms, bond_scale, idx_shift=0):

    natoms = len(atoms)
    atomic_nums = [int(i) for i in atoms.get_atomic_numbers()]

    d = atoms.get_all_distances()
    r = covalent_radii[atomic_nums]
    A = (d / (r[None, :] + r[:, None])) < bond_scale
    np.fill_diagonal(A, False)

    graph = {}
    for i in range(natoms):
        atom_neighbors = {}
        for j in range(natoms):
            if A[i, j]:
                atom_neighbors[j + idx_shift] = 1
        graph[i + idx_shift] = atom_neighbors

    return graph, atomic_nums, natoms


