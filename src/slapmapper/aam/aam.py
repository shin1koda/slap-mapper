from ..core import SlapMapper
from ._smiles import HAS_RDKIT, smiles2lgp, get_numbered_rxn_smiles
from ._geom import HAS_ASE, geoms2lgp
from .._idxmapstr import lgp2idx_map_str


class SlapAAM(SlapMapper):
    """
    Sequential LAP-based approximate mapper for a pair of labeled graphs.

    The algorithm reduces the problem of matching two labeled graphs to a
    sequence of label refinement steps coupled with solving a series of
    linear assignment problems (LAPs). Labels are progressively refined
    using neighborhood information, and branching is introduced when multiple
    symmetrically distinct LAP solutions are available. Symmetry-breaking
    heuristics can be applied to choose a representative mapping among
    multiple symmetrically equivalent ones.

    Parameters
    ----------
<<<<<<< HEAD
    binary : bool, optional (default=False)
=======
    binary : bool, optional (default=True)
>>>>>>> joss
        If True, all edge weights are converted to 0/1 prior to evaluation.
        This yields a purely structural matcher based only on adjacency.

    Attributes
    ----------
    results : list[dict]
        List of matching results. Each dictionary contains:
            - ``"lgp"`` : the pair of refined LabeledGraph objects
            - ``"val"`` : cumulative LAP cost, which coincides with the original graph matching cost when all symmetries are broken.
            - ``"lap_sols"`` : LAP solutions for each label
            - ``"cd"`` : chemical distance computed from cumulative LAP cost
            - ``"base"`` : (if interactive) base index for user prompts
            - ``"choices"`` : (if interactive) record of the user-selected mapping, written in :ref:`the index mapping string<idxmapstr>` format
            - ``"smiles"`` : (if mapped from SMILES) mapped reaction SMILES
            - ``"mapping"`` : (if mapped from 3D structures) mapping in :ref:`the index mapping string<idxmapstr>` format
    minval : int
        Minimum cumulative LAP cost encountered.
    binary : bool
        Whether binary edge processing is used.

    """

    def __init__(self, binary=True):

        super().__init__(binary)
        self._valfactor = 2

    def get_maps(self, lgp, break_sym_targets=None, interactive=False, base=None):
        """
        Almost the same as :meth:`SlapMapper.get_maps` but each element in 
        ``self.results`` additionally contains the chemical distance ('cd')
        value computed from the cumulative LAP cost ('val').
        Users may use :meth:`map_smiles` or :meth:`map_3d` methods instead of
        calling this method directly.
        """
        super().get_maps(lgp, break_sym_targets, interactive, base)

        for r in self.results:
            if r['val']%self._valfactor==0:
                r['cd'] = int(r['val']//self._valfactor)
            else:
                r['cd'] = r['val']/self._valfactor

    def map_smiles(self, rxn_smiles, add_Hs=True, break_sym="heavy", interactive=False):
        """
        Compute AAMs from a reaction SMILES string.

        RDKit is used to convert the SMILES representation of a reaction
        into a pair of labeled graphs. SLAP-based matching is then applied,
        followed by optional symmetry-breaking.

        Parameters
        ----------
        rxn_smiles : str
            Reaction SMILES string of the form ``"A>>B"``. Partial mappings
            annotated using atom map numbers (e.g. ``"[C:1]"``) are allowed
            for imposing constraints on the mapping.
        add_Hs : bool, optional (default=True)
            If True, explicit hydrogens are added before graph construction.
        break_sym : {"heavy", "all"} or list[int], optional (default="heavy")
            Strategy for symmetry-breaking:
            - ``"heavy"`` : break symmetry only on non-hydrogen atoms  
            - ``"all"`` : break symmetry on all atoms  
            - list[int] : explicit atom indices of the reactant graph
        interactive : bool, optional (default=False)
            Enable interactive selection of symmetry-breaking choices.

        Notes
        -----
        - Requires RDKit.
        - After mapping, each result dictionary contains::
            
              r["smiles"] : str
                  Reaction SMILES annotated with mapped atom indices.
        """

        if not HAS_RDKIT:
            raise ImportError("RDKit is required for processing SMILES")

        if not self.binary:
            self._valfactor = 4

        lgp = smiles2lgp(rxn_smiles,add_Hs=add_Hs)

        targets = self._get_targets(break_sym, lgp[0].props["atomic numbers"])

        if interactive:
            natoms = len(lgp[0].labels)
            idxs_1based = list(range(1, natoms + 1))
            smis = get_numbered_rxn_smiles(
                rxn_smiles, [idxs_1based, idxs_1based]
            ).split(">>")
            print("Reaction SMILES with 1-based indexes")
            print(smis[0])
            print(">>")
            print(smis[1])
            print()

        self.get_maps(lgp, break_sym_targets=targets, interactive=interactive, base=1)
        for r in self.results:
            r["smiles"] = get_numbered_rxn_smiles(
                rxn_smiles, [r["lgp"][0].labels, r["lgp"][1].labels]
            )

    def map_3d(
        self,
        files_react,
        files_prod,
        ase_read_kwargs={},
        break_sym="heavy",
        constraints=None,
        base=0,
        bond_scale=1.2,
        interactive=False,
    ):
        """
        Compute AAMs from 3D molecular structures.

        The input structures are read using ASE, converted into labeled graphs.
        SLAP-based matching is then applied, followed by optional symmetry-breaking.

        Parameters
        ----------
        files_react : list[str] or str
            Filename(s) of the reactant structure(s).
        files_prod : list[str] or str
            Filename(s) of the product structure(s).
        ase_read_kwargs : dict, optional (default={})
            Keyword arguments forwarded to ``ase.io.read``.
        break_sym : {"heavy", "all"} or list[int], optional (default="heavy")
            Symmetry-breaking strategy (same as in :meth:`map_smiles`).
        constraints : str or None, optional (default=None)
            :ref:`The index mapping string<idxmapstr>` for imposing constraints on the mapping.
        base : int or None , optional (default=0)
            Base index for atom numbering in ``constraints`` and interactive prompts.
        bond_scale : float, optional (default=1.2)
            Multiplier used when determining adjacency from interatomic
            distances based on covalent radii.
        interactive : bool, optional (default=False)
            Enable interactive selection of symmetry-breaking choices.

        Notes
        Notes
        -----
        - Requires ASE.
        - After mapping, each result dictionary contains::

              r["mapping"] : str
                    Mapping in the index mapping string format.

        """
        if not HAS_ASE:
            raise ImportError("ASE is required for reading structures")

        lgp = geoms2lgp(
            files_react,
            files_prod,
            ase_read_kwargs=ase_read_kwargs,
            idx_map_str=constraints,
            base=base,
            bond_scale=bond_scale,
        )

        targets = self._get_targets(break_sym, lgp[0].props["atomic numbers"])

        self.get_maps(lgp, break_sym_targets=targets, interactive=interactive, base=base)
        for r in self.results:
            r["mapping"] = lgp2idx_map_str(r["lgp"],base)

    def _get_targets(self, break_sym, atomic_nums):

        if break_sym == "heavy":
            return [i for i in range(len(atomic_nums)) if atomic_nums[i] > 1]
        elif break_sym == "all":
            return list(range(len(atomic_nums)))
        else:
            return break_sym

