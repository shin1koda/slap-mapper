from ..core import SlapMapper, simplify_labels
from .smiles import HAS_RDKIT, smiles2lgp, get_numbered_rxn_smiles
from .geom import HAS_ASE, geoms2lgp


class SlapAAM(SlapMapper):

    def map_smiles(self, rxn_smiles, add_Hs=True, break_sym="heavy", interactive=False):
        if not HAS_RDKIT:
            raise ImportError("RDKit is required for processing SMILES")

        if not self.binary:
            self.valfactor = 4

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
        label_string=None,
        base=None,
        bond_scale=1.2,
        interactive=False,
    ):
        if not HAS_ASE:
            raise ImportError("ASE is required for reading structures")

        lgp = geoms2lgp(
            files_react,
            files_prod,
            ase_read_kwargs=ase_read_kwargs,
            label_string=label_string,
            base=base,
            bond_scale=bond_scale,
        )

        targets = self._get_targets(break_sym, lgp[0].props["atomic numbers"])

        self.get_maps(lgp, break_sym_targets=targets, interactive=interactive, base=base)

    def _get_targets(self, break_sym, atomic_nums):

        if break_sym == "heavy":
            return [i for i in range(len(atomic_nums)) if atomic_nums[i] > 1]
        elif break_sym == "all":
            return list(range(len(atomic_nums)))
        else:
            return break_sym

