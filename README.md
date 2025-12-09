# SLAPMapper: A scalable Python implementation of general-purpose atom-to-atom mapping

## Overview

**SLAPMapper** is a Python package for **atom-to-atom mapping (AAM)** tool based on a fast,
Weisfeiler–Lehman-like iterative label refinement coupled with
sequential linear assignment problems (SLAPs).  
This approach yields **high accuracy across diverse chemical domains**
while remaining **scalable to reactions involving thousands of atoms**.

The tool is designed for researchers who require reliable AAM tools
for reaction templating, mechanistic analysis, quantum-chemical simulations, and so on.  
Beyond chemistry, the core algorithm can also serve as a general-purpose
**approximate graph-matching framework**.

## Features

- Chemical-distance (CD) minimization
- WL-like iterative label refinement
- Sequential linear assignment formulation
- Enumeration of alternative mappings
- Partial constraints on AAMs
- RDKit-based SMILES and ASE-based 3D structure handling


## Requirements

- Python 3.8+
- NumPy
- SciPy
- RDKit (optional, for handling SMILES)
- ASE (optional, for handling 3D-structure)


## Installation

### Installing with RDKit

We generally recommend installing this package via **conda**, as RDKit is most reliably installed through conda.

```bash
conda create -n slap python=3.10
conda activate slap
conda install -c conda-forge rdkit
pip install "slapmapper @ git+https://github.com/shin1koda/slap-mapper"
```

### Installing with ASE only

If you only need ASE support, a pure pip installation is sufficient:

```bash
pip install "slapmapper[ase] @ git+https://github.com/shin1koda/slap-mapper"
```

### Installing with both RDKit and ASE

Use conda to install RDKit and ASE, then install SLAPMapper via pip:

```bash
conda create -n slap python=3.10
conda activate slap
conda install -c conda-forge rdkit ase
pip install "slapmapper @ git+https://github.com/shin1koda/slap-mapper"
```


## Example usage

An example script is provided in `sample/sample.py` of [SLAPMapper GitHub repository](https://github.com/shin1koda/slap-mapper):

```python
from slap_mapper.aam import SlapAAM

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
results = mapper.results
```


## Documentation

For more details, please refer to the [API documentation](https://shin1koda.github.io/slap-mapper/).


## Citation

 1. S.-i. Koda and  S. Saito, General and scalable atom-to-atom mapping via Weisfeiler–Lehman-like approximate graph matching, ChemRxiv (2025). [doi:10.26434/chemrxiv-2025-hthwn](https://doi.org/10.26434/chemrxiv-2025-hthwn)

If you use this software, please cite **Ref. 1**.
  

## Community guidelines

### Contributing

Contributions to this project are welcome. If you would like to contribute new features, improvements, or documentation, please open a pull request on GitHub.
Before submitting a PR, we recommend opening a short issue to discuss the proposed change.

### Reporting issues

If you encounter a problem, unexpected behavior, or a potential bug, please report it through the GitHub issue tracker:

https://github.com/shin1koda/slap-mapper/issues

When reporting an issue, please include:
- A clear description of the problem
- Steps to reproduce the issue
- Your environment (Python version, RDKit/ASE version, etc.)
- Any relevant error messages or logs

### Seeking support

If you have questions about the usage of the package, or need help integrating it into your workflow, feel free to open an issue labeled “question” on GitHub.
We will do our best to provide guidance based on availability.


## License

This project is distributed under the MIT License.