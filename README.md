# SLAPMapper

**SLAPMapper** is a Python package implementing a **general and scalable atom-to-atom mapping (AAM)** algorithm based on **Weisfeiler–Lehman-like approximate graph matching**.  

---

## Features

- Based on chemical distance (CD) minimization
- WL-like iterative label refinement  
- Sequential Linear Assignment Problems (LAPs)  
- Enumeration of solutions
- RDKit-based SMILES and ASE-based 3D-structure handling  
- Scalable to large molecules (e.g., polypeptides with 8,000+ atoms)  

---

## Installation

### From source
```
git clone https://github.com/YOURNAME/slap-mapper.git
cd slap-mapper
conda create -n slap python=3.10
conda activate slap
pip install -e .
```

### Dependencies
- Python 3.8+
- NumPy / SciPy  
- RDKit (optional, for handling SMILES)
- ASE (optional, for handling 3D-structure)

---

## Quick Example

```
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

## License

This project is licensed under the **MIT License**.

