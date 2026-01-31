---
title: 'SLAPMapper: A Python package for general and scalable atom-to-atom mapping'
tags:
  - Python
  - chemistry
  - cheminformatics
  - graph matching
  - atom-to-atom mapping
  - AAM
authors:
  - name: Shin-ichi Koda
    orcid: 0000-0003-0993-3678
    affiliation: "1, 2"
    corresponding: true
  - name: Shinji Saito
    orcid: 0000-0003-4982-4820
    affiliation: "1, 2"
    corresponding: true
affiliations:
  - name: Institute for Molecular Science, National Institutes of Natural Sciences, Okazaki, 444-8585, Japan
    index: 1
  - name: Graduate University for Advanced Studies, SOKENDAI, Okazaki, 444-8585, Japan
    index: 2

date: 10 December 2025
bibliography: paper.bib
---

# Summary

Atom-to-atom mapping (AAM) identifies the atomic correspondence between reactants and products and is essential for reaction templating, mechanistic analysis, and data-driven or quantum-chemical workflows. More broadly, AAM is a domain-specific instance of graph matching, a fundamental computational task that arises across many scientific fields. SLAPMapper is a Python package that implements a recently developed polynomial-time approximation to chemical-distance minimization, formulated through sequential linear assignment with iterative label refinement. The underlying algorithm is provided both as a general framework for matching labeled graphs and as a specialized AAM implementation built on top of it [@koda2025generala; @slapmapper]. This approach preserves the generality of optimization-based AAM while achieving runtime performance comparable to machine-learning methods and scaling to systems with thousands of atoms. SLAPMapper supports both SMILES and 3D structural inputs, can enumerate symmetrically distinct minimal mappings, and accepts optional partial constraints, making it a versatile and efficient tool for reaction informatics, synthesis planning, and automated computational chemistry.


# Statement of need

Graph matching is a fundamental problem in computational science, supporting tasks in pattern recognition, computer vision, social-network analysis, and bioinformatics. Many workflows rely on identifying correspondences between two similar graphs, and efficient, scalable matching algorithms therefore have broad methodological impact across scientific domains.

In chemistry, atom-to-atom mapping (AAM) represents a specialized yet widely used instance of graph matching. By determining the atomic correspondence between reactants and products, AAM enables reaction-template extraction, mechanistic interpretation, similarity search, construction of high-quality training data for machine-learning models, and quantum chemical simulations.

However, existing AAM tools face a persistent trade-off between generality and runtime. Rule-based and, in particular, optimization-based approaches are broadly applicable but often slow. Conversely, machine-learning-based (ML-based) methods offer high accuracy and fast inference within trained domains, yet their performance degrades for unseen chemistry.

Our recent work introduced a theory-driven algorithm that resolves this trade-off by relaxing chemical-distance minimization into sequential linear assignment with Weisfeiler–Lehman-like label refinement [@koda2025generala]. Chemical distance measures how many bonds must be broken or formed to transform one molecular graph into another; minimizing it provides a physically reasonable way to identify chemically plausible mappings. By approximating this minimization efficiently, our method preserves the generality of optimization-based approaches while achieving runtimes comparable to ML-based tools and scaling to systems with thousands of atoms. Building on this theoretical foundation, SLAPMapper [@slapmapper] provides an accessible Python implementation that enables cheminformatics researchers, computational chemists, and data scientists to perform accurate and scalable AAM without domain-specific training data.


# Key features

SLAPMapper provides a unified interface for performing AAM from both 2D and 3D molecular representations.

* **Support for SMILES-based mapping via RDKit.**  
  RDKit parses SMILES strings, which are widely used linear notations for molecular graphs. SLAPMapper then uses this parsed information to construct its labeled graph representation. This allows the tool to integrate smoothly into existing SMILES-based workflows.

* **Support for 3D structures via the Atomic Simulation Environment (ASE).**  
  ASE allows SLAPMapper to accept Cartesian molecular geometries, infer bonding from interatomic distances, and handle systems that are difficult or ambiguous to represent using SMILES, including coordination complexes, transition-metal species, and other inorganic reactions.

* **Enumeration of symmetrically distinct CD-minimal mappings.**  
  The algorithm identifies symmetrically distinct solutions during the sequential linear assignments, making it possible to enumerate multiple minimal-distance mappings rather than returning only a single candidate. Such enumeration is useful for detecting potential side reactions.

* **Partial constraints on mappings.**  
  Users may pre-specify correspondences for selected atoms, for example by providing partially annotated SMILES. SLAPMapper incorporates these constraints as fixed labels and completes the remainder of the mapping automatically. This feature can substantially reduce manual annotation effort when constructing training data for ML-based AAM tools.

* **General graph-matching framework.**  
  The core algorithm is implemented as a general framework for matching labeled graphs, and the AAM functionality is built as a specialization of this framework. This design allows SLAPMapper to be used not only for chemical reactions but also for broader graph-matching tasks that require scalable approximate alignment.


# Related Works

A wide variety of AAM tools have been proposed. Rule-based approaches include Automapper [@automapper], Indigo [@pavlov2011indigoun], RDTool [@rahman2016reaction], NameRXN [@namerxn], and recent optimization-based formulations such as the Ising model chemical-distance optimizer [@ali2025enumerat]. ML-based methods, including RXNMapper [@schwaller2021extracti], GraphormerMapper [@nugmanov2022bidirect], LocalMapper [@chen2024precisea], and SAMMNet [@astero2025enhancin], have demonstrated high mapping accuracy and fast inference. LocalMapper additionally incorporates a human-in-the-loop framework, and SAMMNet provides an alternative deep-learning architecture, although its trained parameters are not publicly available. Benchmarks have shown that RXNMapper significantly outperforms rule-based tools in speed while achieving high accuracy within its training domain [@lin2022atomtoat], and that GraphormerMapper and LocalMapper provide further accuracy improvements.

SLAPMapper [@slapmapper] was benchmarked against representative ML-based AAM tools and showed comparable accuracy within their trained chemical domain and substantially higher accuracy for unseen chemistry, while achieving comparable or better runtime performance and substantially better scalability [@koda2025generala]. These results highlight the benefit of a theory-driven and training-free approach that preserves chemical generality while operating efficiently on large molecular graphs.


# Acknowledgements

This work has been supported by JSPS KAKENHI, Grant Number JP22K14652 (S-i.K.), JP21H04676, and JP23K17361 (S.S.). The software development and benchmarking were partially performed using Research Center for Computational Science, Okazaki, Japan (Project: 24-IMS-C193 and 25-IMS-C223).

# References
