# pdbcraft

## Introduction

**pdbcraft** is a lightweight Python package for parsing and manipulating PDB files.

**pdbcraft** has no external dependencies and stores data in simple data structures (nested dictionaries and lists) that can be easily modified with external tools if needed.

With **pdbcraft**, you can:

* Remove models, chains, segments, residues, and atoms.
* Renumber models and residues.
* Rename chains, segments, residues, and atoms.
* Merge structures.
* Move atoms from one residue, segment, chain, or model to another place in the structure or to a new structure.
* Sort the atoms of specific residues according to a custom order.
* Write a PDB file of your structure.
* Write a FASTA file containing the sequences of the chains found in your structure. You can also write multiple FASTA files containing all chains in a specific model, all instances of a specific chain in all models, or each chain of each model separately.
* ... all while your CONECT data, if present in the original structure, get updated accordingly.

However, **pdbcraft** does not perform complex tasks such as structural alignments, hydrogen bond detection, or secondary structure assignment. For these operations, we recommend packages such as [BioPython](https://biopython.org/docs/1.75/api/index.html), which covers a broader scope.

## Documentation

Detailed **documentation** with a description of the API and full-fledged tutorials can be found [here]().

## Bug reports

We welcome **reports for any bugs or problems** you may encounter in the dedicated ["issues" section]() here on GitHub.