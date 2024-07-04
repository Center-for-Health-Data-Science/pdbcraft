<h1 align="center">
<img src="./branding/pdbcraft_logo.png" width="300">
</h1><br>

## Introduction

**pdbcraft** is a lightweight Python package for parsing and manipulating PDB and PDBx/MMCIF files.

**pdbcraft** has no external dependencies and stores data in simple data structures (nested dictionaries) that can be easily modified with external tools.

With **pdbcraft**, you can:

* Build a structure from either a PDB or PDBx/MMCIF file. The structure is built using both atomic coordinates and connectivity information.
* Remove selected models, chains, segments, residues, and atoms.
* Keep only selected models, chains, segments, residues, and atoms.
* Renumber models and residues.
* Rename models, chains, segments, residues, and atoms.
* Modify attributes for models, chains, segments, residues, and atoms.
* Merge structures.
* Sort the atoms of specific residue types according to a custom order.
* Write a PDB file of your structure.
* Write a PDBx/MMCIF file of your structure.
* Write a FASTA file containing the sequences of the chains found in your structure. You can also write multiple FASTA files containing all chains in a specific model, all instances of a particular chain in all models, or each model's chain separately.
* ... all while your connectivity data gets updated accordingly if present in the original structure(s).

However, **pdbcraft** does not perform complex tasks such as structural alignments, hydrogen bond detection, or secondary structure assignment.

We recommend packages such as [BioPython](https://biopython.org/docs/1.75/api/index.html) for these operations.

## Requirements

**pdbcraft** requires Python version 3.11 or higher.

## Installation

Here, we provide a quick way to install the current version of **pdbcraft**.

A more complete installation guide will follow and will be hosted on ReadTheDocs as part of **pdbcraft**'s documentation.

To install **pdbcraft**:

* Download it as a zipped file from [GitHub](https://github.com/Center-for-Health-Data-Science/pdbcraft).

* Rename it `pdbcraft`.

* Unzip it in your preferred location.

* (Optional) Activate the Python virtual environment where you want to install it.

* Install it with ``pip`` by running:

  ```shell
  pip install ./pdbcraft
  ```

## Documentation

Detailed **documentation**, with a description of the API and full-fledged tutorials, is under construction and will soon be available on ReadTheDocs.

For now, the documentation is available as HTML files in the `pdbcraft/doc/build` directory.

## Bug reports

We welcome **reports for any bugs or problems** you may encounter in the dedicated ["issues" section](https://github.com/Center-for-Health-Data-Science/pdbcraft/issues) on GitHub.