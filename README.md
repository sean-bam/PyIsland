# pyIsland
==============================

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A python library for comparative genomics of microbial "islands"

# Installation

Clone the repository and install with pip:
>pip install -e .

Project Organization
------------
    │
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── docs               <- *TBD* A default Sphinx project; see sphinx-doc.org for details
    │
    │── bin               <- Miscellaneous standalone scripts
    │
    ├── requirements.txt   <- The requirements file for building the libary
    │
    ├── setup.py           <- makes project pip installable
    │
    ├── pyIsland                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes pyIsland a Python module
    │   │
    │   │
    │   ├── Drivers       <- Python wrappers around cmd line programs
    │   │
    │   ├── GCP           <- Convenience fxns for working with GCP
    │   │
    │   ├── Parsers       <- Common parsers
    │   │
    │   ├── FastaUtils     <- Wrappers around BioPython for fasta manipulation
    │   │
    │   ├── Islands       <- Extracting genomic regions and annotating them
    │   │
    │   ├── tests         <- Unit tests with pytest
    │   │
    │   ├── TBD       <- TBD

--------

# Usage

TBD

# Testing / code formating

This project comes with `pytest` and `black`. Run pytest using
>pytest pyIsland/tests/ --cov=pyIsland

And format the code using
>black pyIsland/


