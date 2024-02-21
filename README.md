# pyIsland
==============================

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Python library for general bioinformatics, comparative genomics and GCP wrappers

Clone the repository and install with pip:
>pip install -e .

Project Organization
------------
    │
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
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
    │   ├── data           <- Reference data
    │   │
    │   ├── drivers       <- Python wrappers around cmd line programs
    │   │
    │   ├── GCP           <- Convenience fxns for working with GCP
    │   │
    │   ├── Parsers       <- Common parsers
    │   │
    │   ├── FastaIO       <- Wrappers around BioPython for fasta manipulation
    │   │
    │   ├── Islands       <- Extracting genomic regions and annotating them
    │   │
    │   ├── tests         <- Unit tests with pytest
    │   │
    │   ├── TBD       <- TBD
    │
    └── tox.ini            <- **TBD** tox file with settings for running tox; see tox.readthedocs.io


--------

## Testing / code formating

This project comes with `pytest` and `black`. Run pytest using
>pytest pyIsland/tests/ --cov=pyIsland

And format the code using
>black pyIsland/


