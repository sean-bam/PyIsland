#!/bin/bash python3

"""
Constructs a HHsuite DB from a dir of MSAs
"""

# Imports --------------------------------------------------------------------------------------------------
import argparse
from pathlib import Path
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
import shutil

# Args -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input", help="path/to/input", type=str, required=True)

parser.add_argument(
    "-db", "--dbname", help="name of output database", type=str, required=True
)

parser.add_argument(
    "-s",
    "--suffix",
    help="suffix of fasta files. Default = '.afa'",
    type=str,
    required=False,
)

# parser.add_argument('-k',
#                    '--keep_msa_names',
#                    help="dont rename first seq in MSA to input filename",
# type = bool,
#                    action = 'store_true',
#                    required=False)

args = parser.parse_args()

# Constants ------------------------------------------------------------------------------------------------
if args.suffix:
    suffix = args.suffix
else:
    suffix = "afa"


# Functions ------------------------------------------------------------------------------------------------
def fix_fasta_for_hhsuite(fasta, output):
    """
    1. changes the header of the first sequence in a fasta file to the filename
    2. replaces "*" with gaps, because * is ignored by HHsuite.
    Thus, if it isn't present in ALL seqs, HHsuite thinks the MSA is unaligned
    """
    infile = Path(fasta)
    outfile = Path(output)
    i = 1

    with outfile.open("w") as o:
        with infile.open() as handle:
            for title, seq in SimpleFastaParser(handle):
                if i == 1:
                    title = infile.stem
                if "*" in seq:
                    seq = seq.replace("*", "-")
                print(f">{title}", seq, sep="\n", file=o)
                i += 1


def sort_db():
    """
    code to sort the HHsuite DB according to length,
    as suggested in the wiki

    TBD: the HHM index/DB throws an error when trying to sort
    So skipping this
    """
    sort_dat = tmp_dir / Path("sorting.dat")

    hhm_index_sorted = Path(f"{args.dbname}_hhm_ordered.ffindex")
    hhm_data_sorted = Path(f"{args.dbname}_hhm_ordered.ffdata")

    a3m_index_sorted = Path(f"{args.dbname}_a3m_ordered.ffindex")
    a3m_data_sorted = Path(f"{args.dbname}_a3m_ordered.ffdata")

    subprocess.run(
        f"sort -k3 -n -r {args.dbname}_cs219.ffindex | cut -f1 > {sort_dat}",
        shell=True,
        check=True,
    )

    # if not hhm_index_sorted.is_file() and not hhm_data_sorted.is_file():
    #    subprocess.run(f'ffindex_order {sort_dat} {hhm_data} {hhm_index} {hhm_data_sorted} {hhm_index_sorted}',
    #           shell = True,
    #           check = True
    #          )

    if not a3m_index_sorted.is_file() and not a3m_data_sorted.is_file():
        subprocess.run(
            f"ffindex_order {sort_dat} {a3m_data} {a3m_index} {a3m_data_sorted} {a3m_index_sorted}",
            shell=True,
            check=True,
        )

    # hhm_index_sorted.replace(hhm_index)
    # hhm_data_sorted.replace(hhm_data)

    a3m_index_sorted.replace(a3m_index)
    a3m_data_sorted.replace(a3m_data)

    sort_dat.unlink()


# -----------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    # make a temporary directory to write files
    tmp_dir = Path("hhdb_dir/")
    if not tmp_dir.is_dir():
        tmp_dir.mkdir()

    # Change the first name of the sequence in each fasta file
    # When we make the HHsuite DB, this will be the name
    # rather than "consensus" or something else
    for msa in Path(args.input).glob(f"*.{suffix}"):
        # if not args.keep_msa_names:
        fix_fasta_for_hhsuite(msa, tmp_dir / msa.name)
        # else:
        # shutil.copyfile(msa, tmp_dir / msa.name)

    # convert to a3m-format. This is necessary b/c fasta DBs throw segfaults
    # trim header below 33 characters (also causes segfaults)
    if suffix != "a3m":
        for msa in tmp_dir.glob(f"*.{suffix}"):
            subprocess.run(
                f"reformat.pl fas a3m {msa} .a3m -d 33", shell=True, check=True
            )

            # remove fasta file if successful
            assert msa.with_suffix(
                ".a3m"
            ).is_file(), f"something went wrong trying to reformat {msa} to a3m"
            msa.unlink()

    # make a ffindex/data file
    a3m_index = Path(f"{args.dbname}_a3m.ffindex")
    a3m_data = Path(f"{args.dbname}_a3m.ffdata")

    if not a3m_data.is_file() and not a3m_index.is_file():
        print("making MSA ffindex/data")
        subprocess.run(
            f"ffindex_build -s {args.dbname}_a3m.ffdata {args.dbname}_a3m.ffindex {tmp_dir}",
            shell=True,
            check=True,
        )
    # make HMMs
    hmm_dir = tmp_dir / Path("hhms/")
    if not hmm_dir.is_dir():
        hmm_dir.mkdir(parents=True)

    for file in tmp_dir.glob(f"*.a3m"):
        hhm = hmm_dir / f"{file.stem}.hhm"
        if not hhm.exists():
            subprocess.run(
                f"hhmake -i {file} -o {hhm} -name {file.stem} -v 0 -diff 1000",
                shell=True,
                check=True,
            )
            assert (
                hhm.exists()
            ), f"something went wrong trying to build the HHM of {file}"

    # make a ffindex from the HMMs
    hhm_data = Path(f"{args.dbname}_hhm.ffdata")
    hhm_index = Path(f"{args.dbname}_hhm.ffindex")

    if not hhm_data.is_file() and not hhm_index.is_file():
        print("making HHM indexes")

        subprocess.run(
            f"ffindex_build -s {args.dbname}_hhm.ffdata {args.dbname}_hhm.ffindex {hmm_dir}",
            shell=True,
            check=True,
        )

    # run cstranslate
    if (
        not Path(f"{args.dbname}_cs219.ffdata").exists()
        and not Path(f"{args.dbname}_cs219.ffindex").exists()
    ):
        print("running cstranslate")
        subprocess.run(
            f"cstranslate -f -x 0.3 -c 4 -I a3m -i {args.dbname}_a3m -o {args.dbname}_cs219",
            shell=True,
            check=True,
        )

    # sort. This step seems to be optional
    # Strangely, the HHM ffindex/database throws an error
    sort_db()

    # cleanup
    for file in tmp_dir.rglob("*.a3m"):
        file.unlink()
    for file in tmp_dir.rglob("*.hhm"):
        file.unlink()
    hmm_dir.rmdir()
    tmp_dir.rmdir()
