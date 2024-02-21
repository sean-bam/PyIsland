from pathlib import Path

import pandas as pd


def mmseqs_result2flat_decompose(flatfile, output_dir, suffix="faa"):
    """
    expects the output of mmseqs result2flat, which looks like:

    >Q0KJ32
    >Q0KJ32
    MAGA....R
    >C0W539
    MVGA....R
    >E3HQM9
    >E3HQM9
    MCAT...Q
    >Q223C0
    MCAR...Q

    Converts each cluster into an individual fasta file named after the representative

    fasta files are overwritten by default
    """

    if not Path(output_dir).is_dir():
        Path(output_dir).mkdir(parents=True)

    with open(flatfile) as infile:
        lines = infile.readlines()

        for i in range(0, len(lines)):
            # if the line after the current one doesn't have seq data,
            # its a representative.
            if lines[i].startswith(">") and lines[i + 1].startswith(">"):
                # get the first word of the header
                header = lines[i].strip().split(">")[1]

                # make an output file
                output_file = Path(output_dir) / Path(f"{header}.{suffix}")

                # overwrite existing output files
                if output_file.exists():
                    output_file.unlink()

                continue

            else:
                with output_file.open("a") as o:
                    print(lines[i], end="", file=o)


def mmseqs_msadb_split(msaDB, output_dir, suffix=".a3m"):
    """
    Splits each MSA in an MMSeqs msaDB into individual files
    MMSeqs separates MSAs with a null byte, starting with the 2nd MSA

    Note: I tried
    ffindex_unpack msaDB_ca3m.ffdata msaDB_ca3m.ffindex tmp/
    But the output MSAs are in binary format and I'm not sure they work with HHsearch
    """

    data = []
    first_msa = True
    i = 0

    p = Path(output_dir)
    if not p.is_dir():
        p.mkdir(parents=True)

    with open(msaDB, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if first_msa:
                if line.startswith(">"):
                    header = line.strip().split()[0].replace(">", "")

                    if i == 0:
                        # representative doesn't have a null byte in header
                        output_file = p / Path(header + suffix)

                        # overwrite if file exists
                        if output_file.is_file():
                            output_file.unlink()

                        # increment the counter
                        i += 1

                elif line.startswith("\x00"):
                    # we hit a new cluster
                    first_msa = False
                    header = line.strip().split("\x00")[1].replace(">", "")
                    output_file = p / Path(header + suffix)

                    # overwrite if file exists
                    if output_file.is_file():
                        output_file.unlink()

                    # print the header w/o the null byte
                    with output_file.open("a") as o:
                        print(line.strip().split("\x00")[1], file=o)

                    # exit to all other clusters
                    continue

                    # i = 2

                with output_file.open("a") as o:
                    print(line.strip(), file=o)

            # all other clusters
            else:
                if line.startswith("\x00"):
                    header = line.strip().split("\x00")[1].replace(">", "")
                    output_file = p / Path(header + suffix)

                    # overwrite if file exists
                    if output_file.is_file():
                        output_file.unlink()

                    # print the header w/o the null byte
                    with output_file.open("a") as o:
                        print(line.strip().split("\x00")[1], file=o)

                else:
                    with output_file.open("a") as o:
                        print(line.strip(), file=o)

    # not sure why this file gets made, but remove
    badfile = p / Path(suffix)
    badfile.unlink()


def mmseqs_lookup_to_dict(lookup):
    """ """
    d = {}
    i = 0
    headers = set()
    with open(lookup) as infile:
        for line in infile:
            mmseqs_internal_id, header, drop = line.strip().split()
            d[header] = mmseqs_internal_id
            i += 1
    assert len(d) == i, f"Your lookup file has duplicate headers"

    return d


def mmseqs_linear_tsv_to_df(tsv):
    # df_cls = pd.read_csv(tsv, sep="\t", names=["cluster", "members"])
    df_cls = mmseqs_tsv_to_df(tsv)

    # Explode the members column
    df_cls["members"] = df_cls.members.str.split(" ").to_list()
    df_cls = df_cls.explode("members")

    return df_cls


def mmseqs_tsv_to_df(tsv):
    df_cls = pd.read_csv(tsv, sep="\t", names=["cluster", "members"])

    ## Explode the members column
    # df_cls["members"] = df_cls.members.str.split(" ").to_list()
    # df_cls = df_cls.explode("members")

    return df_cls


def cls2fa(fasta, member_to_cluster, outfolder):
    """
    Expects a fasta file and a dictionary of member (key) to rep/cluster (value)
    writes member sequences in fasta file named after the rep
    """
    with open(fasta) as ff:
        for seq in SeqIO.parse(ff, "fasta"):
            if seq.id in member_to_cluster:
                mapping = member_to_cluster[seq.id]
                cluster_fasta_name = outfolder + "/" + mapping + ".faa"

                with open(cluster_fasta_name, "a") as fasta_file:
                    fasta_file.write(seq.format("fasta"))


def get_mmseqs_internal_ids_as_set(index):
    mmseqs_internal_ids = set()
    with open(index) as infile:
        for line in infile:
            internal_id, offset, size = line.strip().split()
            mmseqs_internal_ids.add(internal_id)

    return mmseqs_internal_ids


def get_mmseqs_internal_ids_as_dict(lookup):
    d = {}
    with open(lookup) as infile:
        for line in infile:
            internal_id, header, offset = line.strip().split()
            d[internal_id] = header

    return d


def mmseqs_tsv_to_dict(tsv):
    d = {}
    with open(tsv) as infile:
        for line in infile:
            rep, member = line.strip().split()
            d[member] = rep
    return d
