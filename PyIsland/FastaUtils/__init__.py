from pathlib import Path

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

from . import MSA


def subseq_fasta_biopython(iterable, infasta, outfasta):
    """
    Selects fasta records from one file using a python iterable
    """
    records = (r for r in SeqIO.parse(infasta, "fasta") if r.id in iterable)

    outfile = Path(outfasta)
    SeqIO.write(records, outfile, "fasta")


def slice_fasta_biopython(contig, start, stop, infasta, outfasta):
    """
    Accepts 1-based nucleotide coordinates and a contig ID
    slices the input fasta file to ouput fasta file
    """
    # adjust to 0-based for biopython
    start = start - 1

    seq_record = False

    records = SeqIO.parse(infasta, "fasta")
    for record in records:
        if record.id == contig:
            seq_record = record[start:stop]

    assert (
        seq_record and len(seq_record) > 0
    ), f"""
    could not find {contig} in {infasta}
    """
    SeqIO.write(seq_record, outfasta, "fasta")


def count_seqs_in_fasta(fasta):
    i = 0
    with open(fasta) as infile:
        for line in infile:
            if line.startswith(">"):
                i += 1
    return i


def deduplicate_fasta(fasta, output):
    headers = set()
    seqs = set()

    nr_count = 0
    all_count = 0
    with open(fasta) as infile:
        with open(output, "w") as outfile:
            for title, seq in SimpleFastaParser(infile):
                all_count += 1

                if title not in headers:
                    headers.add(title)
                    print(
                        f">{title}",
                        seq,
                        sep="\n",
                        file=outfile,
                    )

                    nr_count += 1
    print(f"Wrote {nr_count} of {all_count} seqs to {output}")


def fa_strict(fasta_in, fasta_out, min_len=0):
    """
    - Discards any sequences with X
    - Removes asterisk stop codon symbol
    - Removes duplicate fasta entries based on header
    """
    seq_records = []
    i = 0
    headers = set()

    for seq_record in SeqIO.parse(fasta_in, "fasta"):
        i += 1

        if "X" in seq_record.seq.upper():
            continue

        if seq_record.id in headers:
            continue

        if len(seq_record.seq) < min_len:
            continue

        if "*" in seq_record.seq:
            seq_record.seq = seq_record.seq.replace("*", "")

        headers.add(seq_record.id)
        seq_records.append(seq_record)

    with open(fasta_out, "w") as o:
        SeqIO.write(seq_records, fasta_out, "fasta")
        print(f"wrote {len(seq_records)}/{i} sequences")


def get_subseq_coords(query, seq_record):
    """
    reports 1-based start/stop/strand of a query subsequence in a subject seq
    query subsequence is a string
    subject is BioPython SeqRecord

    Returns a tuple of start,stop, strand
    Start = -1 if the subsequence is not found
    """

    qlen = len(query)

    # search the positive strand
    start = seq_record.seq.find(query)
    stop = start + qlen
    strand = "1"

    # if no hit, search the RC
    if start == -1:
        # reverse complement the query
        query = str(Seq(query).reverse_complement())

        # search
        start = seq_record.seq.find(query)
        stop = start + qlen
        strand = "-1"

    # no hits on either strand
    if start == -1:
        start = -2
        stop = -1
        print(f"could not find the qseq on either + or - strand")

    # update to 1-based
    start += 1
    return start, stop, strand
