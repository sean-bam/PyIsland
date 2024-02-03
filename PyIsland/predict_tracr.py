#!/usr/bin/env python
"""
Accepts a Genbank-formatted file with features annotated (ORFs as "CDS" and CRISPRs as "repeat_region")
The ORFs/CRISPRs must be sorted by start position
Identifies candidate TRACR sequences using BLASTn and RNAFold from the ViennaRNA package
Outputs a CSV-formatted table containing likely locations of tracRNAs
"""
# Imports --------------------------------------------------------------------------------------------------

import subprocess
from pathlib import Path
import argparse

# add this folder to my python path, which holds the ViennaRNA
import sys

sys.path.append("/usr/local/lib/python3.7/site-packages/")

import RNA
import forgi

import pandas as pd

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Args -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input", help="/path/to/genbank", type=str, required=True)

parser.add_argument("-o", "--output", help="/path/to/output", type=str, required=True)

parser.add_argument(
    "-start", "--cas9_start", help="start coordinate of cas9", type=int, required=True
)

parser.add_argument(
    "-end", "--cas9_end", help="stop coordinate of cas9", type=int, required=True
)

parser.add_argument(
    "-w",
    "--window_size",
    help="return hits only within W distance from start/stop, default 1500",
    type=int,
    default=1500,
    required=False,
)

# parser.add_argument('-f',
#                    '--features',
#                    help="/path/to/featureTable",
#                    type = str,
#                   required=False)

args = parser.parse_args()

# Constants ------------------------------------------------------------------------------------------------
min_antirepeat_len = 12
max_antirepeat_gaps = 2
min_mfe = -30
min_hloops = 2
# max_dist_from_cas9 = 1500


# Functions ------------------------------------------------------------------------------------------------
def get_interregions_gbk(
    genbank_path, output, output_bed, intergene_length=50, min_orf_size=300
):
    """
    Records the start/stop coordinates of all annotated genes (CDS, tRNA, rRNA, ncRNA),
    extracts regions between them, including contig termini
    """
    # If the outputs exist, remove
    if Path(output).exists():
        Path(output).unlink()
    if Path(output_bed).exists():
        Path(output_bed).unlink()

    for seq_record in SeqIO.parse(genbank_path, "genbank"):
        # Check to make sure there is sequence data
        n_count = seq_record.seq.count("N")
        seq_len = len(seq_record.seq)
        assert (
            n_count < seq_len and seq_len > 0
        ), f"There doesn't seem to be any sequence data in {genbank_path}"

        # Get the name of the contig
        contig = seq_record.id.split("|")[0]

        cds_list = []
        intergenic_records = []

        # Loop over the genome file, get the gene features on each of the strands
        for feature in seq_record.features:
            if (
                feature.type == "tRNA"
                or feature.type == "ncRNA"
                or feature.type == "rRNA"
                or feature.type == "RNA"
                or feature.type == "repeat_region"
            ):
                mystart = feature.location.start.position
                myend = feature.location.end.position
                cds_list.append((mystart, myend))
            if feature.type == "CDS" or feature.type == "misc_feature":
                mystart = feature.location.start.position
                myend = feature.location.end.position
                if myend - mystart >= min_orf_size:
                    cds_list.append((mystart, myend))

        if len(cds_list) == 0:
            print(
                f"There are no features annotated as CDS/tRNA/ncRNA/rRNA in {seq_record.id}, skipping"
            )

        # Get the 5' end sequence
        if len(cds_list) > 0:
            first_feature_start = cds_list[0][0]
            if first_feature_start >= intergene_length:
                intergene_seq = seq_record.seq[0:first_feature_start]
                ign_num = "ign_" + "0"
                ign_start = 1
                ign_stop = first_feature_start
                intergenic_records.append(
                    SeqRecord(
                        intergene_seq,
                        id=f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                        description="",
                    )
                )

        # Get all intergenic loci
        if len(cds_list) > 1:
            for i, pospair in enumerate(cds_list[1:]):
                # Compare current start position to previous end position
                last_end = cds_list[i][1]
                this_start = pospair[0]
                if this_start - last_end >= intergene_length:
                    intergene_seq = seq_record.seq[last_end:this_start]
                    ign_num = "ign_" + str(i)
                    ign_start = last_end + 1
                    ign_stop = this_start

                    intergenic_records.append(
                        SeqRecord(
                            intergene_seq,
                            id=f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                            description="",
                        )
                    )

        # Get the 3' end sequence
        if len(cds_list) > 0:
            last_feature_end = cds_list.pop()[1]
            if len(seq_record.seq) - last_feature_end >= intergene_length:
                intergene_seq = seq_record.seq[last_feature_end:]
                ign_num = "ign_" + str(len(cds_list))
                ign_start = last_feature_end + 1
                ign_stop = len(seq_record.seq)
                intergenic_records.append(
                    SeqRecord(
                        intergene_seq,
                        id=f"{contig}|{ign_num}|{ign_start}|{ign_stop}",
                        description="",
                    )
                )
        # Write the outputs
        with open(output, "a") as f:
            SeqIO.write(intergenic_records, f, "fasta")

        with open(output_bed, "a") as f:
            for seq_record in intergenic_records:
                contig, ign, start, stop = seq_record.id.split("|")
                print(contig, ign, start, stop, sep=",", file=f)


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


def blasttab7_to_df(blasttab7):
    """
    Converts the blast tabular output with headers (outfmt -7) to a dataframe
    if no hits, returns an empty dataframe
    """
    # get the names of the columns
    # if there are no hits, the line starting with "# Fields:" won't exist
    # and names will stay false
    names = False
    with open(blasttab7) as infile:
        for line in infile:
            if line.startswith("# Fields:"):
                columns_raw = line.strip().split(" ", maxsplit=2)[2]

                # sanitize the names
                names = (
                    columns_raw.replace(" ", "")
                    .replace(".", "")
                    .replace("%", "perc")
                    .split(",")
                )

    if names:
        df = pd.read_csv(blasttab7, names=names, comment="#", sep="\t")
    else:
        df = pd.DataFrame()
    return df


def add_features_to_seqrecord(seq_record, df_featuretable):
    """
    This function accepts a nucleotide seq and a dataframe of features
    The dataframe must have the following columns:
        1. start (1-based)
        2. stop (1-based)

    With optional columns:
        3. strand (integer)
        4. feature type (str) (CDS or repeat_region)
        5. product (str)
        6. protein_id (str)
        7. rpt_type (str)
        8. rpt_unit_seq (str)

    returns a seqrecord object
    """

    seq_record.annotations["molecule_type"] = "DNA"
    seq_record.annotations = {"molecule_type": "DNA"}

    for index, row in df_featuretable.iterrows():
        # biopython is 0-based
        start0 = int(row.start) - 1

        # add the feature
        feat = SeqFeature(
            FeatureLocation(int(start0), int(row.stop)),
            # strand = int(row.strand),
            # type = row.feature_type,
            # qualifiers = {"product" : row.annotation,
            #              "protein_id" : row.protein_id}
        )
        seq_record.features.append(feat)

        # decorate features
        try:
            if row.feature_type:
                feat.type = row.feature_type

            if row.strand:
                feat.strand = int(row.strand)

            # add qualifiers
            ##ORFs
            if feat.type == "CDS":
                if row.annotation:
                    feat.qualifiers.update({"product": row.annotation})

                if row.protein_id:
                    feat.qualifiers.update({"protein_id": row.protein_id})

            # REPEATS
            if feat.type == "repeat_region":
                if row.repeat_type:
                    feat.qualifiers.update({"rpt_type": row.repeat_type})

                if row.rpt_unit_seq:
                    feat.qualifiers.update({"rpt_unit_seq": row.rpt_unit_seq})
        except AttributeError:
            pass

    return seq_record


def add_antirepeats_to_ign(df_blast, ign_fasta):
    """ """
    new_records = []
    # get a seqrecord of the intergenic seqs
    seq_records = SeqIO.parse(ign_fasta, "fasta")

    # iterate through the intergenic seqs
    for seq_record in seq_records:
        # get the header
        contig = seq_record.id

        # get the hits on the intergenic seq
        df_repeats = df_blast.query("queryaccver == @contig")

        if not df_repeats.empty:
            # add the repeat as an annotation
            seq_record2 = add_features_to_seqrecord(seq_record, df_repeats)

            new_records.append(seq_record2)

    return new_records


def get_seq_around_feature(seq_records, tail_len=200):
    """
    Accepts a Biopython SeqRecord object annotated with features
    Extracts the local sequence up/downstream from each feature
    returns a dataframe with the coordinates (1-based) and sequence
    """
    data = []
    for seq_record in seq_records:
        for feature in seq_record.features:
            feature_start = feature.location.start.position
            feature_stop = feature.location.end.position
            feature_len = feature_stop - feature_start
            feature_seq = seq_record.seq[feature_start:feature_stop]

            # get the upstream tail
            tail_start = feature_start - tail_len
            if tail_start < 0:
                tail_start = 0
            seq5 = seq_record.seq[tail_start:feature_start]

            # get the downstream tail
            tail_end = feature_stop + tail_len
            seq3 = seq_record.seq[feature_stop:tail_end]

            # if we hit the edge of a contig, tail_end may be larger than the actual end
            if len(seq3) < (tail_end - feature_stop):
                tail_end = feature_stop + len(seq3)

                # make sure math is right!
                assert (
                    len(seq3) == tail_end - feature_stop
                ), f"""
                Check this contig and coordinates:
                {seq_record.id} from {feature_stop} to {tail_end}
                """

            # make sure we didn't get too much seq
            assert len(seq5) <= (tail_len + feature_len) and len(seq3) <= (
                tail_len + feature_len
            ), f"""
            something went wrong.
            The feature is from {feature_start} to {feature_stop}
            The 5' extracted seq is from {tail_start} to {feature_start} and is {seq5}
            The 3' extract seqs is from {feature_stop} to {tail_end} and is {seq3}
            The length should not be longer than {tail_len+feature_len}
            But the length of seq5 is {len(seq5)} and seq3 is {len(seq3)}
            """

            # add the upstream
            data.append(
                {
                    "nucleotide_accession": seq_record.id,
                    "start5": tail_start + 1,
                    "stop5": feature_start,
                    "seq5": str(seq5),
                    "feature_start": feature_start + 1,
                    "feature_stop": feature_stop,
                    "feature_seq": str(feature_seq),
                    "start3": feature_stop + 1,
                    "stop3": tail_end,
                    "seq3": str(seq3),
                }
            )

    df = pd.DataFrame.from_records(data)

    return df


def stop_at_terminator(sequence, blackout=0):
    """
    Accepts a DNA string
    returns the DNA string clipped after the first terminator sequence (TTTT/UUUU)
    optionally excludes the first n based defined by blackout

    """

    if "TTTT" in sequence:
        terminator_idx = sequence.upper().find("TTTT", blackout)

        new_stop = terminator_idx + 4
        sequence_clipped = sequence[:new_stop]

    elif "UUUU" in sequence:
        terminator_idx = sequence.upper().find("UUUU", blackout)

        new_stop = terminator_idx + 4
        sequence_clipped = sequence[:new_stop]
    else:
        sequence_clipped = sequence
    return sequence_clipped


def run_rnafold(sequence):
    """
    quantifies MFE,secondary structure, and number of helixes from a sequence string
    sequence string can be DNA/RNA
    returns a tuple
    """
    data = []

    ss, MFE = RNA.fold(sequence)
    (bg,) = forgi.load_rna(ss)
    num_hloops = len(list(bg.hloop_iterator()))
    data.append(
        {
            "sequence": sequence,
            "ss": ss,
            "MFE": MFE,
            "num_hloops": num_hloops,
        }
    )

    df = pd.DataFrame.from_records(data)
    return ss, MFE, num_hloops


# -----------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"window size is {args.window_size}")

    # check the inputs
    gbk = Path(args.input)

    # get intergenic seqs
    ign = Path(f"{gbk.stem}_ign.fna")
    bed = Path(f"{gbk.stem}_bed.csv")
    get_interregions_gbk(gbk, ign, bed)

    # Blast CRISPR DRs versus intergenic seqs
    DR_fasta = Path(f"{gbk.stem}_DR.fasta")
    blast_output = Path(f"{gbk.stem}_blastn.tab")
    DR_seq = ""
    with DR_fasta.open("w") as o:
        seq_record = SeqIO.read(gbk, "genbank")
        for feature in seq_record.features:
            if feature.type == "repeat_region":
                if feature.qualifiers["rpt_type"] == ["CRISPR"]:
                    DR_seq = str(feature.qualifiers["rpt_unit_seq"][0])

                    print(">repeat", DR_seq, sep="\n", file=o)

    if len(DR_seq) < 1:
        print(f"Could not find any CRISPR repeats annotated in the input, exiting")
        sys.exit()

    subprocess.run(
        f'blastn -query {ign} -subject {DR_fasta} -word_size 4 -outfmt "7 std slen qseq sseq" -dust no -evalue 100 > {blast_output}',
        check=True,
        shell=True,
    )
    df_blast = blasttab7_to_df(blast_output)

    assert (
        len(df_blast) > 0
    ), f"""
    No alignments were generated between the CRISPR direct repeat and intergenic seqs
    Are you sure both are present? 
    You may want to adjust the minimum size of intergenic sequences, if not
    """

    # filter the hits
    cond1 = df_blast.alignmentlength >= min_antirepeat_len
    cond2 = df_blast.alignmentlength < df_blast.subjectlength * 0.9
    cond3 = df_blast.gapopens <= max_antirepeat_gaps
    df_blast2 = df_blast.query("@cond1 and @cond2 and @cond3").rename(
        columns={
            "subjectseq": "crispr_dr",
            "queryseq": "anti-repeat",
        }
    )

    # Add the candidate anti-repeats in the query intergenic sequences as annotations
    df_blast3 = df_blast2.loc[:, ["queryaccver", "qstart", "qend", "crispr_dr"]].rename(
        columns={"qstart": "start", "qend": "stop"}
    )
    df_blast3["strand"] = 1
    df_blast3["feature_type"] = "repeat_region"
    df_blast3["repeat_type"] = "anti-repeat"
    seq_records_with_antirepeats = add_antirepeats_to_ign(df_blast3, ign)

    if len(seq_records_with_antirepeats) < 1:
        print("no antirepeats, exiting")
        sys.exit()

    # extract the adjacent sequences and construct candidate tracrs
    df_window = get_seq_around_feature(seq_records_with_antirepeats)
    df_window = df_window.rename(
        columns={
            "feature_seq": "anti-repeat",
            "feature_start": "anti-repeat_start",
            "feature_stop": "anti-repeat_stop",
        }
    ).drop(
        columns=[
            "start5",
            "stop5",
            "start3",
            "stop3",
        ]
    )

    df_window["tracr5"] = df_window["seq5"] + df_window["anti-repeat"]
    df_window["tracr3"] = df_window["anti-repeat"] + df_window["seq3"]

    # Add the CRISPR repeat hit to the dataframe
    df_tracr = pd.merge(
        df_window,
        df_blast3,
        how="left",
        left_on=["nucleotide_accession", "anti-repeat_start", "anti-repeat_stop"],
        right_on=["queryaccver", "start", "stop"],
    ).drop(
        columns=[
            "queryaccver",
            "start",
            "stop",
            "strand",
            "feature_type",
            "repeat_type",
        ]
    )

    df_tracr["crispr_dr"] = df_tracr.crispr_dr.str.replace(
        "-", ""
    )  # remove blast alignment gaps from the CRISPR DR

    data = []
    for index, row in df_tracr.iterrows():
        crispr_dr = row.crispr_dr
        crispr_dr_rc = str(Seq(crispr_dr).reverse_complement())

        ####evaluate upstream of the anti-repeat###
        candidate_tracr5 = row.tracr5
        candidate_tracr5rc = str(Seq(candidate_tracr5).reverse_complement())
        candidate_tracr5_clipped = stop_at_terminator(candidate_tracr5rc, 40)

        # make a hybrid with the CRISPR DR
        hybrid = crispr_dr + "&" + candidate_tracr5_clipped

        # evaluate with RNA Fold
        ss, MFE, num_hloops = run_rnafold(hybrid)

        # check non-canonical
        hybrid2 = crispr_dr_rc + "&" + candidate_tracr5_clipped
        ss2, MFE2, num_hloops2 = run_rnafold(hybrid2)

        if MFE2 < MFE:
            data.append(
                {
                    "nucleotide_accession": row.nucleotide_accession,
                    "tracr": candidate_tracr5_clipped,
                    "hybrid": hybrid2,
                    "ss": ss2,
                    "MFE": MFE2,
                    "num_hloops": num_hloops2,
                }
            )
        else:
            data.append(
                {
                    "nucleotide_accession": row.nucleotide_accession,
                    "tracr": candidate_tracr5_clipped,
                    "hybrid": hybrid,
                    "ss": ss,
                    "MFE": MFE,
                    "num_hloops": num_hloops,
                }
            )

        ###evaluate downstream of the anti-repeat###
        candidate_tracr3 = row.tracr3
        candidate_tracr3_clipped = stop_at_terminator(candidate_tracr3, 40)

        # make a hybrid with the reverse complement of the CRISPR DR
        hybrid = crispr_dr_rc + "&" + candidate_tracr3_clipped

        # evaluate with RNA Fold
        ss, MFE, num_hloops = run_rnafold(hybrid)

        # check non-canonical
        hybrid2 = crispr_dr + "&" + candidate_tracr3_clipped
        ss2, MFE2, num_hloops2 = run_rnafold(hybrid2)

        if MFE2 < MFE:
            data.append(
                {
                    "nucleotide_accession": row.nucleotide_accession,
                    "tracr": candidate_tracr3_clipped,
                    "hybrid": hybrid2,
                    "ss": ss2,
                    "MFE": MFE2,
                    "num_hloops": num_hloops2,
                }
            )
        else:
            data.append(
                {
                    "nucleotide_accession": row.nucleotide_accession,
                    "tracr": candidate_tracr3_clipped,
                    "hybrid": hybrid,
                    "ss": ss,
                    "MFE": MFE,
                    "num_hloops": num_hloops,
                }
            )

    df_rnafold = pd.DataFrame.from_records(data)

    # get the nt coordinates of the tracrs on the parent molecule
    data = []
    for index, row in df_rnafold.iterrows():
        contig = row.nucleotide_accession
        tracr = row.tracr
        seq_record = SeqIO.read(gbk, "genbank")

        start, stop, strand = get_subseq_coords(tracr, seq_record)
        data.append({"tracr": tracr, "start": start, "stop": stop, "strand": strand})

    df_coords = pd.DataFrame.from_records(data)
    df_rnafold2 = pd.merge(df_rnafold, df_coords, how="left", on="tracr")

    # add cas9 coordinates
    df_rnafold2["cas9_start"] = args.cas9_start
    df_rnafold2["cas9_end"] = args.cas9_end
    df_rnafold2["dist1"] = (df_rnafold2["start"] - df_rnafold2["cas9_start"]).abs()
    df_rnafold2["dist2"] = (df_rnafold2["start"] - df_rnafold2["cas9_end"]).abs()
    df_rnafold2["dist3"] = (df_rnafold2["stop"] - df_rnafold2["cas9_start"]).abs()
    df_rnafold2["dist4"] = (df_rnafold2["stop"] - df_rnafold2["cas9_end"]).abs()
    df_rnafold2["dist_from_cas9"] = df_rnafold2[
        ["dist1", "dist2", "dist3", "dist4"]
    ].min(axis=1)
    df_rnafold2 = df_rnafold2.drop(
        columns=["cas9_start", "cas9_end", "dist1", "dist2", "dist3", "dist4"]
    )

    # filter for likely hits
    cond1 = df_rnafold2.MFE <= min_mfe
    cond2 = df_rnafold2.num_hloops >= min_hloops
    cond3 = df_rnafold2.dist_from_cas9 <= args.window_size
    df_rnafold3 = df_rnafold2.query("@cond1 and @cond2 and @cond3").sort_values(
        "dist_from_cas9"
    )

    df_rnafold3.to_csv(args.output, index=False)

    # cleanup
    DR_fasta.unlink()
    blast_output.unlink()
    ign.unlink()
    bed.unlink()
