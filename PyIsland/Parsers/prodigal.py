import pandas as pd

from pyIsland import Parsers


def prodigal_fasta2df(prodigal_fasta):
    """
    This function takes the fasta output from prodigal and returns a pandas dataframe
    """
    data = []
    with open(prodigal_fasta) as fastafile:
        for line in fastafile:
            if line.startswith(">"):
                newline = line.strip(">")
                protein_id, start, stop, strand, info = newline.split(" # ")
                contig, orf = protein_id.split("_")
                # if int(strand) == 1:
                #    sign = '+'
                # else:
                #    sign = '-'
                data.append([contig, orf, protein_id, start, stop, strand, info])
        df = pd.DataFrame(
            data,
            columns=["contig", "orf", "protein_id", "start", "stop", "strand", "info"],
        )

        return df


def prodigal_gff2df(prodigal_gff):
    """
    DEPRECATED
    This function takes the GFF output from prodigal and returns a datafame
    """
    data = []
    with open(prodigal_gff) as f:
        for line in f:
            if not line.startswith("#"):
                (
                    contig,
                    prodigal_version,
                    CDS,
                    start,
                    stop,
                    number0,
                    strand,
                    number,
                    info,
                ) = line.split()
                (
                    id_raw,
                    partial_raw,
                    start_type_raw,
                    rbs_raw,
                    rbs2_raw,
                    gc_raw,
                    conf_raw,
                    score_raw,
                    cscore_raw,
                    sscore_raw,
                    rscore_raw,
                    uscore_raw,
                    tscore_raw,
                    empty,
                ) = info.split(";")

                orf_id = contig + "_" + id_raw.split("_")[1]
                partial = partial_raw.split("=")[1]
                confidence = round(float(conf_raw.split("=")[1]))
                gc = round(float(gc_raw.split("=")[1]), 2)

                if strand == "+":
                    strand = 1
                else:
                    strand = -1

                data.append(
                    {
                        "contig": str(contig),
                        # "prodigal_version" : prodigal_version,
                        "start": int(start),
                        "stop": int(stop),
                        "strand": strand,
                        "protein_id": str(orf_id),
                        "partial": str(partial),
                        "confidence": confidence,
                        "gc": gc,
                    }
                )

        return pd.DataFrame.from_records(data)


def prodigal_gff_to_df(gff):
    df = Parsers.gff_to_df(gff)

    df["id"] = df.attributes.str.split(";", expand=True)[0]
    df["protein_id"] = df.id.str.split("=", expand=True)[1]

    # df["strand"] = df.mask(df.strand == "+", 1).strand
    # df["strand"] = df.mask(df.strand == "-", -1).strand

    return df.drop(columns=["id", "attributes", "score", "phase", "source"])
