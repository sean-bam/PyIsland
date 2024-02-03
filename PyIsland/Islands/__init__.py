from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from PyIsland import Parsers
import pandas as pd


def make_gbk(fasta, gff, gbk, contig, start=1, stop=1000000000):
    """
    Accepts a fasta file, corresponding GFF from prodigal
    and a contig accession, with optional 1-based start/stop coordinates

    Adds the ORFs as features and prints a genbank-formatted file
    """

    contig = str(contig)

    num_records_output = 0
    # get a seq record of the contig
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if seq_record.id == contig:
            # get a dataframe of ORF start/stops on the contig of interest
            # df = prodigal_fasta2df(gff).query('contig == @contig')
            df = prodigal_gff_to_df(gff).query("contig == @contig")

            if df.empty:
                print(f"no ORFs found on {contig}")
            else:
                df["feature_type"] = "CDS"

                # add the ORFs to the seq record
                seq_record = add_features_to_seqrecord(seq_record, df)

            # window
            start = start - 1  # adjust to 0 based
            assert start >= 0, f"start location must be 1 or greater"
            if stop > len(seq_record):
                stop = len(seq_record)

            seq_record = seq_record[start:stop]

            # update the ID/source to reflect the subsequence

            start = start + 1  # adjust back to 1-based

            # this causes errors with reverse complementing
            # and is inconsistent with the NCBI
            # df_source = pd.DataFrame.from_records([{"contig" : contig,
            #                                       "feature_type": "source",
            #                                       "start" : start,
            #                                       "stop" : stop,
            #                                       "strand" : 1 #arbitrary
            #                                      }]
            #                                     )
            # seq_record = add_features_to_seqrecord(seq_record, df_source)

            # updates, but SeqIO.write secretly sanitizes this field
            # seq_record.annotations["accessions"] = [seq_record.name, "REGION:", f"{start}..{stop}"]

            # this actually gets written, must be a single line
            seq_record.annotations["source"] = [f"{start}..{stop}"]

            seq_record.id = seq_record.id + f":{start}-{stop}"
            seq_record.description = ""

            # write Gbk file
            num_records_output = SeqIO.write(seq_record, gbk, "genbank")

            assert (
                num_records_output == 1
            ), f"""
            Wrote more than 1 record to {gbk}, with {contig} given as an ID
            """
            break

    assert num_records_output > 0, f"Could not find {contig} in {fasta}"


def add_features_to_seqrecord(seq_record, df_featuretable):
    """
    This function accepts a nucleotide seq and a dataframe of features
    The dataframe must have the following columns:
        0. contig
        1. start (1-based)
        2. stop (1-based)
        3. strand (integer)

    With optional columns:
        4. feature type (str)
        5. product (str, if feature type == CDS)
        6. protein_id (str, if feature type == CDS)
        7. rpt_type (str)
        8. rpt_unit_seq (str)

    returns a seqrecord object
    """

    seq_record.annotations["molecule_type"] = "DNA"
    seq_record.annotations = {"molecule_type": "DNA"}

    # check dataframe and seq record match
    contigs = df_featuretable.contig.unique().tolist()
    assert (
        len(contigs) == 1
    ), f""" 
    You provided multiple contig IDs, which isn't supported yet
    """

    assert (
        contigs[0] == seq_record.id
    ), f"""
    The contig ID in the feature table doesnt match the 
    seq record ID
    """

    columns = df_featuretable.columns
    assert (
        "start" in columns and "stop" in columns and "strand" in columns
    ), f"""
    Did not find "start", "stop", and/or "strand" in the input dataframe
    """

    # if start > stop, update to start < stop
    cond1 = df_featuretable.start > df_featuretable.stop
    df_featuretable[["start", "stop"]] = df_featuretable.loc[:, ["start", "stop"]].mask(
        cond1, df_featuretable[["stop", "start"]], axis="index"
    )

    for index, row in df_featuretable.iterrows():
        # biopython is 0-based
        start0 = int(row.start) - 1

        # biopython requires strand be integers
        strand = None
        if row.strand == "+":
            strand = 1
        elif row.strand == "-":
            strand = -1
        else:
            strand = row.strand
        assert type(strand) == int, f"couldnt convert strand {strand} to an integer"

        # add the feature
        feat = SeqFeature(
            FeatureLocation(int(start0), int(row.stop)),
            strand=int(strand),
        )

        seq_record.features.append(feat)

        # decorate features
        if "feature_type" in columns:
            feat.type = row.feature_type

            # add qualifiers
            ##ORFs
            if feat.type == "CDS":
                if "annotation" in columns:
                    feat.qualifiers.update({"product": row.annotation})

                if "protein_id" in columns:
                    feat.qualifiers.update({"protein_id": row.protein_id})

            # REPEATS
            if feat.type == "repeat_region":
                if "repeat_type" in columns:
                    feat.qualifiers.update({"rpt_type": row.repeat_type})

                if "rpt_unit_seq" in columns:
                    feat.qualifiers.update({"rpt_unit_seq": row.rpt_unit_seq})

    return seq_record


class genomicWindows:
    def add_dist_from_seed(df):
        """
        return the minimum number of rows separating a given
        row from a seed row
        """
        df = df.sort_values("start").reset_index(drop=True)
        all_indexes = df.index
        seed_indexes = df.query("seed == True").index
        df["dist_from_seed"] = genomicWindows.get_minimum(all_indexes, seed_indexes)

        return df

    def get_minimum(seriesA, seriesB):
        """
        return the minimum distance for each element in series A from all elements in series B
        """
        return [min([abs(a - b) for b in seriesB]) for a in seriesA]

    def select_neighborhood(gff: str, seedgff: str, window: int, check=True):
        """
        Accepts two GFF files and a window size in bp
        Extracts +/- window size around seeds from main gff

        Optionally checks all contigs given in seed gff
        are also present in main gff

        returns a dataframe
        """

        df_gff = Parsers.gff_to_df(gff)
        df_seed = Parsers.gff_to_df(seedgff)
        df_seed["seed"] = True

        # get a subdataframes of corresponding seed/GFF contigs
        df_seed = df_seed.query("contig in @df_gff.contig")
        df_gff = df_gff.query("contig in @df_seed.contig")

        if df_seed.empty or df_gff.empty:
            if check:
                raise AssertionError(f"{gff} and {seedgff} share no contigs in common")
            else:
                print(f"{gff} and {seedgff} share no contigs in common, exiting")
                sys.exit()

        # combine dataframes
        df = pd.concat([df_gff, df_seed]).sort_values(["contig", "start"])

        # if seeds are identical to a main GFF feature
        # drop the duplicate main GFF feature
        # so that distance is calculated correctly
        df = df.sort_values("seed").drop_duplicates(
            subset=[
                "contig",
                "feature_type",
                "start",
                "stop",
                "attributes",
            ]
        )

        df_list = []
        for index, row in df_seed.iterrows():
            contig = row.contig
            start = row.start - window
            stop = row.stop + window

            # select ORFs that start before the stop limit
            # or stop after the start limit
            ##inclusive##
            df_neighborhood = df.query(
                "contig == @contig and start <= @stop and stop >= @start"
            )
            
            #label island
            isl_start = str(df_neighborhood.start.min())
            isl_stop = str(df_neighborhood.stop.max())
            df_neighborhood["island"] = str(contig) + ":" + isl_start + "-" + isl_stop

            # label distance from seeds
            df_neighborhood2 = genomicWindows.add_dist_from_seed(df_neighborhood)

            df_list.append(df_neighborhood2)

        df2 = (
            pd.concat(df_list)
            .reset_index(drop=True)
            .drop(columns=["seed"])
            .sort_values(
                [
                    "contig",  # in case of ORFs near multiple seeds
                    "start",
                    "dist_from_seed",
                ]
            )
            .drop_duplicates(
                subset=[
                    "contig",  # take the one labelled as closest to a seed
                    "feature_type",
                    "start",
                    "stop",
                    "attributes",
                ]
            )
        )

        return df2
