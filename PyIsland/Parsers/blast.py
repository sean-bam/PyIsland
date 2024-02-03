import pandas as pd


def blasttab7_to_df(blasttab7):
    """
    Converts the blast tabular output with headers (outfmt 7) to a dataframe
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


def blast_outfmt6_to_df(file):
    df = pd.read_csv(
        file,
        sep="\t",
        names=[
            "qaccver",
            "saccver",
            "pid",
            "length",
            "mismatch",
            "gap",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bit",
        ],
    )
    return df


class overlap_filtering(object):
    def add_current_hit_to_df(df, row):
        """
        Appends the current hit info as columns
        to a blast dataframe

        Probably a better way to do this,
        but works for now
        """
        df["qaccver2"] = row.qaccver
        df["saccver2"] = row.saccver
        # df["length2"] = row.length
        df["qstart2"] = row.qstart
        df["qend2"] = row.qend
        df["bitscore2"] = row.bitscore

        return df

    def add_overlap(df):
        """
        Adds columns calculating the overlap of current hit
        against all others
        """

        # add length of alignment
        df["q_aln"] = df["qend"] - df["qstart"]
        df["q_aln2"] = df["qend2"] - df["qstart2"]

        # calc intersection
        df["intersection"] = df[["qend", "qend2"]].min(axis=1) - df[
            ["qstart", "qstart2"]
        ].max(axis=1)

        # calc min length
        df["min_len"] = df[["q_aln", "q_aln2"]].min(axis=1)

        # calc overlap
        df["overlap"] = df["intersection"] / df["min_len"]

        return df

    def get_indexes_of_overlapping_hits(df, threshold):
        """
        selects rows where overlap exceeds threshold and evalue is below current hit

        returns a
        """

        return df.query(
            "overlap > @threshold and bitscore < bitscore2 and saccver != saccver2"
        ).index

    def remove_overlapping_hits(df, threshold):
        # initialize empty list of indexes
        idx_to_remove = df.iloc[0:0].index

        # columns to keep
        columns = df.columns.tolist()

        for index, row in df.iterrows():
            df2 = overlap_filtering.add_current_hit_to_df(df, row)
            df3 = overlap_filtering.add_overlap(df2)
            overlapping_idxs = overlap_filtering.get_indexes_of_overlapping_hits(
                df3, threshold
            )
            idx_to_remove = idx_to_remove.append(overlapping_idxs)

        # remove overlapping hits
        df4 = df.query("index not in @idx_to_remove")

        # remove extra columns
        df5 = df4.loc[:, columns]

        return df5

    def select_nonoverlapping_hits(df, threshold=0.33):
        df_list = []
        
        #qstart must be > qend
        #otherwise, overlap can be falsely calculated
        df[["qstart","qend"]] = (df.loc[:,["qstart","qend"]]
                                   .mask(df.qstart > df.qend, df[["qend","qstart"]], axis = 'index')
                                )

        for name, group in df.groupby("qaccver"):
            df2 = overlap_filtering.remove_overlapping_hits(group, threshold)
            df_list.append(df2)

        return pd.concat(df_list)
