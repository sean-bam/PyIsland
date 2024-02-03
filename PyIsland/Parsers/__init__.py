from . import infernal
from . import mmseqs
from . import minced
from . import prodigal
from . import hhsuite
from . import blast

import pandas as pd
import numpy as np
from Bio import SeqIO


def gff_to_df(gff):
    df = pd.read_csv(
        gff,
        sep="\t",
        comment="#",
        names=[
            "contig",
            "source",
            "feature_type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
        dtype={
            "contig": "str",
            "source": "str",
            "feature_type": "str",
            "start": "int64",
            "stop": "int64",
            "score": "str",  # b/c sometimes is "."
            "strand": "str",
            "phase": "str",
            "attributes": "str",
        },
    )

    return df


def genbank_to_df(gbk):
    data = []
    with open(gbk) as f:
        seq_records = SeqIO.parse(f, "genbank")
        for seq_record in seq_records:
            for feature in seq_record.features:
                protein_id = np.nan
                product = np.nan
                if feature.type == "CDS":
                    try:
                        protein_id = feature.qualifiers.get("protein_id")[0]
                        product = feature.qualifiers.get("product")[0]
                    except TypeError:
                        pass

                data.append(
                    {
                        "contig": seq_record.id,
                        "start": feature.location.start,
                        "stop": feature.location.end,
                        "strand": feature.strand,
                        "feature_type": feature.type,
                        "protein_id": protein_id,
                        "product": product,
                    }
                )

    return pd.DataFrame.from_records(data)
