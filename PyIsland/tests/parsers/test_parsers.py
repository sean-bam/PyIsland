from PyIsland import Parsers
from pathlib import Path
from PyIsland import tests
import pytest


@pytest.fixture
def get_test_gff():
    return Path("PyIsland/tests/lambda.gff")


@pytest.fixture
def get_test_gbk():
    return Path("PyIsland/tests/lambda.gbk")


@pytest.fixture
def get_test_fasta():
    return Path("PyIsland/tests/lambda.fna")


@pytest.fixture
def get_test_orfs():
    return Path("PyIsland/tests/lambda.faa")


def test_gff_to_df(get_test_gff):
    df = Parsers.gff_to_df(get_test_gff)

    # check I can get a contig using the .query method and a string
    assert not df.query('contig == "NC_001416.1"').empty


# def test_gff_to_df_with_prodigal_gff():

#    df = Parsers.gff_to_df('PyIsland/tests/example_prodigal.gff')

# check I can get a contig using the .query method and a string
#    assert not df.query('contig == "10007077.1"').empty


def test_gbk_to_df(get_test_gbk):
    df = Parsers.genbank_to_df(get_test_gbk)

    features = df.feature_type.unique().tolist()
    expected_features = [
        "source",
        "gene",
        "CDS",
        "sig_peptide",
        "mat_peptide",
        "mRNA",
        "variation",
        "misc_recomb",
        "misc_binding",
        "regulatory",
        "misc_feature",
        "unsure",
    ]

    # check I get all features
    assert features == expected_features

    # check I labelled protein_ids
    assert df.protein_id.notna().sum() == 73

    # check products are labelled
    assert len(df.query("product.notna()")) == 73


def _test_gff_to_df(get_test_gff):
    seqs_to_get = ["seq1"]

    outfasta = Path("PyIsland/tests/fastautils/subseq_fasta_biopython.fasta")

    df = Parsers.gff_to_df(get_test_gff)

    output_hash = tests.md5sum(outfasta)
    expected_output_hash = tests.md5sum(expected_output)

    assert output_hash == expected_output_hash
    outfasta.unlink()


def test_prodigal_gff_to_df(get_test_gff):
    df = Parsers.prodigal.prodigal_gff_to_df("PyIsland/tests/example_prodigal.gff")

    # check I can get a contig using the .query method and a string
    assert not df.query('contig == "10007077.1"').empty

    # check the protein_id was parsed
    assert df.protein_id.isna().sum() == 0
