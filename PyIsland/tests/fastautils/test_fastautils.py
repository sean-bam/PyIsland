from PyIsland import FastaUtils
from pathlib import Path
from PyIsland import tests
import pytest

@pytest.fixture
def get_test_fasta():
    return Path('PyIsland/tests/example.fasta')

def test_subseq_fasta_biopython(get_test_fasta):
    seqs_to_get = ['seq1']
    
    outfasta = Path('PyIsland/tests/fastautils/subseq_fasta_biopython.fasta')
    expected_output = Path('PyIsland/tests/fastautils/subseq_fasta_expected_output.fasta')
    
    FastaUtils.subseq_fasta_biopython(seqs_to_get,
                                      get_test_fasta,
                                      outfasta)
    
    output_hash = tests.md5sum(outfasta)
    expected_output_hash = tests.md5sum(expected_output)
    
    assert output_hash == expected_output_hash
    outfasta.unlink()
    
def test_slice_fasta_biopython(get_test_fasta):
    seqs_to_get = 'seq2'
    start = 4
    stop = 9
    outfasta = Path('PyIsland/tests/fastautils/slice_fasta_biopython.fasta')
    expected_output = Path('PyIsland/tests/fastautils/slice_fasta_expected_output.fasta')
    
    FastaUtils.slice_fasta_biopython(seqs_to_get,
                                     start,
                                     stop,
                                     get_test_fasta,
                                     outfasta
                                    )
    
    output_hash = tests.md5sum(outfasta)
    expected_output_hash = tests.md5sum(expected_output)
    
    assert output_hash == expected_output_hash
    outfasta.unlink()
    
def test_fa_strict(get_test_fasta):
    
    outfasta = Path('PyIsland/tests/fastautils/fa_strict.fasta')
    expected_output = Path('PyIsland/tests/fastautils/fa_strict_expected_output.fasta')
    
    FastaUtils.fa_strict(get_test_fasta,
                         outfasta,
                         min_len=5
                        )
    
    output_hash = tests.md5sum(outfasta)
    expected_output_hash = tests.md5sum(expected_output)
    
    assert output_hash == expected_output_hash
    outfasta.unlink()
    
def test_count_seqs_in_fasta(get_test_fasta):
    
    num_seqs = FastaUtils.count_seqs_in_fasta(get_test_fasta)
    
    assert num_seqs == 6
    
def test_deduplicate_fasta(get_test_fasta):
    
    outfasta = Path('PyIsland/tests/fastautils/dedup.fasta')
    expected_output = Path('PyIsland/tests/fastautils/dedup_expected_output.fasta')
    
    FastaUtils.deduplicate_fasta(get_test_fasta,
                                 outfasta,
                                )
    
    output_hash = tests.md5sum(outfasta)
    expected_output_hash = tests.md5sum(expected_output)
    
    assert output_hash == expected_output_hash
    outfasta.unlink()