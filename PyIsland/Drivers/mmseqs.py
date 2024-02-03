from pathlib import Path
import subprocess
import shutil
import tempfile

from pyIsland import Parsers


def run_mmseqs_cls(
    fasta,
    seqdb,
    clsdb,
    tsv,
    cluster_mode=2,
    coverage_mode=1,
    seq_id=0.5,
    cov=0.75,
    linclust=False,
    overwrite=False,
):
    cluster_algorithm = "cluster"
    params = f"--cluster-mode {cluster_mode} --cov-mode {coverage_mode} --min-seq-id {seq_id} -c {cov} "

    if linclust:
        cluster_algorithm = "linclust"

    # seqdb = Path(f'{dbname}seqDB')
    # clsdb = Path(f'{dbname}clsDB')
    tmp_dir = tempfile.TemporaryDirectory(dir='.')
    tmp_dir_path = Path(tmp_dir.name)
    tsv = Path(tsv)

    if overwrite:
        try:
            seqdb.with_suffix(".dbtype").unlink()
            clsdb.with_suffix(".dbtype").unlink()
        except FileNotFoundError:
            pass

    else:
        assert not clsdb.with_suffix(
            ".dbtype"
        ).is_file(), f"""
        ClusteringDB exists already, delete it or set overwrite = True
        """

    subprocess.run(f"mmseqs createdb {fasta} {seqdb}", shell=True, check=True)

    subprocess.run(
        f"mmseqs {cluster_algorithm} {seqdb} {clsdb} {tmp_dir_path} {params}",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
    )

    subprocess.run(
        f"mmseqs createtsv {seqdb} {seqdb} {clsdb} {tsv}", shell=True, check=True
    )

    # for file in Path('tmpmmseqsdb').rglob('*'):
    #    try:
    #        file.unlink()
    #    except IsADirectoryError:
    #        pass
    #shutil.rmtree("tmpmmseqsdb/")
    
def get_reps_from_clustering(seqDB, clsDB, fasta):
    
    tmp_dir = tempfile.TemporaryDirectory(dir='.')
    tmp_dir_path = Path(tmp_dir.name)
    repdb = tmp_dir_path / Path('repDB')
    
    subprocess.run(
        f"mmseqs createsubdb {clsDB} {seqDB} {repdb}",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    
    subprocess.run(
        f"mmseqs convert2fasta {repdb} {fasta}",
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
    )



def cls2pro(
    seqDB, clsDB, proDB, evalue=0.001, qid=0.0, qsc=-20, diff=1000, max_seq_id=0.9
):
    filtering_params = (
        f"--qid {qid} --qsc {qsc} --diff {diff} --max-seq-id {max_seq_id}"
    )

    subprocess.run(
        f"mmseqs createsubdb {clsDB} {seqDB} {clsDB}rep", shell=True, check=True, stdout=subprocess.DEVNULL
    )

    subprocess.run(
        f"mmseqs createsubdb {clsDB} {seqDB}_h {clsDB}rep_h", shell=True, check=True, stdout=subprocess.DEVNULL
    )

    subprocess.run(
        f"mmseqs result2profile {clsDB}rep {seqDB} {clsDB} {proDB} -e {evalue} --e-profile {evalue} {filtering_params}",
        shell=True,
        check=True,
    )


def proVScls(seqDB, clsDB, proDB, alnDB, tsv, evalue=10):
    subprocess.run(
        f"mmseqs align {proDB} {seqDB} {clsDB} {alnDB} -e {evalue} --max-accept 100000 --alt-ali 10 -a",
        shell=True,
        check=True,
    )

    subprocess.run(
        f"mmseqs convertalis {proDB} {seqDB} {alnDB} {tsv} --format-mode 4 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov",
        shell=True,
        check=True,
    )


def filter_proVScls(
    seqDB, proDB, alnDB, alnDBfilt, tsv, qid=0, qsc=-20, diff=0, max_seq_id=1.0, cov=0
):
    filtering_params = (
        f"--qid {qid} --qsc {qsc} --diff {diff} --max-seq-id {max_seq_id} --cov {cov}"
    )

    # subprocess.run(f'mmseqs convertalis {proDB} {seqDB} {alnDB} pro_vs_member.tsv',
    #               shell = True,
    #               check = True)

    subprocess.run(
        f"mmseqs filterresult {proDB} {seqDB} {alnDB} {alnDBfilt} {filtering_params}",
        shell=True,
        check=True,
    )

    subprocess.run(
        f"mmseqs convertalis {proDB} {seqDB} {alnDBfilt} {tsv} --format-mode 4 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov",
        shell=True,
        check=True,
    )


def filter_clusterdb(
    seqDB,
    clsDB,
    alnDB,
    seqDB_filtered,
    clsDB_filtered,
    qid=0,
    qsc=0,
    diff=0,
    max_seq_id=0.9,
):
    """
    Runs MMSeqs filterresult on a clustering database, and outputs a new
    clustering+sequence database

    Note: filterresult automatically removes singleton clusters, even with
    all filtering flags disabled

    Uses this hack to convert a resultDB into a clusterDB
    https://github.com/soedinglab/MMseqs2/issues/316
    """
    clsdb_path = Path(clsDB)
    clsdb_filtered_path = Path(clsDB_filtered)

    cls_dbtype = clsdb_path.with_suffix(".dbtype")
    assert cls_dbtype.is_file(), f"cant find {cls_dbtype}"

    # run the filtering
    filtering_params = (
        f"--qid {qid} --qsc {qsc} --diff {diff} --max-seq-id {max_seq_id}"
    )

    subprocess.run(
        f"mmseqs filterresult {seqDB} {seqDB} {clsDB} {alnDB} {filtering_params}",
        shell=True,
        check=True,
    )

    # trick mmseqs into thinking the resultdb is a clusterdb
    subprocess.run(
        f"mmseqs filterdb {alnDB} {clsDB_filtered} --trim-to-one-column",
        shell=True,
        check=True,
    )

    res_dbtype = clsdb_filtered_path.with_suffix(".dbtype")

    with cls_dbtype.open() as fsource:
        with res_dbtype.open("w") as fdest:
            shutil.copyfileobj(
                fsource, fdest
            )  # overwrite the resultdb type with the original clustering database type

    # repopulate a seqDB from the filtered DB
    ##first, copy all the mmseqs internal IDs from the filtered DB to a file
    tmpfile = clsdb_filtered_path.parent / Path("tmpids.list")

    with tmpfile.open("w") as o:
        for file in clsdb_filtered_path.parent.glob(clsdb_filtered_path.stem + ".*"):
            # loops through all clsDB.[0-9]* files.
            # also ensuring filenames match.
            if file.suffix != ".dbtype" and file.suffix != ".index":
                print(file)
                with file.open() as infile:
                    for line in infile:
                        internal_id = line.strip().replace("\x00", "")
                        if len(internal_id) > 0:  # empty IDs mess things up
                            print(line.strip().replace("\x00", ""), file=o)

    subprocess.run(
        f"mmseqs createsubdb {tmpfile} {seqDB} {seqDB_filtered}", shell=True, check=True
    )


def get_filtered_seqs_as_db(seqDB, unfiltered_alnDB, filtered_alnDB, new_seqDB):
    subprocess.run(
        f"cut -f1 {filtered_alnDB}.* | sort | uniq > passed.ids", shell=True, check=True
    )

    subprocess.run(
        f"cut -f1 {unfiltered_alnDB}.* | sort | uniq > all.ids", shell=True, check=True
    )

    subprocess.run(f"comm -23 all.ids passed.ids > missed.ids", shell=True, check=True)

    subprocess.run(
        f"mmseqs createsubdb missed.ids {seqDB} {new_seqDB} --subdb-mode 1",
        shell=True,
        check=True,
    )

    Path("all.ids").unlink()
    Path("passed.ids").unlink()
    Path("missed.ids").unlink()


def get_missing_seqs_as_fasta_using_indexes(db1, db2, output):
    """
    Accepts paths to MMSeqs DB1 and DB2
    Gets their corresponding .index files
    Returns a fasta file of sequences in DB1 not in DB2
    """

    # set paths
    db1_idx = Path(db1).with_suffix(".index")
    db2_idx = Path(db2).with_suffix(".index")
    tmp_dir = tempfile.TemporaryDirectory(dir=".")
    tmp_dir_path = Path(tmp_dir.name)
    missing = tmp_dir_path / Path("missing.ids")
    subDB = tmp_dir_path / Path("subDB")

    # check index files are present
    assert (
        db1_idx.is_file() and db2_idx.is_file()
    ), f"Can't find {db1_idx} and/or {db2_idx}"

    # get indexes in DB1 missing from DB2
    db1_ids = Parsers.mmseqs.get_mmseqs_internal_ids_as_set(db1_idx)
    db2_ids = Parsers.mmseqs.get_mmseqs_internal_ids_as_set(db2_idx)
    missing_ids = db1_ids - db2_ids
    assert len(missing_ids) > 0, f"No seqs in {db1} are absent from {db2}"

    # run mmseqs to get them as fasta
    with missing.open("w") as o:
        for internal_id in missing_ids:
            print(internal_id, file=o)

        subprocess.run(
            f"mmseqs createsubdb {missing} {db1} {subDB}", shell=True, check=True
        )

        subprocess.run(f"mmseqs convert2fasta {subDB} {output}", shell=True, check=True)


def get_missing_seqs_as_fasta_using_tsv(db1, tsv, output):
    """
    Accepts paths to MMSeqs aminoacid DB and a two-column TSV file (e.g., clustering)
    Gets the corresponding .lookup file from db1
    Returns a fasta file of sequences in DB1 not listed in the TSV
    """

    # set paths
    db1_lookup = Path(db1).with_suffix(".lookup")
    tmp_dir = tempfile.TemporaryDirectory(dir=".")
    tmp_dir_path = Path(tmp_dir.name)
    missing = tmp_dir_path / Path("missing.ids")
    subDB = tmp_dir_path / Path("subDB")

    # check lookup file is present
    assert db1_lookup.is_file(), f"Can't find {db1_lookup}"

    # get headers in DB1 missing from DB2
    db1_dict = Parsers.mmseqs.get_mmseqs_internal_ids_as_dict(db1_lookup)  # internal id : header
    tsv_dict = Parsers.mmseqs.mmseqs_tsv_to_dict(tsv)  # member : representative

    missing_ids = db1_dict.values() - tsv_dict.keys()
    # print(f"num missing ids: {len(missing_ids)}")

    # run mmseqs to get them as fasta
    # incredible: it took me ~4 hours to debug,
    # but if the subprocesses are underneath the 'with' block
    # some sequences get missed!
    # strangest error I've encountered in python to date
    with missing.open("w") as o:
        for internal_id in missing_ids:
            print(internal_id, file=o)

    subprocess.run(
        f"mmseqs createsubdb {missing} {db1} {subDB} --id-mode 1",
        shell=True,
        check=True,
    )

    subprocess.run(f"mmseqs convert2fasta {subDB} {output}", shell=True, check=True)
