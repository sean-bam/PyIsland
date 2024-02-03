from pathlib import Path
import subprocess


def update_hhsuite_db(old_db, new_db, datfile):
    "removes entries from HHsuite DB"

    # for file in Path('../data/interim/profiles/aln1/').glob('aln1HHdb*.ff*'):
    #    print(file)

    old_db_ffindex = Path(old_db).with_suffix(".ffindex")
    old_db_ffdata = Path(old_db).with_suffix(".ffdata")

    assert (
        old_db_ffindex.is_file() and old_db_ffdata.is_file()
    ), f"""
    Cant find {old_db_ffindex} and/or {old_db_ffdata}
    """

    new_db_ffindex = Path(new_db).with_suffix(".ffindex")
    new_db_ffdata = Path(new_db).with_suffix(".ffdata")

    # make a copy of the oldDB
    shutil.copyfile(old_db_ffindex, new_db_ffindex)
    shutil.copyfile(old_db_ffdata, new_db_ffdata)

    # remove entries
    subprocess.run(
        f"ffindex_modify -s -u -f {datfile} {new_db_ffindex}", shell=True, check=True
    )

    # check
    total_ids = 0
    with open(old_db_ffindex) as infile:
        for line in infile:
            total_ids += 1

    num_ids_to_remove = 0
    with open(datfile) as infile:
        for line in infile:
            num_ids_to_remove += 1

    remaining_ids = 0
    with open(new_db_ffindex) as infile:
        for line in infile:
            remaining_ids += 1

    assert (
        total_ids - num_ids_to_remove == remaining_ids
    ), f"""
    There are {total_ids} in the HHsuite DB, and you asked to remove {num_ids_to_remove},
    leaving {total_ids - num_ids_to_remove} IDs in the new HHSuite DB.
    However, there are {remaining_ids} in the new DB, so something is off.
    Did you check the ID extension in the datfile? 
    It should probably end in .a3m or .hhm
    """


def get_hhsuite_neff(msa):
    p1 = subprocess.run(
        f"hhalign -i {msa} -o /dev/null -diff inf -id 100 -cpu 12",
        shell=True,
        capture_output=True,
        text=True,
    )

    neff = float(p1.stdout.strip().split()[10])
    return round(neff, 2)

def align_clusters_with_hhalign(input_dir, output_dir, df_clustering):
    """
    accepts a table of profile, length and cluster assignment
    runs hhalign of the longest query versus all templates in the cluster
    singletons are converted to a3m
    """
    input_dir = Path(input_dir)
    assert input_dir.is_dir(), f"{input_dir} isnt a directory!"

    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True)

    # df = pd.read_csv(clustering_table)
    grouped = df_clustering.groupby("cluster")

    # make a

    for name, group in grouped:
        output_a3m = output_dir / Path(f"cls_{name}").with_suffix(".a3m")

        # sort by profile length≈
        profiles = group.sort_values(by="qlen", ascending=False).qseq.tolist()

        # get the first, longest profile
        query_profile = input_dir / Path(f"{profiles[0]}.afa")

        # for singletons, just convert to a3m
        if len(group) < 2:
            subprocess.run(
                f"reformat.pl fas a3m {query_profile} {output_a3m}",
                shell=True,
                check=True,
            )
        else:
            # make a string of all the remaining profiles to pass to hhalign
            templates = ""
            for profile in profiles[1:]:
                subject = input_dir / Path(f"{profile}.afa")

                assert subject.is_file(), f"cant find {subject}"
                templates = templates + f"-t {subject} "

            # run hhalign
            subprocess.run(
                f"hhalign -i {query_profile} {templates} -M 50 -glob -id 100 -diff inf -oa3m {output_a3m}",
                shell=True,
                check=True,
            )