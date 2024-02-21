from pathlib import Path
import random
import string

from google.cloud import bigquery
from atxlib import GCP


def upload_gff_to_bq(gff, table_id):
    """
    Accepts a GFF file and a BQ location
    """

    client = bigquery.Client()

    # convert GFF to CSV
    gff = Path(gff)
    csv = gff.with_suffix(".csv")
    df = utils.gff_to_df(gff)
    df.to_csv(csv, index=False, header=False)

    # Set table_id to the ID of the table to create.
    # table_id = self.project + "." + self.dataset + "." + table_id

    job_config = bigquery.LoadJobConfig(
        schema=[
            bigquery.SchemaField("contig", "STRING"),
            bigquery.SchemaField("source", "STRING"),
            bigquery.SchemaField("feature_type", "STRING"),
            bigquery.SchemaField("start", "INT64"),
            bigquery.SchemaField("stop", "INT64"),
            bigquery.SchemaField("score", "STRING"),
            bigquery.SchemaField("strand", "STRING"),
            bigquery.SchemaField("phase", "STRING"),
            bigquery.SchemaField("attributes", "STRING"),
        ],
    )

    load_job = client.load_table_from_file(
        csv.open("rb"), table_id, job_config=job_config
    )  # Make an API request.

    load_job.result()  # Waits for the job to complete.
    destination_table = client.get_table(table_id)  # Make an API request.
    print("Loaded {} rows.".format(destination_table.num_rows))

    # cleanup
    csv.unlink()


class bq_tools:
    def __init__(self, project, dataset, bucket):
        self.project = project
        self.dataset = dataset
        self.bucket = bucket

    def check_inputs(self):
        try:
            client.get_dataset(self.dataset)  # Make an API request.
            print("Dataset {} already exists".format(self.dataset))
        except NotFound:
            print("Dataset {} is not found, creating it".format(self.dataset))
            dataset = bigquery.Dataset(self.project + "." + self.dataset)
            dataset.location = "us-central1"
            dataset = client.create_dataset(dataset, timeout=30)


def load_table(uri: str, table_id: str):
    """
    Accepts a bucket URI pointing to a single column table
    of protein accessions

    Uploads to bigquery as the given table_id

    returns the number of uploaded rows
    """

    client = bigquery.Client()

    job_config = bigquery.LoadJobConfig(
        schema=[bigquery.SchemaField("protein_id", bigquery.SqlTypeNames.STRING)],
        write_disposition="WRITE_TRUNCATE",  # overwrite table if exists
        # skip_leading_rows=1,
    )

    # Make an API request.
    load_job = client.load_table_from_uri(uri, table_id, job_config=job_config)

    load_job.result()  # Waits for the job to complete.

    destination_table = client.get_table(table_id)  # Make an API request.
    num_rows = destination_table.num_rows

    return num_rows


def left_join(left_table: str, right_table: str, output_table: str):
    client = bigquery.Client()

    job_config = bigquery.QueryJobConfig(destination=output_table)

    sql = f"""
    SELECT
        seq.protein_id,
        seq.seq,
    FROM 
        {left_table}
        JOIN {right_table}
        USING (protein_id)
    """

    query_job = client.query(sql, job_config=job_config)
    query_job.result()


def csv2fasta(csv, output):
    i = 0
    with open(output, "w") as o:
        with open(csv) as infile:
            for line in infile:
                if line.startswith("protein_id"):
                    continue
                header, seq = line.strip().split(",")

                print(f">{header}", seq, sep="\n", file=o)

                i += 1
    return i


def write_table_to_gcs(input_table: str, uri: str, location="us-central1"):
    client = bigquery.Client()
    job_config = bigquery.ExtractJobConfig()
    # job_config.destination_format = bigquery.DestinationFormat.CSV
    job_config.print_header = False

    extract_job = client.extract_table(
        input_table,
        uri,
        location=location,
        job_config=job_config,
    )  # API request
    extract_job.result()  # Waits for job to complete.


def get_seqs_from_bq(
    uri_of_ids: str,
    bq_seq_table: str,
    bq_tmp_dataset: str,
    output_uri: str,
    location="us-central1",
    keep_local=True,
):
    """
    !!! DEPRECATED !!!
    Accepts
     - a URI of a single-column text file with protein IDs (e.g., gs://bucket/path/to/file.txt)
     - a string of a BigQuery table containing two columns of "protein_id" and "seq" (e.g., 'project.dataset.seq')
     - a string of a BigQuery dataset to write temporary tables (e.g., 'project.TMP')
     - a URI of the output fasta

    Merges the protein_IDs with their Seqs, writes a fasta file to current working directory and bucket

    The local file can be deleted with keep_local=False
    """

    client = bigquery.Client()

    # set paths
    random_str = "".join(random.choices(string.ascii_uppercase, k=4))
    bq_tmp_table1 = bq_tmp_dataset + "." + random_str
    bq_tmp_table2 = bq_tmp_dataset + "." + random_str + "2"

    bucket_name, output_blob_path = GCP_utils.extract_bucketname_and_blob_path_from_uri(
        output_uri
    )
    output_blob_name = output_blob_path.split("/")[-1]  # e.g. 'output.faa'
    tmp_blob_name = output_blob_name + ".tmp"  # e.g. 'output.faa.tmp'
    tmp_blob_path = output_blob_path + ".tmp"  # e.g. 'path/to/output.faa.tmp'
    tmp_blob_uri = (
        "gs://" + bucket_name + "/" + tmp_blob_path
    )  # e.g. gs://bucket/path/to/output.faa.tmp

    # load table to temporary BQ
    print(f"Loading {uri_of_ids} to BQ as {bq_tmp_table1}")
    num_rows = load_table(uri_of_ids, bq_tmp_table1)

    # join with seq table
    print(
        f"Joining {bq_tmp_table1} with {bq_seq_table} and writing the output to {bq_tmp_table2}"
    )
    left_join(bq_tmp_table1, bq_seq_table, bq_tmp_table2)

    # write CSV to bucket
    print(f"Extracting {bq_tmp_table2} to {tmp_blob_uri}")
    extract_job = client.extract_table(
        bq_tmp_table2,
        tmp_blob_uri,
        location=location,
    )
    extract_job.result()

    # download CSV
    print(f"Downloading {tmp_blob_uri} to {tmp_blob_name}")
    GCP_utils.download_blob(bucket_name, tmp_blob_path, tmp_blob_name)

    # convert to fasta
    print(f"Converting {tmp_blob_name} to {output_blob_name}")
    num_seqs = csv2fasta(tmp_blob_name, output_blob_name)

    if num_seqs < num_rows:
        print(
            f"""
        Warning: {uri_of_ids} has {num_rows} rows in BigQuery, 
        but only {num_seqs} sequences were written to {output_blob_name}.
        Some sequences may be missing in the output file
        """
        )

    # push to bucket
    print(f"Uploading {output_blob_name} to {output_uri}")
    GCP_utils.upload_blob(bucket_name, output_blob_name, output_blob_path)

    # cleanup
    print(f"Deleting temporary BQ tables")
    client.delete_table(bq_tmp_table1)
    client.delete_table(bq_tmp_table2)

    print(f"Deleting temporary CSVs {tmp_blob_uri} and {tmp_blob_name}")
    GCP_utils.delete_blob(bucket_name, tmp_blob_path)
    Path(tmp_blob_name).unlink()

    if not keep_local:
        print(f"Deleting local output {output_blob_name}")
        Path(output_blob_name).unlink()


def blast_rep_to_member(
    seq_table: str, clustering_table: str, blast_table: str, output_table: str
):
    """
    Accepts a bigquery table of sequences, a bigquery table of clustering information, and a bigquery table of blast hits
    and returns a table of sequences that have a blast hit to a representative sequence
    """
    client = bigquery.Client()

    # extract seqs to table
    job_config = bigquery.QueryJobConfig(
        dry_run=False,
        destination=output_table,
        write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
    )

    sql = f"""
        SELECT *
        FROM `{seq_table}`
        WHERE `protein_id` IN (
            SELECT `member`
            FROM `{clustering_table}`
            WHERE `representative` IN (
                SELECT `qseqid`
                FROM `{blast_table}`
            )
        )
    """

    query_job = client.query(sql, job_config=job_config)
    query_job.result()

    destination_table = client.get_table(output_table)
    print("Loaded {} rows.".format(destination_table.num_rows))


def get_seqs_from_bq_as_fasta(seq_table: str, output_uri: str, keep_local=True):
    """
    Accepts a bigquery table with two rows "protein_id" and "sequence" and returns a fasta file
    """

    # set paths
    (
        bucket_name,
        output_blob_path,
    ) = GCP.storage.extract_bucketname_and_blob_path_from_uri(output_uri)
    output_blob_name = output_blob_path.split("/")[-1]  # e.g. 'output.faa'
    output_blob_dir = (
        "/".join(output_blob_path.split("/")[:-1]) + "/"
    )  # e.g. 'path/to/folder/'
    tmp_blob_dir = output_blob_dir + "tmp_blobs/"  # e.g. 'path/to/folder/tmp_blobs/'
    tmp_blobs = (
        "gs://" + bucket_name + "/" + tmp_blob_dir + "split*"
    )  # e.g. 'gs://bucket/path/to/folder/tmp_blobs/split*'

    print(f"Writing {seq_table} to {tmp_blobs}")
    write_table_to_gcs(seq_table, tmp_blobs)

    # download CSV files from the directory
    print(f"Downloading files from {tmp_blobs} to current directory")
    blobs = GCP.storage.list_blobs_with_prefix(bucket_name, tmp_blob_dir)
    for blob in blobs:
        tmp_blob_name = blob.name.split("/")[-1]
        tmp_blob_path = blob.name
        GCP.storage.download_blob(bucket_name, tmp_blob_path, tmp_blob_name)

    # convert to fasta
    num_seqs = 0
    with open(output_blob_name, "w") as outfile:
        for split in Path(".").glob("split*"):
            # print(f"Converting {split} to {split}.faa")
            i = csv2fasta(split, split.with_suffix(".faa"))
            num_seqs += i
            with open(split.with_suffix(".faa")) as infile:
                for line in infile:
                    print(line.strip(), file=outfile)
    print(f"Converted {num_seqs} rows to fasta sequences")

    # push to bucket
    print(f"Uploading {output_blob_name} to {output_uri}")
    GCP.storage.upload_blob(bucket_name, output_blob_name, output_blob_path)

    # cleanup
    print(f"Cleaning up tmp files")
    for split in Path(".").glob("split*"):
        split.unlink()

    if not keep_local:
        Path(output_blob_name).unlink()

    # remove temporary CSV files from bucket
    print(f"Removing {tmp_blobs} from bucket")
    blobs = GCP.storage.list_blobs_with_prefix(bucket_name, tmp_blob_dir)
    for blob in blobs:
        GCP.storage.delete_blob(bucket_name, blob.name)
