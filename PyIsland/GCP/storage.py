from pathlib import Path
import random
import string

from google.cloud import storage

storage_client = storage.Client()


def get_blobs_in_bucket_as_set(bucket, prefix):
    """
    Accepts strings pointing to a bucket name and a "prefix" (i.e., a "folder")
    The prefix should end with a "/"
    Returns a python set of "filenames"

    """
    client = storage.Client()
    bucket = client.get_bucket(bucket)

    blobs = bucket.list_blobs(prefix=prefix, delimiter="/")

    blob_set = set()
    for blob in blobs:
        f = Path(blob.name)  # .e.g, path/to/file.txt
        blob_set.add(f.name)  # e.g., file.txt

    return blob_set


def download_blob(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the bucket.
    Take from
    https://github.com/googleapis/python-storage/blob/main/samples/snippets/storage_download_file.py
    """
    # The ID of your GCS bucket
    # bucket_name = "your-bucket-name"

    # The ID of your GCS object
    # source_blob_name = "storage-object-name"

    # The path to which the file should be downloaded
    # destination_file_name = "local/path/to/file"

    bucket = storage_client.bucket(bucket_name)

    # Construct a client side representation of a blob.
    # Note `Bucket.blob` differs from `Bucket.get_blob` as it doesn't retrieve
    # any content from Google Cloud Storage. As we don't need additional data,
    # using `Bucket.blob` is preferred here.
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)


def upload_blob(bucket_name, source_file_name, destination_blob_name):
    """Uploads a file to the bucket."""
    # The ID of your GCS bucket
    # bucket_name = "your-bucket-name"
    # The path to your file to upload
    # source_file_name = "local/path/to/file"
    # The ID of your GCS object
    # destination_blob_name = "storage-object-name"

    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)

    blob.upload_from_filename(source_file_name)


def delete_blob(bucket_name, blob_name):
    """Deletes a blob from the bucket."""
    # bucket_name = "your-bucket-name"
    # blob_name = "your-object-name"

    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.delete()


def list_blobs_with_prefix(bucket_name, prefix, delimiter="/"):
    """Lists all the blobs in the bucket."""
    # bucket_name = "your-bucket-name"

    from google.cloud import storage

    storage_client = storage.Client()

    # Note: Client.list_blobs requires at least package version 1.17.0.
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix, delimiter=delimiter)

    return blobs


def extract_bucketname_and_blob_path_from_uri(uri):
    bucket_name = uri.split("/")[2]
    blob_path = "/".join(uri.split("/")[3:])

    return bucket_name, blob_path


def delete_directory(bucket, prefix):
    """
    Deletes a folder in a bucket
    """
    print(f"cleaning temporary blobs in {prefix}")
    for blob in bucket.list_blobs(prefix=prefix):
        blob.delete()


class Concatenate:
    def pdb_id_generator(size=4, chars=string.ascii_uppercase + string.digits):
        """
        Generates a random 4 character string
        from:
        https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits
        """
        return "".join(random.choice(chars) for _ in range(size))

    def _run_compose(bucket, blob_name, blobs_to_compose):
        composed_blob = bucket.blob(blob_name)
        composed_blob.compose(blobs_to_compose, timeout=600)
        return composed_blob

    def concatenate_blobs(input_bucket, input_folder, output_blobname="cat"):
        """
        Concatenates all blobs the given input_folder to output_blobname

        Example usage:
        concatenate_blobs('startlit-myth-iscb_screen', "00db/jgi/prot/", 'JGI_041322.faa')

        will combine all files directly under /prot/ into "JGI_041322.faa"
        """
        client = storage.Client()
        bucket = client.get_bucket(input_bucket)

        blobs_to_compose = []
        running_size = 0

        for i, blob in enumerate(
            bucket.list_blobs(prefix=input_folder, delimiter="/"), 1
        ):
            running_size += blob.size
            blobs_to_compose.append(blob)

            if len(blobs_to_compose) == 31:
                random4letters = Concatenate.pdb_id_generator()
                blob_name = input_folder + "compose/" + random4letters

                composed_blob = Concatenate._run_compose(
                    bucket, blob_name, blobs_to_compose
                )
                blobs_to_compose = [composed_blob]

        # compose the remaining, if any.
        if len(blobs_to_compose) != 0:
            random4letters = Concatenate.pdb_id_generator()
            blob_name = input_folder + "compose/" + random4letters
            composed_blob = Concatenate._run_compose(
                bucket, blob_name, blobs_to_compose
            )

        # final operations -- rename blob name & delete temporary directory.
        assert (
            composed_blob.size == running_size
        ), f"""
        Error - the concatenated file size does not equal the sum of the input files
        """
        blob_copy = bucket.copy_blob(
            composed_blob,
            bucket,
            input_folder + output_blobname,
        )

        delete_directory(bucket, prefix=input_folder + "compose/")
