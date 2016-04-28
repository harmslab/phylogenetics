import time

# import server exceptions
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2

# Use the biopython entrez api
from Bio import Entrez

def download(accession_list,
    email,
    db="protein",
    batch_download_size=50):
    """
    Download the Entrez XML metadata from a list of accessions.

    Parameters
    ----------
    accession_list : list
        list of ncbi accesion numbers
    out_file :
        file in which to write output in fasta/xml format
    db :
        database to use for accession
    batch_download_size :
        size of individual download packets
    write :  True/False.
        Overwrite existing download file. If False, the program
        throws a notice that an old file is being used rather than re-
        downloading.
    """
    # Biopython is imported here...  I realize this is a bit overkill for now.

    Entrez.email = email

    # Posting too many GIDs to Entrez will break it, so chunks
    # of 500 is a reasonable size.
    total_xml = ""
    gi_size = 500
    n_gids = len(accession_list)

    for post_start in range(0,n_gids, gi_size):
        post_end = min(n_gids, post_start+gi_size)

        #first get GI for query accesions
        query  = " ".join(accession_list[post_start:post_end])
        handle = Entrez.esearch( db=db,term=query,retmax=10**9 )
        giList = Entrez.read(handle)['IdList']
        handle.close()

        # post first chunk of GID list to NCBI
        search_handle = Entrez.epost(db=db, id=",".join(giList))
        search_results = Entrez.read(search_handle)
        search_handle.close()
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]

        # Now download the sequences (in fasta/xml format).
        count = len(giList)

        # Step through this epost to Entrez and fetch all sequences in smaller batches
        for start in range(0,count,batch_download_size):
            end = min(count, start+batch_download_size)

            # Try to fetch from Entrez, If it fails three times, raise error
            try:
                fetch_handle = Entrez.efetch(db=db, rettype="fasta",
                                         retmode="xml",retstart=start,
                                         retmax=batch_download_size,
                                         webenv=webenv,query_key=query_key)

            except HTTPError:
                raise

            data = fetch_handle.read()
            fetch_handle.close()
            total_xml += data + "\n"

    return total_xml
