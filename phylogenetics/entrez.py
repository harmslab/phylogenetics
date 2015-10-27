from Bio import Entrez

def download(accession_list, email, out_file, db="protein",
                      batch_download_size=50, write=False):
    """
    Download the Entrez sequences from a list of accessions.

    Arguments :
    ---------
    accession_list: list
        list of ncbi accesion numbers
    out_file:
        file in which to write output in fasta/xml format
    db:
        database to use for accession
    batch_download_size:
        size of individual download packets
    force:  True/False.
        Overwrite existing download file. If False, the program
        throws a notice that an old file is being used rather than re-
        downloading.
    """
    # Biopython is imported here...  I realize this is a bit overkill for now.
    Entrez.email = email

    #first get GI for query accesions
    query  = " ".join(accession_list)
    handle = Entrez.esearch( db=db,term=query,retmax=10**9 )
    giList = Entrez.read(handle)['IdList']

    #post GID list to NCBI
    search_handle = Entrez.epost(db=db, id=",".join(giList))
    search_results = Entrez.read(search_handle)
    webenv,query_key = search_results["WebEnv"], search_results["QueryKey"]

    # Now download the sequences (in fasta/xml format).
    count = len(accession_list)
    total_xml = ""
    for start in range(0,count,batch_download_size):
        end = min(count, start+batch_download_size)
        fetch_handle = Entrez.efetch(db=db, rettype="fasta",
                                     retmode="xml",retstart=start,
                                     retmax=batch_download_size,
                                     webenv=webenv,query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        total_xml += data + "\n"

    # Write to file is interested.
    if write:
        out_handle = open(out_file, "w")
        out_handle.write(data)
        out_handle.close()

    # parse the xml and get sequence data as list
    sequence_data = parse_entrez_xml(total_xml)

    # Map accession list to their sequences
    mapping = dict([(accession_list[i], sequence_data[i]["sequence"]) for i in range(len(accession_list))])

    return mapping

def parse_entrez_xml(xml_string):
    """
        Parse Blast's Fasta/XML formatted file, returning a list of each
        sequence's data.

        Args:
        ----------
        xml_string: str
            Fasta/XML formatted string from Blast Output.

        Returns:
        -------
        sequences: list of dicts
            List of sequences data in their own lists.
    """
    # Fix screwed up XML because sequences downloaded and output concatenated
    sequence_input = flatten_concatenated_XML(xml_string, "TSeqSet")

    # Now we should have valid XML...
    sequence_input = ET.XML(sequence_input)

    # XML Tag prefix to strip
    prefix = "TSeq_"

    sequences = []
    for i, sequence in enumerate(sequence_input):
        # Rip out all properties of a sequence
        properties = dict([(s.tag[len(prefix):],str(s.text).strip()) for s in sequence])

        # Append to sequences.
        sequences.append(properties)

    return sequences
