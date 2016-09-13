from __future__ import absolute_import

import re
import requests

__doc__ = """Ping the NCBI Taxonomy web service and retrieve data about a sequence.

Stupidly simple module... probably needs more intelligent/safe parsing of data.
"""
def query(keyword, type="name", **kwargs):
    """Send a query to retrieve taxonomic data about keyword from NCBI Taxonomy
    Web service

    Parameters
    ----------
    keyword : str
        query string for BLAST Taxonomy. Can be common names, scientific names,
        taxonomic ID, etc.
    type : string
        type of query. either 'name' or 'id'.
    **kwargs are used to add extra options to the HTTPS request.

    Returns
    -------
    taxonomy : dict
        dictionary of all taxonomic keys returned by query.
    """
    main_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?"

    # Construct arguments for query.
    kwargs[type] = keyword
    end_url = "&".join([key+"="+kwargs[key] for key in kwargs])
    total_url = main_url + end_url

    # Send HTTPS request to NCBI
    data = requests.get(total_url)

    # initialize taxonomy datatype
    taxonomy = {}

    #########  Strip out the HTML with the taxonomic lineage ##########
    regex = re.compile('<[aA][^<]+</[aA]>')
    lineage = regex.findall(data.text)
    # Get taxonomic tags
    ##### Result looks like a list of: #####
    # <a ALT="superkingdom"
    # href="/Taxonomy/Browser/wwwtax.cgi?mode=Undef&amp;id=2759&amp;lvl=3&amp;keep=1&amp;srchmode=1&amp;unlock"
    # TITLE="superkingdom">Eukaryota
    #</a>

    for level in lineage:
        # Get the classification

        try:
            classification = re.compile('TITLE="([A-Za-z\s]+)"').search(level).group(1)

            # Get the named classification
            try:
                label = re.compile('>([\w ]+)<').search(level).group(1)

                # Add classification to taxonomy dictionary
                taxonomy[classification] = label
            except AttributeError:
                pass

        # If the tag pulled from last REGEX is not a classification (but
        # another <a> tag), skip it.
        except AttributeError:
            pass

    return taxonomy
