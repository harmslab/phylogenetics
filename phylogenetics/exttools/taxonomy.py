"""Ping the NCBI Taxonomy web service and retrieve data about a sequence.

Stupidly simple module... probably needs more intelligent/safe parsing of data.
"""
import re
import requests

def query(keyword, type="name", **kwargs):
    """Send a query to retrieve taxonomic data about a set of sequences.
    """
    main_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?"

    # Construct arguments for query.
    kwargs[type] = keyword
    end_url = "&".join([key+"="+kwargs[key] for key in kwargs])

    total_url = main_url + end_url

    data = requests.get(total_url)

    # Strip out the HTML with the taxonomic lineage
    regex = re.compile('<[aA][^<]+</[aA]>')
    lineage = regex.findall(data.text)
    # Get taxonomic tags
    ##### Result looks like a list of: #####
    # <a ALT="superkingdom"
    # href="/Taxonomy/Browser/wwwtax.cgi?mode=Undef&amp;id=2759&amp;lvl=3&amp;keep=1&amp;srchmode=1&amp;unlock"
    # TITLE="superkingdom">Eukaryota
    #</a>

    taxonomy = {}
    for level in lineage:
        # Get the classification
        try:
            classification = re.compile('TITLE="([A-Za-z\s]+)"').search(level).group(1)

            # Get the named classification
            try:
                label = re.compile('>([\w ]+)<').search(level).group(1)
            except AttributeError:
                raise Exception("Found classification, but cannot find ")

            # Add classification to taxonomy dictionary
            taxonomy[classification] = label

        # If the tag pulled from last REGEX is not a classification (but
        # another <a> tag), skip it.
        except AttributeError:
            pass

    return taxonomy
