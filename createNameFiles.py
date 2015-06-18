#!/usr/bin/env python3
__description__ = \
"""
Create a name file and sanitized fasta file from a uniprot fasta file. 

A name file can be used in conjunction with editNames.py to convert the strings
within text files between human-readable and phyml readable formats.  The
sanitized fasta file is ready to be edited with editNames.py.  These are written
out as fileroot_names.txt and fileroot_names.fasta
"""
__author__ = "Michael J. Harms"
__date__ = "2014-07-26"
__usage__ = "createNameFiles.py fasta_file_from_uniprot"

import sys, re

def parseUniprotLine(line):
    """
    Parse a uniprot fasta file line, returning the uniprot id and species.  This
    is a rather hacked parser that is brittle to changes in the default uniprot
    fasta file header.  
    """

    cols = line.split("|")
    uniprot_id = cols[1]

    description = cols[2]
    species_tmp = description.split("OS=")[1]
    species_tmp = species_tmp.split("=")[0]
    species = " ".join(species_tmp.split(" ")[:-1])

    return uniprot_id, species

def parseGenericLine(line):
    """
    Grab the name of a sequence without any (known) structure in the name.  
    Clean it up so downstream newick files don't choke. 
    """

    out_line = line[1:].strip()
    out_line = re.sub(":","-",out_line)
    out_line = re.sub(",","-",out_line)
    out_line = re.sub("\(","\[",out_line)
    out_line = re.sub("\)","\]",out_line)
   
    return out_line, "unk" 
   

def createMasterFile(fasta_file,delim="\t",parse_type="generic"):
    """
    Walk through ">" entries in the fasta file, extract the uniprot id and
    source species name, and then create a name file line from that.
    """

    parsers = {"generic":parseGenericLine,
               "uniprot":parseUniprotLine}

    parser = parsers[parse_type]

    # Create output list    
    name_out = [delim.join(["number","id","species","internal_name","pretty_name"])]
    fasta_out = []

    # read every line in the fasta file
    counter = 0
    with open(fasta_file) as f:
        for line in f:

            # If this is a header line
            if line.startswith(">"):

                id_string, species = parser(line)

                internal_name = "XX%s" % (str(counter).zfill(8))
                pretty_name = "%s-%s" % (id_string,species)
  
                name_out.append(delim.join([str(counter),id_string,species,
                                       internal_name,pretty_name]))
                
                fasta_out.append(">%s" % pretty_name)

                counter += 1

            else:
                fasta_out.append(line.strip())

    return name_out, fasta_out
    

def main(argv=None):
    """
    Main function, parses command line and creates master file.
    """
 
    if argv == None:
        argv = sys.argv[1:]
    
    try:
        uniprot_fasta = argv[0]
    except IndexError:
        err = "Incorrect arguments. Usage:\n\n%s\n\n" % __usage__
        raise IndexError(err)

    try:
        parser_type = argv[1]
    except IndexError:
        parser_type = "generic"

    name_out, fasta_out = createMasterFile(uniprot_fasta)

    # Strip extension
    file_root =".".join( uniprot_fasta.split(".")[:-1])

    # Write out name database and fasta file with those names
    f = open("%s_name.txt" % file_root,'w')
    f.write("\n".join(name_out))
    f.close()

    f = open("%s_name.fasta" % file_root,'w')
    f.write("\n".join(fasta_out))
    f.close()


# If this is called from the command line
if __name__ == "__main__":
    main()

