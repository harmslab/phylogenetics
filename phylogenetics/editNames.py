#!/usr/bin/env python
__description__ =\
"""
editNames.py

Goes through a file looking for a set of strings, replacing each one with a
specific counterpart.  The strings are defined in a delimited text file
with columns naemd "key" and "value".
"""
__author__ = "Michael J. Harms"
__date__ = "091205"
__usage__ = "editTreeNames.py file_to_modify master_file key_col value_col"

import os, sys, re


class EditNamesError(Exception):
    """
    General error class for this module.
    """

    pass

class NameObject:
    """
    """

    def __init__(self,line,column_names,column_delimiter):
        """
        Hold all column values for this sequence, keyed to column_name.
        """

        column_values = [c.strip() for c in line.split(column_delimiter)]

        if len(column_values) > len(column_names):
            warning = "There are more data columns than data names for this\n"
            warning += "line.  This usually occurs if you added a '%s' within\n" \
                    % column_delimiter
            warning += "one of the column entries.  This may lead to wonkiness.\n"
            warning += "The offending line is: \n\n"
            warning += "%s\n" % line
            warning += "which, when split, yields:\n\n"
            warning += "\n".join(["%s: %s" % (column_names[i],column_values[i])
                                 for i in range(len(column_names))])
            warning += "\n\n"

            sys.stderr.write(warning)

        try:
            self.columns = dict([(column_names[i],column_values[i])
                                 for i in range(len(column_names))])
        except IndexError:
            err = "There is an error in the following line:\n\n"
            err += line
            err += "\n\nIncorrect number of columns?\n"

            raise EditNamesError(err)

        # make it so all column values can be accessed by a simple
        # s.internal_name nomenclature
        self.__dict__.update(self.columns)


def checkUniqueness(some_list):
    """
    Verifies that every entry in a list is unique.  If it is not, it returns
    non-unique values.
    """

    unique_list = dict([(x,[]) for x in some_list]).keys()
    repeated_entries = [u for u in unique_list
                        if len([s for s in some_list if s == u]) > 1]

    return repeated_entries

def readMasterFile(name_file,column_delimiter="\t"):
    """
    Read a delimited (usually comma-delimited) file that has a set of sequence
    attributes under a unique internal_name.
    """

    f = open(name_file,'r')
    lines = f.readlines()
    f.close()

    # Parse the file
    lines = [l for l in lines if l[0] != "#" and l.strip() != ""]

    # Create a dictionary that keys column names to column numbers
    column_names = [c.strip() for c in lines[0].split(column_delimiter)]
    column_dict = dict([(c,i) for i, c in enumerate(column_names)])

    # Make sure file will be useful...
    if "internal_name" not in column_names:
        err = "\nYou must have an 'internal_name' column in this file!\n\n"
        raise EditNamesError(err)

    # Make sure column names are not repeated more than once
    for k in column_dict.keys():
        num_col_in_file = len([c for c in column_names if c == k])
        if num_col_in_file > 1:
            err = "column '%s' occurs more than once in the file!\n" % k
            raise EditNamesError(err)

    # Load all names
    names = []
    for l in lines[1:]:
        names.append(NameObject(l,column_names,column_delimiter))

    # make sure all internal_names are unique
    internal_names = [n.internal_name for n in names]
    repeated_names = checkUniqueness(internal_names)
    if len(repeated_names) != 0:
        err = "internal_name column must have unique entry for every line!\n"
        err += "The following entries are repeated:\n\n"
        err += "\n".join(repeated_names)
        err += "\n\n"

        raise EditNamesError(err)
    return names


def modifyFile(file_to_modify,names,key_column,value_column):
    """
    Read a file and replace all instances of the key_column with value_column
    where key_column and value_column are defined uniquely for each sequence
    in names.
    """

    f = open(file_to_modify)
    contents = f.read()
    f.close()

    # Grab keys and values from every sequence
    keys = []
    values = []
    for n in names:
        try:
            keys.append(n.columns[key_column])
            values.append(n.columns[value_column])
        except KeyError:
            err = "Sequences '%s' does not have '%s' or '%s'!\n\n" % \
                (n.internal_name,key_column,value_column)
            raise EditNamesError(err)

    # Check keys and values to make sure they are unique
    repeated_keys = checkUniqueness(keys)
    if len(repeated_keys) != 0:
        err = "Column '%s' has non-unique entries!\n" % key_column
        err += "The following entries are repeated:\n\n"
        err += "\n".join(repeated_keys)
        err += "\n\n"
    repeated_values = checkUniqueness(values)
    if len(repeated_values) != 0:
        err = "Column '%s' has non-unique entries!\n" % value_column
        err += "The following entries are repeated:\n\n"
        err += "\n".join(repeated_values)
        err += "\n\n"

    # Use regular experessions to replace all instances of key with value in the
    # file.
    name_dictionary = dict(zip(keys,values))
    for key in name_dictionary.keys():
        k = re.compile(key)
        num_counts = len(k.findall(contents))
        contents = k.sub("%s" % name_dictionary[key],contents,count=num_counts)

    f = open(file_to_modify, "w")
    contents = f.write(contents)
    f.close()
    return contents

def main(argv=None):
    """
    Read the command line and master file, then alter contents of
    file_to_modify and print.
    """

    if argv == None:
        argv = sys.argv[1:]
    try:
        file_to_modify = argv[0]
        master_file = argv[1]
        key_column = argv[2]
        value_column = argv[3]
    except IndexError:
        err = "Incorrect number of arguments!\n\nUsage:\n\n%s\n\n" % __usage__
        raise EditNamesError(err)

    names = readMasterFile(master_file)
    out = modifyFile(file_to_modify,names,key_column,value_column)

    print out

if __name__ == "__main__":
    main()
