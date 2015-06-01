#!/usr/bin/env python
__description__ = \
"""
Converts the aLTR scores in a .newick file to likelihood ratios.
"""
__author__ = "Michael J. Harms"
__date__ = "091219"

class ConvertALTRError(Exception):
    """
    General error class for this module.
    """

    pass


def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        tree_file = argv[0]
    except IndexError:
        err = "Incorrect number of arguments!\n\n%s\n\n"
        raise ConvertALTRError(err)

if __name__ == "__main__":
    main()
