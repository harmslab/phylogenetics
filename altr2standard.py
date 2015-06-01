#!/usr/bin/env python

import sys

class Altr2StandardError(Exception):
    """
    General error class for this module.
    """

    pass


def altr2standard(lines):

    final_out = []
    for l in lines:
        out = []
        record = True
        for c in l:
            if c == ")":
                out.append(c)
                record = False
            elif c == ":" and not record:
                record = True

            if record:
                out.append(c)

        final_out.append("%s\n" % "".join(out))

        return "".join(final_out)


def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]
    
    try:
        filename = argv[0]
    except IndexError:
        print __usage__
        sys.exit()

    f = open(filename)
    lines = f.readlines()
    f.close()

    print altr2standard(lines),

if __name__ == "__main__":
    main()
