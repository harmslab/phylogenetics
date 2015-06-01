#!/usr/bin/env python
__description__ = \
"""
compareAncestor.py

"""
__author__ = "Michael J. Harms"
__usage__ = "comapreAncestors.py ancestor_file1 ancestor_file2"
__date__ = "100726"

import sys, phyloBase

class CompareAncestorError(Exception):
    """
    General error class for this module.
    """
    
    pass

def readAncestorFile(ancestor_file):
    """
    """

    f = open(ancestor_file,'r')
    lines = f.readlines()
    f.close()

    # Skip comments and blank lines
    lines = [l for l in lines if l.strip() != "" and l[0] != "#"]

    out = []

    num_states = (len(lines[0].split())-2)/2

    for l in lines[1:]:
        position = int(l[7:12])
        tmp_out = []
        for i in range(num_states):
            aa = l[12+12*i:18+12*i].strip()
            pp = float(l[18+12*i:24+12*i])

            tmp_out.append((aa,pp))

        out.append((position,tmp_out))

    return out

def compareAncestors(ancestor1_file,ancestor2_file,ambiguous_cutoff=0.8):
    """
    """

    anc1 = readAncestorFile(ancestor1_file)
    anc2 = readAncestorFile(ancestor2_file)

    anc1_pos = [p[0] for p in anc1]
    anc2_pos = [p[0] for p in anc2]

    only_in_anc1 = [p for p in anc1_pos if p not in anc2_pos]
    only_in_anc2 = [p for p in anc2_pos if p not in anc1_pos]

    if len(only_in_anc1) > 0:
        print "# Warning: some sites only in ancestor 1:"
        print "".join(["# %i\n" % p for p in only_in_anc1]),
    if len(only_in_anc2) > 0:
        print "# Warning: some sites only in ancestRr 2:"
        print "".join(["# %i\n" % p for p in only_in_anc2]),

    all_pos = [p for p in anc1_pos if p not in only_in_anc1] 
    all_pos.extend([p for p in anc2_pos if p not in only_in_anc2 and p not in all_pos])

    anc1_dict = dict([a for a in anc1 if a[0] in anc1_pos])
    anc2_dict = dict([a for a in anc2 if a[0] in anc2_pos])

    out = []
    out.append("# pos new_state old_state same? state_type?")
    out.append(" ambiguity pp_new pp_old\n")
    out.append("#\n# same?\n")
    out.append("#    \'*\' -> changed\n")
    out.append("#    \' \' -> no change\n")
    out.append("# flipped_with_alternate?\n")
    out.append("#    \'*\' -> took new state\n")
    out.append("#    \'~\' -> took alternate state\n")
    out.append("#    \' \' -> no change in state\n")
    out.append("# ambig_state key:\n")
    out.append("#    \'~\' -> ambiguous in both\n")
    out.append("#    \'-\' -> newly ambiguous\n")
    out.append("#    \'+\' -> newly well supported\n")
    out.append("#    \' \' -> well suppported in both\n")

    for p in all_pos:
        s1 = anc1_dict[p]
        s2 = anc2_dict[p]

        # See if the new reconstruction has the same residue at this position
        same = "*" 
        if s1[0][0] == s2[0][0]:
            same = " "

        # Check to see if new state existed as less likely state in original
        # reconstruction
        flipped = " " 
        if same == "*":
            if s1[0] in [a[0] for a in s2[1:]]:
                flipped = "~"
            else:
                flipped = "*"

        # Remained ambiguous
        if s1[0][1] <= ambiguous_cutoff and s2[0][1] <= ambiguous_cutoff:
            ambig_state = "~"

        # Newly ambiguous
        elif s1[0][1] <= ambiguous_cutoff and s2[0][1] > ambiguous_cutoff:
            ambig_state = "+"

        # Became well supported 
        elif s1[0][1] > ambiguous_cutoff and s2[0][1] <= ambiguous_cutoff:
            ambig_state = "-"

        # Remained well supported
        else:
            ambig_state = " "

        check_me = " "
        if ambig_state == "-" or \
            (same == "*" and ambig_state == " "):
            check_me = "!"

        out.append("%5i %s %s %s %s %s %6.2f%6.2f %s\n" % (p,s1[0][0],s2[0][0],
                   same,flipped,ambig_state,s1[0][1],s2[0][1],check_me))

    return "".join(out)
    

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        ancestor1_file = argv[0]
        ancestor2_file = argv[1]
    except IndexError:
        err = "Incorrect number of arguments!\n\n%s\n\n" % __usage__
        raise CompareAncestorError(err)

    out = compareAncestors(ancestor1_file,ancestor2_file)

    print out


if __name__ == "__main__":
    main() 
