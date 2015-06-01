#!/usr/bin/env python
__description__ = \
"""
Print the length of each sequence in a fasta file.
"""


import sys

f = open(sys.argv[1])
lines = f.readlines()
f.close()

lines = [l for l in lines if l.strip() != ""]

out_lines = []
names = []
for l in lines:
    if l[0] == ">":
        out_lines.append(">")
        names.append(l[1:].strip())
    else:
        t = "".join([a for a in l if a != "-"])
        out_lines.append(t.strip())

combined = "".join(out_lines)
seq = combined.split(">")

lengths = [len(s) for s in seq]

for i in range(len(names)):
    print names[i], lengths[i]
