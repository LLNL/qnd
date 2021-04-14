'''
example_read.py opens and lists the contents of files produced by
example_write_pdb.py and example_write_h5.py
'''
import numpy as np
from qnd.lazy import openb
from qnd.frontend import QList

fnames = ["foo*.pdb", "foo.h5"]

for fname in fnames:
    print("\nOpening file(s): {}".format(fname))

    with openb(fname, "r") as f:
        for k in list(f):
            val = f[k]
            if type(f[k]) is QList:
                val = val[:]
            print(f"{k}: {val}")
