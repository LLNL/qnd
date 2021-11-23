'''
example_write_pdb.py demos how to write single values and multiple records 
into a family of pdb files.
It writes to the current directory, so run it from one where you have permission to write.
'''
import numpy as np
from qnd.pdbf import openpdb

with openpdb("foo000.pdb", "w") as f:

    # write a single value:
    f.x = 3

    # turn on recording and write two records:
    f.recording(1)
    f.tm = 0
    f.y = 0.1
    f.z = np.arange(3)

    f.tm = 1
    f.y = 0.2
    f.z = np.arange(3)+1

# append a single new variable to the file:
with openpdb("foo000.pdb", "a") as f:
    f.a = 10

# open a new file in the family and write one new record:
with openpdb("foo001.pdb", "w") as f:
    f.recording(1)
    f.tm = 2
    f.y = 0.3
    f.z = np.arange(3)+2

