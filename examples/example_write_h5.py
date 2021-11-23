'''
example_write_h5.py demos how to write single values and multiple records 
into an hdf5 file.
It writes to the current directory, so run it from one where you have permission to write.
'''
import numpy as np
from qnd.h5f import openh5

fname = "foo.h5"

with openh5(fname, "w") as f:

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
    
# reopen file to append data:
with openh5(fname, "a") as f:
    
    # append a single new variable to the file
    f.a = 10

    # write one new record:
    f.recording(1)
    f.tm = 2
    f.y = 0.3
    f.z = np.arange(3)+2

