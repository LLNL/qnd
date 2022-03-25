# Misc. qnd example, runs for pdb, h5, and nc3:
# though some things don't work for specific file types

from pylab import *
from qnd.h5f import openh5
from qnd.ncf import opennc
from qnd.pdbf import openpdb
from qnd.frontend import QList

########## Write a new file ##########
iftype = 'pdb'

if iftype == 'h5':
    fname = 'example_misc.h5'
    f = openh5(fname, 'w')
elif iftype == 'nc':
    fname = 'example_misc.nc3'
    f = opennc(fname, 'w')
elif iftype == 'pdb':
    fname = 'example_misc.pdb'
    f = openpdb(fname, 'w')

# non-record vars
x = linspace(0, 10, 201)   # numpy array
y = sin(x)

f.x    = x     # like an object attribute
f['y'] = y     # like a dict

# record vars
f.recording(1)

f.time = 0
f.y    = 1.2
f.time = 1
f.y    = 2.67
# z was not declared before f.recording(), so
# netCDF unhappy: all vars must be declared prior
if iftype != 'nc':  f.z = linspace(-1, 2, 21)

f.close()

########## Read from it ##########
if iftype == 'h5':
    f = openh5(fname)
elif iftype == 'nc':
    f = opennc(fname)
elif iftype == 'pdb':
    f = openpdb(fname)

# find filename programmatically.  Long email from Munro in June 21
# esp. about multi-file family
# 24mar22: Both work for pdb and nc3, neither for h5:
if iftype != 'h5':
    print('filename method 1: ', f._qnd_group.handle.filename())
    print('filename method 2: ', f.root()._qnd_group.handle.filename())

print('Variables in f: list(f):')
print(list(f))

# non-record var: just an array
print('x: non-record', f.x)
print('type(x) :', type(f.x))

# record var: QList, not array
print('time: record', f.time)
print('type(time) :', type(f.time))
print('Is it a QList?', type(f.time) is QList)

# record var --> plain list
t = f.time[:]

f.close()   # comment so you can examine
