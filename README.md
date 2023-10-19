Quick and Dirty (QnD): Python Binary Storage Interface
======================================================

This very simple end-user-friendly API for storing numpy ndarrays in
binary files is written as a frontend capable of supporting multiple
backend file formats.  Current backends include an h5py wrapper for
HDF5 files and pure python implementations of the PDB and netCDF3 file
formats.  The qnd API concentrates on declaring, writing, querying,
and reading numpy ndarrays (excepting dtype 'O' arrays of arbitrary
python objects).  Additionally, qnd maps python dicts (with str-valued
keys) to HDF5 groups or PDB directories, and defines its own mapping
from python lists to groups/directories with automatically generated
member names of the form _0, _1, _2, etc.  The support for both dict
and list objects means that qnd easily handles any data structures
that can be mapped to the JSON data interchange format, which is
pretty nearly anything.

Another qnd innovation is a simplified interface for storing the
history of changing arrays.  In HDF5 or PDB or netCDF3, these are
presented as arrays with an unlimited dimension, but are in fact
stored in discontiguous chunks or blocks in all three cases.  The
qnd API, by contrast, introduces a recording mode in which writing
to a variable creates its next chunk rather than overwriting its
previous value.  A complementary goto mode allows you to retrieve
chunks in the corresponding way.

Not all backends support all of the API features; for example, the
netCDF3 backend cannot handle groups or lists beyond the root group of
the entire file and the lists qnd associates with record variables.

Although the qnd interface provides a very general framework for
binary storage, informed by the excellent netCDF data model, in detail
qnd is highly pythonic.  As much as possible, qnd file, group, or list
handles behave like python dicts or lists, with the same operators,
iterators, and important method names as the workhorse python and
numpy objects they represent.

The QND documentation is available at [https://qnd.readthedocs.io/en/latest/](https://qnd.readthedocs.io/en/latest/)

Authors
-------
QnD was created by David H. Munro.


License
-------

QnD is distributed under the terms of the BSD-3 License.

All new contributions must be made under this license.

See LICENSE and NOTICE for details.

SPDX-License-Identifier: BSD-3-Clause

Lawrence Livermore tracking ID: LLNL-CODE-807802
