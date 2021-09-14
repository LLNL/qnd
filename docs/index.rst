.. QND documentation master file, created by
   sphinx-quickstart on Fri Dec  7 12:06:10 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

+------------------------+------------------------+------------------------+
|    :ref:`genindex`     |    :ref:`modindex`     |     :ref:`search`      |
+------------------------+------------------------+------------------------+

QnD package
===========

The qnd module (quick and dirty) provides a frontend for reading and
writing self-describing binary files.  Backends exist for HDF5 and PDB
file formats (the former via h5py, the latter via pure python code).
Adding backends is not very difficult; the interface is well-defined
and relatively small.

This manual describes the design philosophy behind this user interface;
in a nutshell, the idea is to keep things as simple as possible, but
no simpler (as Einstein said).  The problem we are setting out to
solve is to store collections of scientific data, which means for the
most part arrays of numbers.


Contents
========

.. toctree::

   qnd
   adict
   frontend
   generic
   h5f
   lazy
   ncf
   pdbdump
   pdbf
   pdbparse
   utils

.. toctree::
  :caption: Example

  examples/example_write_h5
  examples/example_write_pdb
  examples/example_read

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
