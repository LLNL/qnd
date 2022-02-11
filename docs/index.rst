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
writing self-describing binary files.  Backends exist for HDF5 (via h5py),
netCDF3, and PDB (via pure python code contained in QnD) file formats.
netCDF4 is not specifically supported, but since it's based on HDF5 
such files can be treated as HDF5.  Adding backends is not very difficult; 
the interface is well-defined and relatively small.

QnD is free and open-source, and hosted on `LLNL's public github <https://github.com/LLNL/qnd/>`_.

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
