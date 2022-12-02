Example Read Files
==================

This script reads the contents of both ``example_write*.py`` scripts
(which should both be run prior to running this script) and prints the file contents::

  > python example_read.py

  Opening file(s): foo*.pdb
  x: 3
  tm: [0, 1, 2]
  y: [0.1, 0.2, 0.3]
  z: [array([0, 1, 2]), array([1, 2, 3]), array([2, 3, 4])]
  a: 10

  Opening file(s): foo.h5
  a: 10
  tm: [0, 1, 2]
  x: 3
  y: [0.1, 0.2, 0.3]
  z: [array([0, 1, 2]), array([1, 2, 3]), array([2, 3, 4])]

It is located at ``[QND repo]/qnd/examples/example_read.py``.

.. literalinclude:: ../../examples/example_read.py
  :linenos:
  :language: python
