"""Low level functions to create and write a PDB file."""
from __future__ import absolute_import


def flusher_for(root):

    def _flusher(f):
        return flusher(f, root)

    return _flusher


def initializer_for(root):

    def _initializer(f):
        return initializer(f, root)

    return _initializer


def initializer(f, root):
    pass


def flusher(f, root):
    pass
