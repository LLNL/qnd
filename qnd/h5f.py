"""QnD wrapper for h5py HDF5 interface."""
from __future__ import absolute_import

# HDF5 conventions:
#
# HDF5 Group --> QnD QGroup, except when __class__ item present and
#   __class__ = 'list' is a QnD Qlist, with item names _0, _1, _2, ...
#   We use the qnd.QnDList class to implement the QList backend.
# HDF5 Dataset --> QnD QLeaf, except leading UNLIMITED dimension is QList
# HDF5 Dataset with UNLIMITED leading dimension --> QnD Qlist
#   A dataset of type u1 with __dtype__ attribute 'b1' is boolean.
#   The __dtype__ attribute of a dataset in general, if present, is the
#   name of its Lading typedef.  (For primitive types, this will only
#   be present for boolean.)
# HDF5 Reference Dataset --> QnD Qlist, recognized but never created
# HDF5 Named Datatype --> Lading typedef, except if __class__ attribute:
#   __none__ = 1   --> None
#   __shape__ = [len0, len1, ...] represents 0-length ndarray of given type

import sys
import re
from os.path import expanduser, expandvars

from h5py import File, Group, Dataset, Datatype

from numpy import zeros, dtype as npdtype

from .frontend import QGroup, QnDList

__all__ = ['openh5']

PY2 = sys.version_info < (3,)
_family_pattern = re.compile(r'%\d*d')


def openh5(filename, mode='r', auto=1, **kwargs):
    """Open HDF5 file using h5py, but wrapped as a QnD QGroup.

    Parameters
    ----------
    filename : str
       Name of file to open.  If filename contains %d format directive,
       open a family of files (the memb_size keyword controls file size).
    mode : str
       One of 'r' (default, read-only), 'r+' (read-write, must exist),
       'a' (read-write, create if does not exist), 'w' (create, clobber if
       exists), 'w-' (create, fail if exists).
    auto : int
       The intial state of auto-read mode.  If the QGroup handle returned
       by openh5 is `f`, then ``f.varname`` reads an array variable, but not
       a subgroup when auto=1, the default.  With auto=0, the variable
       reference reads neither (permitting later partial reads in the case
       of array variables).  With auto=2, a variable reference recursively
       reads subgroups, bringing a whole tree into memory.
    **kwargs
       Other keywords passed to h5py.File constructor.  Note that the
       driver='family' keyword is implicit if filename contains %d.

    Returns
    -------
    f : QGroup
       A file handle implementing the QnD interface.

    """
    driver = kwargs.pop('driver', None)
    if driver is None and _family_pattern.search(filename):
        driver = 'family'
    if driver is not None:
        kwargs['driver'] = driver
    if driver == 'family':
        # Make default family size 100 MB instead of 2 GB.  (??)
        # Also accept maxsize keyword in addition to memb_size.
        kwargs.setdefault('memb_size', kwargs.pop('maxsize', 134217728))
    filename = expanduser(expandvars(filename))
    return QGroup(H5Group(File(filename, mode, **kwargs)), auto=auto)


class H5Group(object):
    __slots__ = 'h5item', '__weakref__'

    def __init__(self, h5item):
        self.h5item = h5item

    @staticmethod
    def isgroup():
        return 1

    @staticmethod
    def islist():
        return 0

    isleaf = islist

    def root(self):
        h5item = self.h5item
        if h5item.name == '/':
            return self
        return H5Group(h5item['/'])

    def close(self):
        self.root().h5item.close()

    def flush(self):
        self.root().h5item.flush()

    def __len__(self):
        return len(self.h5item)

    def __iter__(self):
        return iter(self.h5item)

    def lookup(self, name):
        h5item = self.h5item
        item = h5item.get(name)
        if item is None:
            return None
        if isinstance(item, Group):
            return QnDList.fromgroup(H5Group(item))
        leaf = H5Leaf(item, h5item)
        if isinstance(item, Dataset):
            maxshape = item.maxshape
            if (maxshape and maxshape[0] is None and
                    not any(n is None for n in maxshape[1:])):
                return QnDList(leaf)
        return leaf

    def declare(self, name, dtype, shape, unlim=None):
        h5item = self.h5item
        if dtype == dict:
            return H5Group(h5item.create_group(name))
        if dtype == list:
            return QnDList(H5Group(h5item.create_group(name)), 1)
        if dtype is None or (shape and not all(shape)):
            h5item[name] = npdtype('u1') if dtype is None else dtype
            item = h5item[name]
            if dtype is None:
                item.attrs['__none__'] = True
            else:
                item.attrs['__shape__'] = shape
        elif unlim:
            item = h5item.create_dataset(name, (1,)+shape, dtype=dtype,
                                         maxshape=(None,)+shape)
        else:
            item = h5item.create_dataset(name, shape, dtype=dtype)
        item = H5Leaf(item, h5item)
        return QnDList(item, 1) if unlim else item

    def attget(self, vname):
        item = self.lookup(vname) if vname else self
        if isinstance(item, QnDList):
            item = item.parent()
        return _WrapAttributeManager(item.h5item.attrs)

    def attset(self, vname, aname, dtype, shape, value):
        item = self.lookup(vname) if vname else self
        if isinstance(item, QnDList):
            item = item.parent()
        attrs = item.h5item.attrs
        attrs.create(aname, value, shape, dtype)


class _WrapAttributeManager(object):
    __slots__ = 'attrs',

    def __init__(self, attrs):
        self.attrs = attrs

    # qnd.QAttribute uses only __iter__, get, items, __len__, __contains__

    def __contains__(self, name):
        return name in self.attrs

    def __len__(self):
        # This functionality missing entirely in underlying AttributeManager.
        return len(list(self.attrs))

    def __iter__(self):
        return iter(self.attrs)

    def get(self, key, default=None):
        return self.attrs.get(key, default)

    def items(self):
        return self.attrs.iteritems() if PY2 else self.attrs.items()


class H5Leaf(object):
    __slots__ = 'h5item', 'parent'

    def __init__(self, h5item, parent):
        self.h5item = h5item
        self.parent = parent

    @staticmethod
    def isleaf():
        return 1

    @staticmethod
    def isgroup():
        return 0

    islist = isgroup

    def root(self):
        return H5Group(self.parent['/'])

    def query(self):
        # return dtype, shape, sshape
        h5item = self.h5item
        if isinstance(h5item, Datatype):
            shape = h5item.attrs.get('__shape__')
            if shape is None:
                if h5item.attrs.get('__none__'):
                    return None, (), ()
                else:
                    return type, (), ()
        else:
            shape = h5item.shape
        return h5item.dtype, shape, shape

    def read(self, args=()):
        h5item = self.h5item
        if isinstance(h5item, Datatype):
            shape = h5item.attrs.get('__shape__')
            if shape is None:
                if h5item.attrs.get('__none__'):
                    return None
                else:
                    return h5item.dtype
            value = zeros(shape, h5item.dtype)
        else:
            value = h5item[args]
        return value

    def write(self, value, args=()):
        item = self.h5item
        if isinstance(item, Datatype):
            return
        maxshape = item.maxshape
        if (args and maxshape and maxshape[0] is None and
                not any(n is None for n in maxshape[1:])):
            # Create next element of an UNLIMITED array now.
            try:
                i = int(args[0])
            except (TypeError, ValueError):
                pass
            else:
                shape = item.shape
                if shape[0] == i:
                    item.resize((i+1,)+shape[1:])
        item[args] = value
