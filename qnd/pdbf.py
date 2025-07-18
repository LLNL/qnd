"""Pure python QnD wrapper for PDB files.

Information on the PDB file format is somewhat hard to come by.  
Try the SILO github repo https://github.com/LLNL/Silo
esp. src/pdb/ and src/score/ dirs.
Also the end of the QND file pdbparse.py has a long comment with a detailed description.

Note that yorick-generated PDB files are version II, not version III.  Also,
yorick pointers are written (by default) in a yorick-specific format.  Since
yorick readability has held back many application codes (the LEOS library
and LLNL rad-hydro codes used for ICF design), most of the dwindling legacy
PDB files are version II.  Hence, this implementation focuses on the
version III format, and the version I format is supported only for reading.

Furthermore, this implementation only supports IEEE 754 4 and 8 byte
floating point formats, since those are the only unambiguous floating
point formats supported by numpy.  Fortunately, this covers all modern
PDB files likely to show up in practice, so we have no significant
incentive to do the work required to support exotic formats.

"""
from __future__ import absolute_import

import sys
import weakref
from numbers import Integral
from collections import OrderedDict
from warnings import warn

from numpy import (zeros, arange, fromfile, prod, array, ascontiguousarray,
                   dtype as npdtype)
from numpy.core.defchararray import encode as npencode, decode as npdecode

from .frontend import QGroup, QnDList
from .generic import opener
from .pdbparse import parser, PDBChart
from .pdbdump import flusher_for, initializer_for
from .utils import leading_args

__all__ = ['openpdb']

PY2 = sys.version_info < (3,)
if PY2:
    range = xrange


def openpdb(filename, mode='r', auto=1, hooks=None, **kwargs):
    """Open PDB file or family, and wrap it in a QnD QGroup.

    Parameters
    ----------
    filename : str
       Name of file to open.  See notes below for file family.
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
    hooks : object
       Methods of this object are called to modify default openpdb behavior.
       The only hook currently supported is hooks.pre_parse(chart, symtab,
       version), where version is the PDB version (usually 2), and
       which should accept the raw byte strings for the chart and symbol
       table, modify them, and return the modified (char, symtab).
    **kwargs
       Other keywords.  The `maxsize` keyword sets the size of files in a
       family generated in ``recording==1`` mode; a new file will begin when
       the first item in a new record would begin beyond `maxsize`.  The
       default maxsize is 128 MiB (134 MB).  The `order` keyword can be '>'
       or '<' to force the byte order in a new file; by default the byte
       order is the native order.  File families always have the same order
       for every file, so `order` is ignored if any files exist.

    Returns
    -------
    f : QGroup
       A file handle implementing the QnD interface.

    Notes
    -----
    The `filename` may be an iterable, one string per file in order.  The
    sequence may extend beyond the files which actually exist for 'r+', 'a',
    'w', or 'w-' modes.

    Alternatively `filename` specifies a family if it contains shell globbing
    wildcard characters.  Existing matching files are sorted first by length,
    then alphabetically (ensuring that 'file100' comes after 'file99', for
    example).  If there is only a single wildcard group, it also serves to
    define a sequence of future family names beyond those currently existing
    for 'r+', 'a', 'w', or 'w-' modes.  A '?' pattern is treated the same as
    a '[0-9]' pattern if all its matches are digits or if the pattern
    matches no existing files.  Similarly, a '*' acts like the minimum number
    of all-digit matches, or three digits if there are no matches.

    """
    maxsize = kwargs.pop('maxsize', 134217728)
    order = kwargs.pop('order', None)
    allrecs = kwargs.pop('allrecs', None)
    atnames = kwargs.pop('atnames', None)
    if order:
        if order not in '<>':
            raise ValueError("order must be either > or <")
        order = 1 if order == '>' else 2
    kwargs['nextaddr_mode'] = 1  # tell opener to initialize nextaddr to 0
    handle, n = opener(filename, mode, **kwargs)
    root = PDBGroup(handle, maxsize)
    for i in range(n):
        try:
            parser(handle, root, i, allrecs, atnames, hooks)
        except IOError:
            # Something went terribly wrong.  If this is first file, we die.
            name = handle.filename(i)
            if not i:
                raise IOError("Fatal errors opening PDB file "
                              "{}".format(name))
            handle.open(i-1)
            warn("file family stopped by incompatible {}".format(name))
    if not n and order:
        root.chart.byteorder = order
    # If file was freshly created, setting initializer calls it.
    handle.callbacks(flusher_for(root), initializer_for(root))
    # If any files exist, parser has set nexaddr to the chart address
    # of the last existing file in the family.  If there are no record
    # variables, we let this stand.  However, if there are record variables,
    # and the family is writable, we set nextaddr to the zero address of
    # the next file beyond all existing files.  This causes any new records
    # to be placed in a new file, leaving all existing files in the
    # family undisturbed.
    mode = mode.lower()
    if ((mode.startswith('a') or mode.startswith('r+')) and
        handle.state[4] is not None and _has_records(root)):
            handle.declared(handle.zero_address(len(handle.state[2])), None, 0)
    return QGroup(root, auto=auto)


def _has_records(root):
    for name in root:
        item = root.items[name]
        if item.islist() == 1:
            return True
        # Recurse into groups, but not into lists or objects of any type.
        if item.isgroup() and '__class__' not in item and _has_records(item):
            return True
    return False


class PDBGroup(object):
    """A directory in a PDB file, or a whole file or file family.

    """
    def __init__(self, parent, maxsize=134217728):
        if not isinstance(parent, PDBGroup):  # this is root group
            self.root = weakref.ref(self)
            self.maxsize = maxsize
            self.maxblocks = 0
            self.handle = parent
            self.chart = PDBChart(self)
            self.has_attributes = False
        else:
            self.root = parent.root
        self.items = OrderedDict()
        self.attrs = None

    @staticmethod
    def isgroup():
        return 1

    @staticmethod
    def islist():
        return 0

    isleaf = islist

    def close(self):
        self.root().handle.close()

    def flush(self):
        self.root().handle.flush()

    def __len__(self):
        return len(self.items)

    def __iter__(self):
        return iter(self.items)

    def lookup(self, name):
        item = self.items.get(name)
        if isinstance(item, PDBGroup):
            item = QnDList.fromgroup(item)
        return item

    def declare(self, name, dtype, shape, unlim=None, addr=-1):
        current = self.items.get(name)
        if dtype == dict:
            if current is not None:
                if current.isgroup():
                    return current
                raise KeyError("already a non-group item {}".format(name))
            item = PDBGroup(self)
        elif dtype == list:
            if current is not None:
                if current.islist() == 2:
                    return current
                raise KeyError("already a non-list item {}".format(name))
            item = QnDList(PDBGroup(self), 1)
        else:
            if current is not None:
                raise KeyError("attempt to redeclare {}".format(name))
            if dtype is None and name == '_':
                # Assume we are creating a QList.
                dtype = npdtype('u1')
            elif isinstance(dtype, npdtype) and dtype.kind == 'S':
                # frontend never passes 'U' dtype
                shape = shape + (dtype.itemsize,)
                dtype = npdtype('S1')
            item = PDBLeaf(self, addr, dtype, shape, unlim)
            if unlim:
                item = QnDList(item, None if hasattr(addr, '__iter__') or
                               addr != -1 else 1)
        self.items[name] = item
        return item

    # This is used in pdbparse._endparse to declare or check symbols
    # against declarations from previous files in a family.
    def _leaf_declare(self, name, dtype, shape, addr):
        item = self.items.get(name)
        unlim = isinstance(addr, list)
        if item is not None:
            tsa = None
            if not isinstance(item, PDBLeaf):
                item = None if isinstance(item, PDBGroup) else item.parent()
            if item is not None:
                tsa = item.tsa
                if (unlim != isinstance(tsa[2], list) or
                        tsa[:2] != (dtype, shape)):
                    item = None
                elif unlim:
                    tsa[2].extend(addr)
            if item is None:
                raise IOError("incompatible redeclaration of {}".format(name))
            return
        self.declare(name, dtype, shape, unlim, addr)

    def attget(self, vname):
        item = self.lookup(vname) if vname else self
        if isinstance(item, QnDList):
            item = item.parent()
        attrs = item.attrs
        if attrs is None:
            item.attrs = attrs = PDBAttrs()
        return attrs

    def attset(self, vname, aname, dtype, shape, value):
        item = self.lookup(vname) if vname else self
        if isinstance(item, QnDList):
            item = item.parent()
        attrs = item.attrs
        if attrs is None:
            item.attrs = attrs = PDBAttrs()
            self.root().has_attributes = True
        if value.dtype != dtype or value.shape != shape:
            v = zeros(shape, dtype)
            v[()] = value
            value = v
        attrs[aname] = value

    def _read_attr(self, vname, aname, dtype, shape, addr):  # for pdbparse
        value = PDBLeaf(self.root(), addr, dtype, shape, None).read()
        if dtype.kind == "S" and not PY2:
            try:
                value = npdecode(value, "utf8")
            except UnicodeDecodeError:
                value = npdecode(value, "latin1")
        self.attset(vname, aname, value.dtype, value.shape, value)

    def _read_shape(self, vname, dtype, shape, addr):  # for pdbparse
        value = PDBLeaf(self.root(), addr, dtype, shape, None).read()
        item = self.lookup(vname)
        if value.dtype.kind!="i" or value.ndim!=1 or all(value) or not item:
            return True  # signal error
        tsa = list(item.tsa)
        tsa[1] = tuple(value.astype(int))  # convert to shape
        item.tsa = tuple(tsa)
        return False  # signal success

    def _detached_subgroup(self):  # for pdbdump (holds attributes)
        return PDBGroup(self)

    def _write_attr(self, name, value):  # for pdbdump
        if value.dtype.kind == "U":
            value = npencode(value, "utf8")  # convert to 'S'
        item = self.declare(name, value.dtype, value.shape)
        item.write(value)


class PDBLeaf(object):
    """An ndarray in a PDB file.

    (Eventual stretch goal is to implement None and zero-length arrays.)

    """
    def __init__(self, parent, addr, dtype, shape, unlim):
        root = parent.root()
        if not isinstance(dtype, tuple):
            # Construct full data type: (dtype, stype, align, typename)
            dtype = (dtype,) + root.chart.find_or_create(dtype)
        self.parent = weakref.ref(parent)
        self.attrs = None
        if shape and not all(shape):
            root.has_attributes = True  # pdbdump will write pseudo-attribute
        if hasattr(addr, '__iter__'):
            unlim = 1
            if not isinstance(addr, list):
                addr = list(addr)
        elif addr == -1:
            stype, align = dtype[1:3]
            handle = root.handle
            addr = _align(handle.next_address(), align)
            handle.declared(addr, stype, prod(shape) if shape else 1)
            if unlim:
                addr = [addr]
        elif unlim:
            addr = [int(addr)]
        else:
            addr = int(addr)
        self.tsa = dtype, shape, addr

    @staticmethod
    def isleaf():
        return 1

    @staticmethod
    def isgroup():
        return 0

    islist = isgroup

    def root(self):
        return self.parent().root()

    def query(self):
        # return dtype, shape, sshape
        dtype, shape, addr = self.tsa
        if isinstance(addr, list):
            # Do this for consistency with treatment of h5py chunked data.
            shape = (len(addr),) + shape
        return dtype[0], shape, shape

    def read(self, args=()):
        dtype, shape, addr = self.tsa
        dtype, stype, _, typename = dtype
        if typename == b'NoneType':
            return None
        if shape and not all(shape):
            return zeros(shape, dtype)
        istext = typename == b'text'
        if isinstance(addr, list):
            arg0 = args[0] if args else slice(None)
            args = args[1:]
            if not isinstance(arg0, Integral):
                arg0 = arange(len(addr))[arg0]
                if arg0.ndim == 1:
                    return array([self.read((a,) + args) for a in arg0], dtype)
                elif arg0.ndim:
                    raise TypeError("block variable leading index too complex")
            addr = addr[arg0]
        root = self.root()
        chart = root.chart
        nopartial = chart.nopartial(typename)
        if nopartial is None:
            typename = None
        if typename and nopartial:
            offset = 0
        else:
            args, shape, offset = leading_args(args, shape)
        if offset:
            addr += dtype.itemsize * offset
        f = root.handle.seek(addr)
        value = fromfile(f, stype, prod(shape) if shape else 1)
        if not nopartial:
            value = value.reshape(shape)[args]
        if typename:
            value = chart.read_special(f, typename, value)
            stype = dtype = value.dtype
        if nopartial:
            value = value.reshape(shape)[args]
        if istext and value.shape:
            return value.view('S' + str(value.shape[-1]))[..., 0]
        return value if stype is dtype else value.astype(dtype)

    def write(self, value, args=()):
        dtype, shape, addr = self.tsa
        dtype, stype, align, typename = dtype[:4]
        if dtype is None:
            if value is not None:
                raise TypeError("attempt to write value other than None to "
                                "variable declared as NoneType")
            value = array(0, stype)
        if shape and not all(shape):
            # Act out assignment to catch broadcasting errors.
            a = zeros(shape, dtype)
            a[args] = value  # If this doesn't throw an error, good to go.
            shape = ()
            value = array(0, dtype)  # dtype or stype okay here
        arg0 = args[0] if args else slice(None)
        args = args[1:]
        root = self.root()
        handle = root.handle
        if root.chart.nopartial(typename) is not None:
            raise TypeError("write to pointer type {} unsupported"
                            "".format(typename.decode('latin1')))
        if isinstance(addr, list):
            # This variable has blocks.
            if not isinstance(arg0, Integral):
                arg0 = arange(len(addr))[arg0]
                if arg0.size > 1:
                    raise TypeError("can only write block variables one "
                                    "block at a time")
                arg0 = arg0.reshape(())
            newfile = arg0 == len(addr)
            if newfile:
                # This is a new block for this variable, but not first block.
                # TODO: Should prevent partial writes here?
                selfaddr = addr
                addr, faddr = handle.next_address(both=1)
                if addr is None:
                    pass  # TODO: issue warning here and below?
                if faddr >= root.maxsize and arg0 >= root.maxblocks:
                    a = handle.next_address(newfile=1)
                    if a is not None:
                        addr = a  # Next file in family has been created.
                    else:
                        # No next filename, and current file exceeds maxsize.
                        pass  # TODO: issue warning here and above?
                addr = _align(addr, align)
                selfaddr.append(addr)
                handle.declared(addr, stype, prod(shape) if shape else 1)
            else:
                addr = addr[arg0]
        else:
            newfile = False
        args, shape, offset = leading_args(args, shape)
        if offset:
            addr += stype.itemsize * offset
        seeker = handle.seek
        f = seeker(addr)
        if args:
            # Must do read-modify-write for potentially non-contiguous write.
            v = fromfile(f, stype, prod(shape) if shape else 1).reshape(shape)
            v[args] = value
            value = v
            f = seeker(addr)
        else:
            if stype.kind == 'S' and shape:
                value = value.astype('S' + str(shape[-1]))
                value = value.reshape(value.shape + (1,)).view('S1')
            else:
                value = ascontiguousarray(value, stype)
            if value.shape != shape:
                # Avoid the recent (numpy 1.10) broadcast_to function.
                v = zeros(shape, stype)
                v[()] = value
                value = v
        value.tofile(f)

    def shifted_copy(self, delta):
        # Special helper for copying non-record variables for first file
        # to later files in a family.
        dtype, shape, addr = self.tsa
        if isinstance(addr, list):
            raise TypeError("cannot make shifted copy of record variable")
        parent = self.parent()
        if array(addr, 'u8') >> array(parent.root().handle.abits, 'u8'):
            raise TypeError("expecting non-record vars to be in first file")
        return PDBLeaf(self.parent(), addr+delta, dtype, shape, 0)


def _align(addr, align):
    if align > 1:
        rem = addr & (align - 1)
        if rem:
            addr += align - rem
    return addr


class PDBAttrs(dict):
    """Variable attributes are not a standard feature of PDB.

    We implement a poor man's version here as follows: Attributes are
    held in memory until the metadata is flushed, at which point they
    are written with name 'variable_path:attribute_name' immediately
    before the metadata.  If the file is extended, new data overwrites
    old attributes, which are rewritten just before the metadata once
    again.

    Hence, in memory, a dict suffices.

    """
    __slots__ = ()
    # qnd.QAttribute uses only __iter__, get, items, __len__, __contains__
    # PDBGroup uses __setitem__
    # Only thing that needs fixing is mapping items to iteritems for python2.
    if PY2:
        def items(self):
            return self.iteritems()
    else:
        pass
