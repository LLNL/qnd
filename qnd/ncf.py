"""QnD netCDF3 interface."""
from __future__ import absolute_import

import sys
import weakref
from collections import OrderedDict
from warnings import warn

from numpy import (dtype, prod, fromfile, asarray, array, zeros, concatenate,
                   ascontiguousarray, int64)
from numpy.core.defchararray import decode as npdecode, encode as npencode

from .frontend import QGroup
from .generic import opener
from .utils import leading_args

__all__ = ['opennc']

PY2 = sys.version_info < (3,)
if PY2:
    range = xrange  # noqa

    def itemsof(d): return d.iteritems()  # noqa
else:
    basestring = str

    def itemsof(d): return d.items()  # noqa


def opennc(filename, mode='r', auto=1, **kwargs):
    """Open netCDF-3 file returning a QnD QGroup.

    A netCDF-3 file differs from other self-describing binary file formats
    because no data addresses can be known until every variable to be
    stored is declared.  Therefore, when writing a netCDF-3 file, you
    must declare every variable before you can begin writing anything.

    The qnd API is somewhat at odds with this semantics because it encourages
    you to declare and write each variable in a single step.  The native
    netCDF-3 API forces you to declare everything, then call an `enddef`
    method to complete all variable declarations and permit you to begin
    writing data.  The qnd.ncf backend uses the first call to the ordinary
    qnd `flush` method to emulate the netCDF-3 `enddef` mode switch -- thus
    nothing will be written to the file until the first call to `flush`.
    To minimize the difference between ncf and other qnd backends, if you
    do use the usual qnd declare-and-write idiom, the ncf backend will save
    the variable value in memory until the first `flush` call, which will
    trigger the actual writing of all such saved values.

    Note that closing the file flushes it, so that is also a viable way to
    finish a netCDF-3 file.  Furthermore, when you overwrite any record
    variable in `recording` mode, ncf will implicitly `flush` the file,
    since no new variables can be declared after that.

    Note that you use the standard QnD API, a copy of every variable
    you write to the file until you begin the second record will be
    kept in memory, which could potentially be a problem.  If you wish
    to declare all variables before writing anything, so that your
    code is aligned with the netCDF API, do something like this::

       f = opennc("myfile??.nc", "w")  # wildcards expand to 00, 01, 02, ...
       # declare non-record variables from in-memory arrays
       f.nrvar1 = nrvar1.dtype, nrvar1.shape
       f.nrvar2 = nrvar2.dtype, nrvar2.shape
       # declare record variables from in-memory arrays
       f.recording(1)
       f.rvar1 = rvar1.dtype, rvar1.shape
       f.rvar2 = rvar2.dtype, rvar2.shape
       # flushing the file is equivalent to netCDF ENDDEF mode switch
       f.flush()
       # now write the current values of all the variables
       f.nrvar1 = nrvar1
       f.nrvar2 = nrvar2
       # writing the record variables writes their values for first record
       f.rvar1 = rvar1
       f.rvar2 = rvar2
       # change values of record variables and write the second record
       f.rvar1 = rvar1
       f.rvar2 = rvar2
       # when you've written all records, close the file
       f.close()

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
    **kwargs
       Other keywords.  The maxsize keyword sets the size of files in a
       family generated in recording==1 mode; a new file will begin when
       the first item in a new record would begin beyond `maxsize`.  The
       default maxsize is 128 MiB (134 MB).  The v64 keyword, if provided
       and true, causes new files to be created using the 64-bit netCDF
       format; the default is to create 32-bit files.  (But a file family
       always uses a single format.)
       The nextaddr_mode keyword can be used to indicate whether the next
       new record in 'a' or 'r+' mode should go into a new file.  The
       default behavior is that it should, which is the pdbf module default;
       this is nextaddr_mode true.  Use nextaddr_mode=0 to continue filling
       the final existing file until maxsize.

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
    v64 = kwargs.pop('v64', False)
    mode = mode.lower()
    if mode.startswith('a') or mode.startswith('r+'):
        nextaddr_mode = kwargs.pop('nextaddr_mode', 2) or 1
    else:
        nextaddr_mode = 1
    kwargs['nexaddr_mode'] = nextaddr_mode
    handle, n = opener(filename, mode, **kwargs)
    root = NCGroup(handle, maxsize, v64)
    for i in range(n):
        try:
            ncparse(handle, root, i)
        except IOError:
            # Something went terribly wrong.  If this is first file, we die.
            name = handle.filename(i)
            if not i:
                raise IOError("Fatal errors opening netCDF file {}"
                              "".format(name))
            handle.open(i-1)
            warn("file family stopped by incompatible {}".format(name))
    handle.callbacks(root.flusher, root.initializer)  # may call initializer
    return QGroup(root, auto=auto)


# https://www.unidata.ucar.edu/software/netcdf/docs/
#                              file_format_specifications.html
# All numbers are in XDR (big-endian) format.
#
# header    = magic  numrecs  dim_list  gatt_list  var_list
# magic     = 'C'  'D'  'F'  version
# version   = '\x01' (32-bit offset) | '\x02' (64-bit offset)
# numrecs   = NON_NEG | STREAMING
# dim_list  = ABSENT | 0x00 00 00 0A  NON_NEG  dim*
# gatt_list = att_list
# var_list  = ABSENT | 0x00 00 00 0B  NON_NEG  attr*
# att_list  = ABSENT | 0x00 00 00 0C  NON_NEG  var*
# ABSENT    = 0x00 00 00 00  0x00 00 00 00
# STREAMING = 0xFF FF FF FF
# dim       = name  NON_NEG   (0 length means record dimension)
# name      = NON_NEG  namestring   (0 padded to 4 byte boundary, _.@+-)
# attr      = name  nc_type  NON_NEG  values  (0 padded to 4 byte boundary)
# nc_type   = 1|2|3|4|5|6   (byte|char|short|int|float|double)
# var       = name  NON_NEG  dimid*  att_list  nc_type  vsize  OFFSET
# dimid     = 0-origin index into dim_list
# vsize     = >i4 number of bytes, or 2**32-1 if more than 4GiB
#             write vsize as if padded, but if only 1 record variable of
#             nc_type byte, char, or short, do not use padding
#             - for record variables, byte size of entire record (as if padded)
# OFFSET    = >i4 for version 1, >i8 for version 2
#
# Default fill values:
# char \x00, byte \x81, short \x80 01, int \x80 00 00 01
# float \x7C F0 00 00, double \x47 9E 00 00 00 00 00 00  =9.969209968386869e36
#
# The netCDF-3 header _almost_ has a simple XDR description; the only
# problem is that an attribute attr definition may have a value which is
# a counted array of short (2 byte integers), which XDR does not support.
# (The 64-bit format requires a hack to represent the offset values, and
# its own XDR specification using that hack.)

def ncparse(handle, root, ifile):
    i4be = _netcdf_stypes[3]
    if ifile:
        if not root.nrecs:
            raise IOError("first file in apparent family has no record vars")
        f = handle.open(ifile - 1)
        headsize = root.headsize
        f.seek(0)
        static0 = f.read(headsize)
        f = handle.open(ifile)
        magic = f.read(4)
        if magic == static0[:4]:
            nrecs = int(fromfile(f, i4be, 1)[0])
            static1 = f.read(headsize - 8)
        else:
            static1 = nrecs = None
        if static1 != static0[8:]:
            raise IOError("static variables do not match previous file")
        if nrecs == -1:
            f.seek(0, 2)
            nrecs = (f.tell() - headsize) // root.recsize
        root.nrecs.append(nrecs)
        return
    f = handle.open(ifile)
    magic = fromfile(f, 'S4', 1)[0]
    version = magic[3:]  # in python3, magic[3] is int(1) != b'\x01'
    if magic[:3] != b'CDF' or version not in b'\x01\x02':
        raise IOError("bad magic in netCDF-3 header")
    v64 = version != b'\x01'
    iobe = dtype('>i8') if v64 else i4be
    nrecs = int(fromfile(f, i4be, 1)[0])  # -1 indicates STREAMING
    tag, count = fromfile(f, i4be, 2)
    if tag != 10 and (count or tag):
        raise IOError("bad dim_list in netCDF-3 header")
    dims, recid = [], None
    while count > 0:
        count -= 1
        name = _get_name(f)
        size = int(fromfile(f, i4be, 1)[0])
        if not size:
            recid = len(dims)
        dims.append((name, size))
    attrs = [(None, _get_attrs(f))]
    tag, count = fromfile(f, i4be, 2)
    if tag != 11 and (count or tag):
        raise IOError("bad dim_list in netCDF-3 header")
    variables, recsize, special_case = OrderedDict(), 0, 0
    recaddr = lastaddr = None
    nrecvar = 0
    while count > 0:
        count -= 1
        name = _get_name(f)
        ndim = int(fromfile(f, i4be, 1)[0])
        shape = tuple(fromfile(f, i4be, ndim).astype(int)) if ndim else ()
        attrs.append((name, _get_attrs(f)))
        nctype = int(fromfile(f, i4be, 1)[0])
        if nctype < 1 or nctype > 6:
            raise IOError("bad nc_type (not in 1-6) in netCDF-3 header")
        stype = _netcdf_stypes[nctype - 1]
        fromfile(f, i4be, 1)  # ignore vsize
        offset = int(fromfile(f, iobe, 1)[0])
        # Note: offset is the byte address of the variable in the file
        #       - byte address of first block of a record variable
        if offset < 0:
            raise IOError("bad variable offset in netCDF-3 header")
        unlim = shape and shape[0] == recid
        if unlim:
            shape = shape[1:]
        try:
            sshape = tuple(dims[i][0] for i in shape)
        except IndexError:
            raise IOError("bad dimension index in netCDF-3 header")
        shape = tuple(dims[i][1] for i in shape)
        item = NCLeaf(root, len(variables), offset, stype, shape, sshape)
        variables[name] = itemx = NCList(root, item) if unlim else item
        if unlim:
            itemx.count += nrecs
            nrecvar += 1
            if nrecvar == 1:
                nbytes = stype.itemsize
                if nbytes & 3:
                    if shape:
                        nbytes *= prod(shape) if shape else 1
                    if nbytes & 3:
                        special_case = nbytes
            recsize += _measure_item(item)
            if recaddr is None or offset < recaddr:
                recaddr = offset
        elif lastaddr is None or offset >= lastaddr:
            lastaddr = offset + _measure_item(item)
    if nrecvar == 1 and special_case:
        # Implement special rule for byte, char, or short single record
        # variable; such records are not forced to 4 byte boundaries.
        recsize = special_case
    headsize = f.tell()
    if nrecs == -1 and recsize:
        # Handle special streaming record count by using file size.
        f.seek(0, 2)
        size = f.tell()
        f.seek(headsize)
        nrecs = (size - recaddr) // recsize
    root.variables = variables
    root.dims = OrderedDict(dims)
    root.attrs = OrderedDict(attrs)
    root.headsize = headsize
    root.recaddr = recaddr or lastaddr or headsize
    root.recsize = recsize
    root.nrecs.append(nrecs)
    root.v64 = v64


def _get_name(f):
    nchar = int(fromfile(f, '>i4', 1)[0])
    rem = nchar & 3
    ntot = nchar + 4 - rem if rem else nchar
    name = fromfile(f, 'S1', ntot)[:nchar].view('S' + str(nchar))
    return _bytes_as_str(name)


def _bytes_as_str(text):
    if hasattr(text, 'ravel'):
        text = text.ravel()[0]
    if isinstance(text, bytes):
        need_unicode = False
        if PY2:
            try:
                text.decode('ascii')
            except UnicodeDecodeError:
                need_unicode = True
        else:
            need_unicode = True
        if need_unicode:
            try:
                text = text.decode('utf8')
            except UnicodeDecodeError:  # ignore, but violates netCDF-3 spec
                text = text.decode('latin1')
    return text


def _text_as_bytes(text):
    if hasattr(text, 'ravel'):
        text = text.ravel()[0]
    return text if isinstance(text, bytes) else text.encode('utf8')


def _get_attrs(f):
    i4be = _netcdf_stypes[3]
    tag, count = fromfile(f, i4be, 2)
    if tag != 12 and (count or tag):
        raise IOError("bad attr_list in netCDF-3 header")
    attrs = []
    while count > 0:
        count -= 1
        name = _get_name(f)
        nctype = int(fromfile(f, i4be, 1)[0])
        if nctype < 1 or nctype > 6:
            raise IOError("bad nc_type (not in 1-6) in netCDF-3 header")
        if nctype == 2:
            values = _get_name(f)
        else:
            nvalues = int(fromfile(f, i4be, 1)[0])
            stype = _netcdf_stypes[nctype - 1]
            values = fromfile(f, stype, nvalues)
            rem = values.nbytes & 3
            if rem:
                fromfile(f, 'u1', 4 - rem)
            if values.size == 1:
                values = values[0]
            if not stype.isnative:
                values = values.astype(stype.newbyteorder('='))
        attrs.append((name, values))
    return OrderedDict(attrs)


class NCGroup(object):
    def __init__(self, handle, maxsize=134217728, v64=False):
        self.handle = handle  # a generic.MultiFile
        self.variables, self.dims, self.attrs = {}, {}, {}
        self.headsize = self.recaddr = self.recsize = 0
        self.nrecs = []  # list of record counts in files of family
        self.maxsize = maxsize
        self.v64 = v64
        self.pending = None  # holds pre-flush variable values

    @staticmethod
    def isgroup():
        return 1

    @staticmethod
    def islist():
        return 0

    isleaf = islist

    def root(self):
        return self  # no such thing as directories in netCDF3

    def close(self):
        self.handle.close()

    def flush(self):
        self.handle.flush()

    def __len__(self):
        return len(self.variables)

    def __iter__(self):
        return iter(self.variables)

    def lookup(self, name):
        return self.variables.get(name)

    def declare(self, name, dtype, shape, unlim=None):
        if self.headsize:
            raise RuntimeError("netCDF file defined, no more declarations")
        if shape and not all(shape):
            raise TypeError("netCDF does not support 0-length dimensions")
        stype = _get_stype(dtype)
        sshape = tuple('_' + str(s) for s in shape) if shape else ()
        dims, variables = self.dims, self.variables
        if unlim:
            dims.setdefault('_0', 0)
        for s, n in zip(sshape, shape):
            dims.setdefault(s, n)
        # Set offset to unlim for now, will be set in initializer.
        item = NCLeaf(self, len(variables), unlim, stype, shape, sshape)
        if unlim:
            item = NCList(self, item)
        variables[name] = item
        return item

    # qnd.QAttribute uses only __iter__, get, items, __len__, __contains__
    # In PY2, the dict returned here has an inefficient items() method,
    # but it is not worth fixing that here.

    def attget(self, vname):
        return self.attrs.get(vname if vname else None)

    def attset(self, vname, aname, dtype, shape, value):
        if self.headsize:
            raise RuntimeError("netCDF file defined, no setting attributes")
        stype = _get_stype(dtype)
        strtype = _netcdf_stypes[1]
        if stype == strtype:
            if shape:
                raise TypeError("netCDF does not support array of strings"
                                "as an attribute value")
            value = _bytes_as_str(value)
        else:
            value = asarray(value, stype)
            if shape:
                if len(shape) > 1:
                    raise TypeError("netCDF does not support "
                                    "multi-dimensional attribute values")
                if value.shape != shape:
                    value = value.reshape(shape)
            if not stype.isnative:
                value = value.astype(stype.newbyteorder('='))
        if not vname:
            vname = None
        attrs = self.attrs.get(vname)
        if not attrs:
            self.attrs[vname] = attrs = OrderedDict()
        attrs[aname] = value

    def record_delta(self, irec):
        """Compute delta to add to record variable offset to reach irec."""
        handle, nrecs, maxsize = self.handle, self.nrecs, self.maxsize
        if not self.headsize:
            if not nrecs:  # This is first record variable.
                nrecs.append(1)
            if not irec:
                return 0  # First record is being declared, delta unknown.
            # Beginning to write second record may force first flush,
            # freezing the netCDF file structure, an implicit ENDDEF.
            # Effectively, this flush is writing the first record, even
            # though this call has irec==1 and declaring the first variable
            # of the second record.
            self.flusher(self.handle.open(0))
        recsize = self.recsize
        rec0 = array(nrecs).cumsum()
        # searchsorted needs strictly monotonic array
        # However, because of the 0.5 offset and the fact that irec is
        # an integer, this apparently can never cause a problem here
        # (the monotonicity problem only arises if irec matches two
        # consecutive equal values of rec0-0.5, which could happen if
        # some file has no records).
        ifile = (rec0 - 0.5).searchsorted(irec)
        if ifile >= rec0.size:
            if handle.nextaddr:
                # Handle special case of the first record written after a
                # family is opened in 'a' or 'r+' mode.
                maxsize = 0
            # This is a new record.  We check if maxsize has been exceeded,
            # and force a new file in the family to be created if so.
            n = nrecs[-1]
            if n and (self.recaddr + recsize*n >= maxsize):
                f = handle.open(ifile)  # implicit flush during open
                nrecs.append(0)
                self.initializer(f)
                irec -= rec0[-1]
            else:
                ifile -= 1  # add record to last existing file
                if ifile:
                    irec -= rec0[ifile - 1]
            handle.nextaddr = int64(0)  # special case only triggers once
            nrecs[-1] += 1
        elif ifile:
            irec -= rec0[ifile - 1]
        return handle.zero_address(ifile) + recsize * irec

    def flusher(self, f):
        # The flush method has to serve as ENDDEF for newly created netCDF
        # families, see comments for initializer method below.
        if not self.headsize:
            # Only get here for first file of newly created family.
            self.headsize = 1  # impossible since magic number is 4 bytes
            self.initializer(f)
        # The only metadata that may need to be written is nrecs.
        # The file handle f is the last file in the family.
        if self.nrecs:
            f.seek(4)
            array(self.nrecs[-1], '>i4').tofile(f)

    def initializer(self, f):
        # Called indirectly by handle.callbacks during ncopen in "w" mode.
        # This is the only case in which this point is reachable with zero
        # self.headsize, because ncparse would have filled it in in "r" or
        # "r+" mode, and it would have been written for the first file in
        # the family if this is not the first file.
        # For the first file of a newly created netCDF family, we want to
        # wait for an explicit call to flush() to write the first file header,
        # which is the QnD implementation of the ENDDEF call in the netCDF API.
        first_flush = self.headsize == 1  # impossible value set in flusher
        if first_flush:
            self.headsize = 0
        else:
            return
        # The file f positioned at address 0.
        i4be = _netcdf_stypes[3]
        v64 = self.v64
        array(b'CDF' + (b'\x02' if v64 else b'\x01')).tofile(f)
        array(0, i4be).tofile(f)
        handle = self.handle
        ifile = handle.current_file()
        if ifile:
            # Just copy header and non-record variables to new file.
            f = handle.open(0)
            f.seek(8)
            value = f.read(self.recaddr)
            f = handle.open(ifile)
            f.seek(8)
            f.write(value)
            return

        # This is first file of family.
        dims, variables, attrs = self.dims, self.variables, self.attrs
        if not dims:
            zeros(2, i4be).tofile(f)
        else:
            array((10, len(dims)), i4be).tofile(f)
            for name, size in itemsof(dims):
                _put_name(f, name)
                array(size, i4be).tofile(f)
        _put_attrs(f, attrs.get(None))
        if not variables:
            zeros(2, i4be).tofile(f)
        else:
            array((11, len(variables)), i4be).tofile(f)
        headsize = f.tell()  # including vars tag and count
        if first_flush:
            # The offsets in the variables array are unknown until the
            # symbol table is written, which makes it hard to write for
            # the first file in a family.  We make a clumsy two passes
            # to compute the length of the var_list if the offsets have
            # not yet been set.
            # add space for name length, ndim, nctype, vsize, and offset
            headsize += (20 + 4*v64) * len(variables)
            nrecs = self.nrecs
            for name, item in itemsof(variables):
                unlim = isinstance(item, NCList)
                if unlim:
                    item = item.leaf
                    if not nrecs:
                        nrecs.append(0)
                ndim = len(item.shape or ()) + unlim  # so shape None okay
                vattrs = attrs.get(name)
                namelen = _put_name(None, name)
                headsize += namelen + 4*ndim + _measure_attrs(vattrs)
            offset = self.headsize = headsize
            # Now we can fill in all the offsets and find recaddr.
            recitems = []
            for name, item in itemsof(variables):
                if isinstance(item, NCList):
                    item = item.leaf
                if item.offset:  # This is unlim, see NCGroup.declare.
                    recitems.append(item)
                    continue
                item.offset = offset
                offset += _measure_item(item)
            self.recaddr = offset
            for item in recitems:
                item.offset = offset
                offset += _measure_item(item)
            self.recsize = offset - self.recaddr
        recaddr, recsize = self.recaddr, self.recsize
        dimids = {name: i for i, name in enumerate(dims)}
        recid = None
        for i, (_, n) in enumerate(itemsof(dims)):
            if not n:
                recid = i
                break
        recid = [] if recid is None else [recid]
        iobe = dtype('>i8') if v64 else i4be
        rem = recsize & 3
        if rem:
            recsize += 4 - rem  # used only for vsize
        if recsize > 0xffffffff:
            recsize = 0xffffffff  # vsize overflow convention
        for name, item in itemsof(variables):
            if isinstance(item, NCList):
                item = item.leaf
            stype, offset = item.stype, item.offset
            nctype = _netcdf_stypes.index(stype) + 1
            sshape = item.sshape or ()
            unlim = offset >= recaddr
            _put_name(f, name)
            array(len(sshape) + unlim, i4be).tofile(f)
            sshape = (recid if unlim else []) + [dimids[s] for s in sshape]
            array(sshape, i4be).tofile(f)
            _put_attrs(f, attrs.get(name))
            vsize = recsize if unlim else _measure_item(item)
            array([nctype, vsize], i4be).tofile(f)
            array(offset, iobe).tofile(f)
        headsize = f.tell()
        if headsize != self.headsize:
            raise IOError("netCDF header size mismatch (BUG?)")
        # Header finished, write any pending variables now.
        pending = self.pending
        self.pending = None
        if pending:
            byindex = {}
            for _, item in itemsof(variables):
                if isinstance(item, NCList):
                    item = item.leaf
                byindex[item.index] = item
            for index, value in itemsof(pending):
                byindex[index].write(value)


def _put_name(f, name):
    name = _text_as_bytes(name)
    nchar = len(name)
    rem = nchar & 3
    if f is None:
        rem = (4 - rem) if rem else 0
        return nchar + rem  # not including 4 byte nchar count
    if rem:
        name = name + b'\0'*(4 - rem)
    array(nchar, _netcdf_stypes[3]).tofile(f)
    f.write(name)
    return None


def _put_attrs(f, attrs):
    i4be = _netcdf_stypes[3]
    if not attrs:
        zeros(2, i4be).tofile(f)
        return
    array((12, len(attrs)), i4be).tofile(f)
    for name, value in itemsof(attrs):
        if isinstance(value, basestring):
            nctype = 2
            value = _text_as_bytes(value)
            n = len(value)
            rem = n & 3
            if rem:
                value += b'\0' * (4 - rem)
            value = array(value)
        else:
            value = value.asarray(value)
            dtype = value.dtype
            size = dtype.itemsize
            if dtype.kind == 'f':
                nctype = 5 + (size == 8)
            elif size == 1:
                nctype = 1
            else:
                nctype = 3 + (size == 4)
            stype = _netcdf_stypes[nctype - 1]
            if dtype != stype:
                value = value.astype(stype)
            n = value.size
            if nctype == 3 and (value.size & 1):
                value = concatenate((value.ravel(), zeros(1, stype)))
        _put_name(f, name)
        array((nctype, n), i4be).tofile(f)
        value.tofile(f)


def _measure_attrs(attrs):
    size = 8
    if attrs:
        for name, value in itemsof(attrs):
            size += 24  # name length, nctype, value count
            size += ((len(_text_as_bytes(name)) + 3) >> 2) << 2
            if isinstance(value, basestring):
                size += len(_text_as_bytes(value))
            else:
                size += value.asarray(value).nbytes
            size = ((size + 3) >> 2) << 2
    return size


def _measure_item(item):
    size = item.shape
    size = prod(size) if size else 1
    nbytes = item.stype.itemsize * size
    return ((nbytes + 3) >> 2) << 2


def _get_stype(dtype):
    # bewawre misfeature numpy (1.16.4) dtype('f8') tests == None
    kind = 'X' if dtype is None or dtype in (dict, list,
                                             object) else dtype.kind
    stype = None
    if kind in 'bui':
        size = dtype.itemsize
        sizes = (1, 2, 4, 8)
        if size in sizes:
            stype = _netcdf_stypes[(0, 2, 3, 3)[sizes.index(size)]]
    elif kind == 'f':
        size = dtype.itemsize
        sizes = (2, 4, 8, 12, 16)
        if size in sizes:
            stype = _netcdf_stypes[(4, 4, 5, 5, 5)[sizes.index(size)]]
    elif kind in 'SU':
        stype = _netcdf_stypes[1]
    if stype is None:
        raise TypeError("netCDF-3 does not support this dtype")
    return stype


_netcdf_stypes = [dtype('i1'), dtype('S1'), dtype('>i2'), dtype('>i4'),
                  dtype('>f4'), dtype('>f8')]


class NCLeaf(object):
    __slots__ = 'parent', 'index', 'offset', 'stype', 'shape', 'sshape'

    def __init__(self, parent, index, offset, stype, shape, sshape, _wrp=None):
        self.parent = parent if _wrp else weakref.ref(parent)
        self.index = index
        self.offset = offset
        self.stype = stype
        self.shape = shape
        self.sshape = sshape

    @staticmethod
    def isleaf():
        return 1

    @staticmethod
    def isgroup():
        return 0

    islist = isgroup

    def shift_by(self, delta):
        state = [getattr(self, nm) for nm in self.__slots__]
        state[2] += delta
        return NCLeaf(*state, _wrp=1)

    def root(self):
        return self.parent()

    def _dtype(self):
        dtype = self.stype
        return dtype if dtype.isnative else dtype.newbyteorder('=')

    def query(self):
        # return dtype, shape, sshape
        shape, sshape = self.shape or (), self.sshape
        return self._dtype(), shape, sshape if sshape else shape

    def read(self, args=()):
        parent = self.parent()
        if not parent.headsize:
            raise RuntimeError("cannot read from netCDF file in 'w' mode"
                               " before first flush")
        stype, shape = self.stype, self.shape
        args, shape, offset = leading_args(args, shape)
        f = parent.handle.seek(self.offset + stype.itemsize * offset)
        size = prod(shape) if shape else 1
        value = fromfile(f, stype, size).reshape(shape)[args]
        if not stype.isnative:
            value = value.astype(stype.newbyteorder('='))
        if stype == _netcdf_stypes[1]:
            # Present this as a str or array of str.
            # Note that final netCDF dimension is really length of string.
            shape = value.shape
            if shape:
                shape, strlen = shape[:-1], shape[-1]
                value = value.view('S' + str(strlen)).reshape(shape)
            if PY2:
                try:
                    npdecode(value, 'ascii')
                    need_unicode = False
                except UnicodeDecodeError:
                    need_unicode = True
            else:
                need_unicode = True
            if need_unicode:
                try:
                    value = npdecode(value, 'utf8')
                except UnicodeDecodeError:
                    value = npdecode(value, 'latin1')
            if not shape:
                value = value[()]
        return value

    def write(self, value, args=()):
        parent = self.parent()
        if not parent.headsize:
            if args:
                raise IndexError("no partial writes during declaration")
            pending = parent.pending
            if pending is None:
                pending = parent.pending = {}
            pending[self.index] = value
            return
        offset, stype, shape = self.offset, self.stype, self.shape
        args, shape, off = leading_args(args, shape)
        if off:
            offset += stype.itemsize * off
        value = asarray(value)
        kind = value.dtype.kind
        if kind in 'SU':
            if kind == 'U':
                value = npencode(value, 'utf8')
            shape = value.shape
            value = value.reshape(shape + (1,)).view('S1')
        f = parent.handle.seek(offset)
        if args:
            # Must do read-modify-write for potentially non-contiguous write.
            addr = f.tell()
            v = fromfile(f, stype, prod(shape) if shape else 1).reshape(shape)
            v[args] = value
            value = v
            f.seek(addr)
        else:
            value = ascontiguousarray(value, stype)
            if value.shape != shape:
                # Avoid the recent (numpy 1.10) broadcast_to function.
                v = zeros(shape, stype)
                v[()] = value
                value = v
        value.tofile(f)


class NCList(object):
    """NCLeaf wrapper for record variables."""
    __slots__ = 'parent', 'leaf', 'count'

    def __init__(self, parent, leaf):
        self.parent = weakref.ref(parent)
        self.leaf = leaf
        self.count = 0  # record count needed to know when new record created

    @staticmethod
    def islist():
        return 1

    @staticmethod
    def isgroup():
        return 0

    isleaf = isgroup

    def root(self):
        return self.parent()

    # len, iter, index, declare are list methods called by QList

    def __len__(self):
        return sum(self.parent().nrecs)

    def __iter__(self):
        for i in range(len(self)):
            yield self.index(i)

    def index(self, ndx):
        nrecs = len(self)
        if ndx < 0:
            ndx = ndx + nrecs
        if ndx < 0 or ndx >= nrecs:
            return None  # out of range, let caller raise any exception
        parent = self.parent()
        delta = parent.record_delta(ndx)
        return self.leaf.shift_by(delta)

    def declare(self, dtype, shape):
        # Ignore dtype and shape here; conformability with the NCLeaf
        # dtype and shape will be enforced during NCLeaf.write.
        parent = self.parent()
        delta = parent.record_delta(self.count)
        self.count += 1  # nrecs in NCGroup incremented in record_delta
        return self.leaf.shift_by(delta)
