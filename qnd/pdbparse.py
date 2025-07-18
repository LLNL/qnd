"""Pure python PDB file format parsing.

Parse PDB metadata.  The PDB file format was designed by Stewart Brown
at Lawrence Livermore National Laboratory in the 1990s.

PDB is widely used at LLNL for restart and post-processing data
produced by large radiation-hydrodynamics simulation codes.  The PDB
format can describe and store named arrays of any data type
representable in the C programming language that is derived from one
of the primitive C data types char, short, int, long, float, or double
or pointer to a representable type.  The format was extended to handle
primitive integer or floating point numbers of any size, such as long
long or long double or 16 bit floating point data.  The format was
also extended to organize the named arrays into groups like HDF5.  The
metadata is text with a few embedded ASCII control characters, written
at the end of the file, and intended to be parsed and held in memory
when the file is opened.

"""
from __future__ import absolute_import

import re
import sys
from collections import OrderedDict
from warnings import warn
import weakref

from numpy import (prod, array, arange, concatenate, zeros, dtype as npdtype,
                   fromfile, int64)

PY2 = sys.version_info < (3,)
if PY2:
    from itertools import izip as zip
    range = xrange

    def itemsof(d): return d.iteritems()  # noqa
else:
    from numpy.core.defchararray import decode as npdecode  # noqa
    basestring = str

    def itemsof(d): return d.items()  # noqa


class PDBChart(object):
    """PDB structure chart contains all information about data types."""
    # stype --> dtype on disk
    # All typename and membername values are bytes, not unicode (str).
    # The assumed codec is latin1; legacy PDB files could support UTF-8,
    # but without any explicit mandate, latin1 is the safer option for
    # reading files of unknown origin, since unlike UTF-8 there are no
    # undefined bytes.
    #
    # hasint64   True if any int primitive is 8 bytes, else False.
    #            False means that i8 dtypes will be mapped to i4.
    # byteorder  Determines new stype corresponding to new primitive dtype.
    # structal   Initial struct alignment, before adding any members.
    # csaddr     Chart and symtab address in header.
    # haspointers  Pointer types used: 0 none, 1 PDB style, 2 yorick style
    # primitives [name] --> stype, dtype, align, fpbits or None
    #            fpbits present only if stype is Vn for non-IEEE float.
    # structs    [name] --> stype, dtype, align, members
    #            members is OrderedDict membername --> typename, shape
    # by_dtype   [dtype] --> typename, table (either primitives or structs)
    #            This maps dtype to a PDB typename for interpreting newly
    #            declared variables.  The names _anon_0, _anon_1, ... are
    #            created for dtypes which are previously undeclared.
    # nanons     Count of _anon_X typenames already used.
    # pointers   [typename] --> otype corresponding to stype, dtype
    # yortypes   list of primitives+structs typenames
    #
    # Strategy for converting pointers:
    # 1. Convert stype --> dtype as usual; yorick can do partial reads
    # 2. Create final otype = dtype('O') variant of dtype
    # 3. Read the stype array.
    # 4. Create empty otype array and assign stype array.  In the case of
    #    yorick pointers, numpy converts the pointer values to type 'O'.
    # 5. Iterate through the otype one array item, then one struct item
    #    at a time, converting the 'O' members from pointer to pointee.
    #    Recurse during this loop to handle any pointers encountered in
    #    the pointees.

    def __init__(self, root):
        self.root = weakref.ref(root)  # root() is top level PDBGroup
        self.hasint64 = False
        self.structal = self.byteorder = self.csaddr = None
        self.haspointers = 0
        self.primitives = OrderedDict()
        self.structs = OrderedDict()
        self.by_dtype = {}
        self.nanons = 0
        self.pointers = {}
        self.yortypes = []

    _orders = {1: '>', 2: '<'}
    # Create distinguishable pointer and string dtypes.
    _dtypeo = npdtype('O')
    _dtypes = npdtype('O', copy=1)

    def find_or_create(self, dtype):
        """Find or create stype --> (stype, align, typename)."""
        dtype = self._saveable(dtype)
        typename, table = self.by_dtype.get(dtype, (None, None))
        if typename is not None:
            stype, _, align, _ = table[typename]
        else:
            stype, align, typename = self._add_dtype(dtype)
        return stype, align, typename

    def nopartial(self, typename):
        if self.pointers.get(typename) is not None:
            return bool(self.haspointers == 1)
        return None

    def use(self, dtype):
        """Find or create stype --> (stype, align, typename, nopartial)."""
        stype, align, typename = self.find_or_create(dtype)
        # The native PDB pointer types require special treatment; the
        # entire variable and all of its indirect references must be
        # read at once (no partial reads).
        # The yorick string and pointer types require special post-processing
        # after a read in order to interpret the pointer values.
        # Any struct which has a member requiring special treatment itself
        # requires special treatment.
        return stype, align, typename, self.nopartial()

    def _saveable(self, dtype):
        if dtype is None:
            return dtype
        if dtype.shape:
            raise TypeError("cannot store subarray dtype in PDB file")
        if not dtype.isnative:
            dtype = dtype.newbyteorder('=')
        if dtype.kind in 'ui' and dtype.itemsize > 4:
            # If file does not have int64 type, map int64 to int32.
            if not self.hasint64:
                dtype = npdtype(dtype.kind + '4')
        return dtype

    def _add_dtype(self, dtype):
        if dtype is None:
            stype = self.add_primitive(b'NoneType', (1, 0, 0))
            return stype, 0, b'NoneType'
        by_dtype, nonenone = self.by_dtype, (None, None)
        # Need to declare anonymous type now.
        size = dtype.itemsize
        if size == 1 and dtype.kind == 'S':
            # Declare QnD-specific text primitive.
            stype = self.add_primitive(b'text', (1, 0, 0))
            return stype, 0, b'text'
        if dtype.kind == 'c':
            fsize = size // 2
            ffmt = npdtype('f{}'.format(fsize))
            fields = dict(re=(ffmt, 0), im=(ffmt, fsize))
            names = ['re', 'im']
            if fsize == 8:
                name = b'complex'
            elif fsize == 4:
                name = b'fcomplex'
            else:
                name = 'complex{}'.format(size*8).encode('latin1')
        else:
            fields = dtype.fields
            names = dtype.names
            n = self.nanons
            self.nanons += 1
            name = '_anon_{}'.format(n).encode('latin1')
        if not names:
            order = self.byteorder
            if not order:
                order = dtype.str[0]
                order = 1 if order == '>' else (2 if order == '<' else 0)
            align = size  # WAG, could create align=1 record dtype to check
            desc = size, 0 if size == 1 else order, align
            adder = self.add_primitive
            if dtype.kind == 'f':
                if size == 4:
                    desc += (_binary32,)
                elif size == 8:
                    desc += (_binary64,)
                else:
                    raise NotImplementedError(
                        "only f4 and f8 floating point dtypes supported")
        else:
            align = self.structal
            saveable = self._saveable
            desc = OrderedDict()  # really members list
            _dtype_shape = self._dtype_shape
            for nm in names:
                dtyp, shape = _dtype_shape(fields[nm])
                dtyp = saveable(dtyp)
                typename, table = by_dtype.get(dtyp, nonenone)
                if typename is None:
                    # Even though we are just building the members list here,
                    # we need to recursively create any member stypes which
                    # have not yet been encountered.  Note that the recursion
                    # produces the anonymous types in order of dependency.
                    _, al, typename = self._add_dtype(dtyp)
                else:
                    al = table[typename][2]
                if al > align:
                    align = al
                desc[nm.encode('latin1')] = typename, shape
            desc = size, desc
            adder = self.add_struct
        stype = adder(name, desc)
        return stype, align, name

    @staticmethod
    def _dtype_shape(field):
        dtype = field[0]
        dtype, shape = dtype.base, dtype.shape
        s = dtype.shape
        while s:  # numpy does not automatically condense subarrays
            shape += s
            dtype, s = dtype.base, dtype.shape
        return dtype, shape

    def add_primitive(self, name, desc):
        """desc is (size, order, align) or (size, order, align, fpbits)"""
        ispdbptr = name == b'*'
        if ispdbptr:
            name = b'char *'
        size, order, align = desc[:3]
        if not order:
            order = 0
        elif not hasattr(order, '__getitem__'):
            order = int(order)
        else:
            perm = array(order)
            if perm[0] == 1:
                order = 1
            elif perm[-1] == 1:
                order = 2
                perm = perm[::-1]
            if perm.size != size or (arange(1, 1+size) != perm).any():
                order = 0
        fpbits = desc[3] if len(desc) > 3 else None
        if align & (align - 1):
            # Error if alignment is not a power of 2.
            raise IOError("alignment {} for primitive {} not power of 2"
                          "".format(align, name.decode('latin1')))
        order = self._orders.get(order, '')
        if not order and size > 1:
            stype = npdtype('V{}'.format(size))
        elif fpbits:
            if size == 4 and fpbits == _binary32:
                stype, fpbits = npdtype(order + 'f4'), None
            elif size == 8 and fpbits == _binary64:
                stype, fpbits = npdtype(order + 'f8'), None
            else:
                stype = npdtype('V{}'.format(size))
        else:
            lname = name.lower()
            if lname.startswith(b'u') or name == b'char':
                kind = order + 'u{}'
            elif size == 1:
                if lname in (b'bool', b'boolean'):
                    kind = order + 'b{}'
                elif name == b'text':
                    kind = 'S{}'
                else:
                    kind = 'i{}'
            else:
                kind = order + 'i{}'
            if size in (1, 2, 4, 8):
                stype = npdtype(kind.format(size))
                if size == 8:
                    self.hasint64 = True
            else:
                stype = npdtype('V{}'.format(size))
        if stype.isnative:
            dtype = stype
        else:
            dtype = stype.newbyteorder('=')
        primitives = self.primitives
        current = primitives.get(name)
        if current:
            if current == (stype, dtype, align, fpbits):
                return current[0]
            else:
                return None
        primitives[name] = stype, dtype, align, fpbits
        self.by_dtype.setdefault(dtype, (name, primitives))
        if ispdbptr:
            dtype = npdtype('V0')
            primitives[b'char*'] = dtype, dtype, 0, None
            pointers = self.pointers
            pointers[b'char *'] = pointers[b'char*'] = self._dtypeo
        elif self.haspointers == 2 and name == b'double':
            slong, dlong, along, _ = primitives[b'long']
            primitives[b'string'] = slong, dlong, along, None
            primitives[b'pointer'] = slong, dlong, along, None
            pointers = self.pointers
            pointers[b'string'] = self._dtypes
            pointers[b'pointer'] = self._dtypeo
        return stype

    def add_struct(self, name, members):
        """members is OrderedDict, name --> typename, shape"""
        if not isinstance(name, bytes):
            name = name.encode('latin1')
        primitives, structs = self.primitives, self.structs
        align = self.structal
        pointers = self.pointers
        size, members = members
        names, formats, oformats = [], [], []
        offsets, off, special = [], 0, False
        for mname, (tname, shape) in itemsof(members):
            # mname must be bytes, str_mname must be str
            # hash(dtype) fails if mname is unicode in PY2
            if isinstance(mname, bytes):
                str_mname = mname if PY2 else mname.decode('latin1')
            elif PY2:
                str_mname = mname = mname.encode('latin1')
            else:
                str_mname, mname = mname, mname.encode('latin1')
            if not isinstance(tname, bytes):
                tname = tname.encode('latin1')
            if tname.endswith(b'*'):
                tname = b'char *'
            stype = primitives.get(tname)
            isstruct = stype is None
            if isstruct:
                stype = structs.get(tname)
                if stype is None:
                    raise IOError("struct {} refers to undefined type {}"
                                  "".format(name.decode('latin1'),
                                            tname.decode('latin1')))
            stype, dtype, al, _ = stype
            if stype.shape:
                shape = stype.shape + shape
                stype = stype.base
            if al > 1:
                rem = off & (al - 1)
                if rem:
                    off += al - rem
            if al > align:
                align = al
            offsets.append(off)
            formats.append(npdtype((stype, shape)) if shape else stype)
            names.append(str_mname)
            off += stype.itemsize * (prod(shape) if shape else 1)
            ptrtype = pointers.get(tname)
            if ptrtype is not None:
                dtype, special = ptrtype, True
            oformats.append(npdtype((dtype, shape)) if shape else dtype)
        # if align > 1:
        #     rem = off & (align - 1)
        #     if rem:
        #         off += align - rem
        lname = name.lower()
        iscomplex = b'complex' in lname or b'plx' in lname and len(names) == 2
        if iscomplex and formats[0] == formats[1] and formats[0].kind == 'f':
            re, im = [n.lower().replace('part', '').replace('_', '')
                      for n in names]
            iscomplex = (re, im) in _complex_members
        if iscomplex:
            stype = npdtype(formats[0].byteorder + 'c{}'.format(size))
        else:
            stype = npdtype(dict(names=names, formats=formats,
                                 offsets=offsets, itemsize=size))
        dtype = stype if stype.isnative else stype.newbyteorder('=')
        current = structs.get(name)
        if current:
            if current == (stype, dtype, align, members):
                return current[0]
            else:
                return None  # Attempt to redefine struct type.
        if special:
            pointers[name] = npdtype(dict(names=names, formats=oformats))
        structs[name] = stype, dtype, align, members
        self.by_dtype.setdefault(dtype, (name, structs))
        return stype

    def read_special(self, f, typename, value):
        """Recursively read pointer values."""
        # Convert value from as-stored stype to final caller dtype.
        # This will fill all the dtype('O') fields corresponding to pointers
        # with whatever was stored in the corresponding slot in the stype.
        # (In the case of either yorick or PDB itag pointers, this is
        # a scalar int value.)
        # value can be non-ndarray in v[mname][0] call a few lines down:
        shape = value.shape if hasattr(value, 'shape') else ()
        v = zeros(shape, self.pointers[typename])
        v[()] = value  # works even for PDB-II char* 0-length void value
        value = v
        by_address = {}
        dtype = value.dtype
        names = dtype.names
        # Iterate through each individual pointer in the order which the
        # native PDB pointer mechanism would process them.  This isn't
        # terribly efficient, but reading the individual pointees is
        # likely to be even more expensive than all these nested lookups.
        if names:
            read_special = self.read_special
            pointers = self.pointers
            members = self.structs[typename][3]
            for v in value.reshape(value.size, 1):
                for mname, (tname, shape) in itemsof(members):
                    if tname in pointers:
                        mname = mname.decode('latin1')
                        v[mname][0] = read_special(f, tname, v[mname][0])
                        # If the member has record dtype, it will not become
                        # a recarray here.  This is consistent with the
                        # behavior of non-pointer members.
                        # An additional problem is that the numpy interface
                        # does not handle O scalars as you might wish - when
                        # you extract such a member it remains a scalar
                        # object array, rather than the object itself.
        else:
            haspointers = self.haspointers
            reader = (self.read_yorickptr if haspointers == 2 else
                      self.read_pdbptr)
            for v in value.reshape(value.size, 1):
                v[0] = reader(f, v, by_address)
            if dtype is self._dtypes:
                # Convert from dtype O to string or dtype S
                value = array(value.tolist())
                if not PY2:
                    try:
                        value = npdecode(value, 'utf8')
                    except UnicodeDecodeError:
                        value = npdecode(value, 'latin1')
            if not value.shape:
                value = value[()]
        return value

    def read_pdbptr(self, f, value, by_address):
        addr0 = f.tell()
        header = block = f.read(1024)
        while True:
            m = _line_re.search(block)
            if m:
                break
            header += block
            if len(block) < 1024:
                break
            f.read(1024)
        nhead = m.start if m else len(header) - 1
        addr1 = addr0 + nhead + 1
        f.seek(addr1)  # position file to byte after header
        m = _ptrheader.match(header)
        if m:
            nitems, typename, addr, here = m.group(1, 2, 3, 4)
            addr = int64(-1) if addr is None else int64(addr)
            if addr == -1:
                return None  # This is a NULL pointer.
            if here is not None and not int64(here):
                # We need to get this pointee from elsewhere.
                v = by_address.get(addr, Ellipsis)
                if v is not Ellipsis:
                    return v  # already read and cached this pointee
                # Read pointee from elsewhere, then snap back to here.
                f.seek(addr)
                v = self.read_pdbptr(f, value, by_address)
                f.seek(addr1)
                return v
            # Data is here at addr0, read it now.
            # Unknown what happens if addr != addr0 -- if actual header
            # address isn't used when data not here, this algorithm will
            # fail (e.g.- if addr for one not here points to another not
            # here in a chain).
            nitems = int64(nitems)
            typename = typename.strip()
            stype = self.primitives[typename]
            if stype is None:
                stype = self.structs[typename]
            if stype is not None:
                stype, dtype = stype[:2]
                v = fromfile(f, stype, nitems)
                if typename in self.pointers:
                    # Recurse because pointee contains pointers.
                    v = self.read_special(f, typename, v)
                elif stype is not dtype:
                    v = v.astype(dtype)
                if nitems == 1:
                    # Presume intent is a scalar, PDB offers no way to find
                    # the original shape of a pointee.
                    v = v.reshape(())
                by_address[addr0] = v
                # File is positioned immediately after block of pointees
                # which has just been read.
                return v
        raise IOError("Unable to track PDB pointees at {}".format(addr0))

    def read_yorickptr(self, f, value, by_address):
        addr = int64(value[0])
        v = by_address.get(addr, Ellipsis)
        if v is not Ellipsis:
            return v
        yortypes, primitives = self.yortypes, self.primitives
        if not yortypes:
            # Create yorick type list now.
            yortypes = list(primitives) + list(self.structs)
            self.yortypes = yortypes
        null = addr < 0
        if not null:
            f.seek(addr)
        longstype = primitives[b'long'][0]
        if value.dtype is self._dtypes:
            # This is a yorick string.
            n = 0 if null else fromfile(f, longstype, 1).astype(int64)[0]
            return fromfile(f, npdtype('S{}'.format(n)), 1)[0] if n else b''
        # This is a general yorick pointee.
        if null:
            return None
        ytype, ndim = fromfile(f, longstype, 2).astype(int64)
        if ytype < 0 or ytype >= len(yortypes) or ndim > 10:
            IOError("bad yorick pointee header at {}".format(addr))
        if ndim:
            shape = tuple(fromfile(f, longstype, ndim).astype(int64)[::-1])
            size = prod(shape)
        else:
            shape = ()
            size = 1
        typename = yortypes[ytype]
        if ytype < len(primitives):
            stype, dtype, align = primitives[typename][:3]
        else:
            stype, dtype, align = self.structs[typename][:3]
        addr += (2 + ndim) * longstype.itemsize
        if align > 1:
            rem = addr & (align - 1)
            if rem:
                addr += align - rem
                f.seek(addr)
        v = fromfile(f, stype, size).reshape(shape)
        if typename in self.pointers:
            # Recurse because pointee contains more pointers.
            v = self.read_special(f, typename, v)
        elif stype is not dtype:
            v = v.astype(dtype)
        return v


_complex_members = (('r', 'i'), ('re', 'im'), ('real', 'imag'),
                    ('real', 'imaginary'))


def parser(handle, root, index=0, allrecs=None, atnames=None, hooks=None):
    """Parse PDB file with the given MultiFile handle and PDBGroup group."""
    f = handle.open(index)

    # Begin by reading header (or trailer for version 3).
    header = f.read(1024)
    m = _magic2.match(header)
    if m:
        # Most common version 2 PDB format.
        root.pdb_version = version = 2
        n = 13 + (ord(header[13]) if PY2 else header[13])
        priminfo = map(ord, header[14:n]) if PY2 else [h for h in header[14:n]]
        ma = _magic2a.match(header[n:])
        if ma:  # save address of chart and symtab addresses for pdbdump
            root.chart.csaddr = n + ma.end(2) + 2
        try:
            nums = list(map(int, ma.group(1, 2, 3, 4)))
        except (AttributeError, ValueError):
            nums = [-1]
        if any(n < 0 for n in nums):
            raise IOError("PDB version 2 header illegible")
        bias = nums[:2]  # float, double exponent bias
        chart, symtab = nums[2:]
        # Unpack tentative primitive information from header.
        psz, ssz, isz, lsz, fsz, dsz = priminfo[0:6]
        sord, iord, lord = priminfo[6:9]  # 1 for big, 2 for little endian
        i = 9 + fsz
        j = i + dsz
        ford = tuple(priminfo[9:i])  # 1...sz permutation, position as in data,
        dord = tuple(priminfo[i:j])  # ...value is position in big endian order
        fbits = priminfo[j:j+7]      # [ N  e#  s#  -&  e&  s&  1? ]
        dbits = priminfo[j+7:j+14]
        # (itemsize, order, align)
        # order == 1 for big-endian, 2 for little-endian, 0 unknown
        primtypes = ((b'char', (1, 0, 0)),
                     (b'short', (ssz, sord, 0)),
                     (b'integer', (isz, iord, 0)),
                     (b'long', (lsz, lord, 0)),
                     (b'float', (fsz, ford, 0, tuple(fbits+bias[0:1]))),
                     (b'double', (dsz, dord, 0, tuple(dbits+bias[1:2]))),
                     (b'*', (psz, lord, 0)))

    else:
        m = _magic1.match(header)
        if not m:
            try:
                f.seek(-4096, 2)
            except IOError:
                f.seek(0)  # assume file shorter than 4096 bytes
            trailer = f.read(4096)
            i = trailer.rfind(b'StructureChartAddress:')
            m = _magic3.match(trailer[i:])
            if not m:
                raise IOError("file is not any version PDB file")
            # Version III PDB format, not legible by yorick.
            root.pdb_version = version = 3
            try:
                chart, symtab = map(int64, m.group(1, 2))
            except ValueError:
                chart, symtab = -1, -1
            if chart <= 0 or symtab <= 0:
                raise ValueError("PDB version 3 trailer illegible")
        else:
            # Legacy version 1 PDB format, no support for writing.
            root.pdb_version = version = 1
            try:
                nums = map(int64, m.group(1, 2, 3))
            except (AttributeError, ValueError):
                nums = [-1]
            if any(n < 0 for n in nums):
                raise IOError("PDB version 1 header illegible")
            platform, chart, symtab = nums
            primtypes = _1layouts.get(platform)
            if primtypes is None:
                primtypes = _1layouts[5]

    # Next, read chart and symbol table from end of file.
    if symtab < chart:
        # PDBLib (any version) fails unless symbol table immediately after
        # structure chart -- effectively symtab just gives length of chart.
        raise IOError("PDB file chart must precede symtab")
    f.seek(chart)
    handle.declared(handle.zero_address() | int64(chart), None, 0)
    chart_contents = f.read(symtab - chart)  # "chart" is PDB type table
    symtab_contents = f.read()
    if getattr(hooks, "pre_parse", None):
        chart_contents, symtab_contents = hooks.pre_parse(
            chart_contents, symtab_contents, version)
    if version == 3:
        _parse3(f, root, chart_contents, symtab_contents)
        return
    errors = []

    # Parse chart_contents into an OrderedDict of compounds,
    # typename --> size, members
    #   members is an OrderedDict, membername --> typename, shape
    structs = OrderedDict()
    pointertypes = set()
    chartprims, dundertype = [], None
    for line in _line_re.finditer(chart_contents):
        line = line.group(1)
        if line == b'\002':
            break
        members = line.split(b'\001')[:-1]  # must end with \001
        try:
            name = members[0].strip()
            if (not name) or (name in structs):
                errors.append("empty or repeated type name: {}"
                              "".format(name.decode('latin1')))
                raise ValueError
            size = int64(members[1])
            members = members[2:]
            if members:
                mlist, members = members, OrderedDict()
                for mem in mlist:
                    m = _memberdef.match(mem)
                    typ, ind, nm, dims = m.group(1, 2, 3, 4)
                    if nm in members:
                        errors.append("repeated member name: {}"
                                      "".format(nm.decode('latin1')))
                        raise AttributeError
                    if ind:
                        typ += b'*' * ind.count(b'*')
                        pointertypes.add(name)
                    elif typ in pointertypes:
                        pointertypes.add(name)
                    shape = _parse_shape(dims)
                    members[nm] = typ, shape
                structs[name] = size, members
                if name == b'__':
                    dundertype = members
            else:
                # Ultimately we will ignore these, but go ahead and
                # collect them here for informational purposes.
                chartprims.append((name, (size, 0, 0)))
        except (AttributeError, ValueError, IndexError):
            continue

    # Parse the symtab_contents into an OrderedDict,
    # symbolname --> address, typename, shape
    symtab = []
    haspointers = maybe_pointers = 0
    has_dirs = False
    recordsym = None
    basisnames, basismap, basisrecs = set(), {}, set()
    symtab_contents = _line_re.finditer(symtab_contents)
    for line in symtab_contents:
        line = line.group(1)
        if not line:
            break
        sym = line.split(b'\001')
        name = sym[0].strip()
        try:
            lens = list(map(int64, sym[2:-1]))  # always ends with ''
            count, addr = lens[0:2]
        except (IndexError, ValueError):
            errors.append("{} has bad symtab entry"
                          "".format(name.decode('latin1')))
            continue
        shape = lens[3::2]  # ignore index origins
        if not shape and count > 1:
            shape = [count]
        if prod(shape) == count:
            typ = sym[1].strip()
            if typ.endswith(b'*'):
                # Remove optional whitespace between pointer * characters.
                typ = sym[1].replace(b' ', b'').replace(b'\t', b'')
                haspointers |= 1
            elif (dundertype and typ == b'__' and name.endswith(b'__@history')
                  and not shape):
                recordsym = name
            elif name.count(b'@') == 1 and not atnames:
                bname, pkg = name.split(b'@')
                if pkg in (b'macro', b'funct'):
                    continue  # skip basis macros or functions
                # Strip @pkg if no bname conflict.
                if bname not in basisnames:
                    basisnames.add(bname)
                    basismap[name] = bname
                    name = bname
                if pkg == b'history':
                    basisrecs.add(name)
            if allrecs and name not in basisrecs:
                basisrecs.add(name)
            if not has_dirs:
                has_dirs |= typ == b'Directory' or name.startswith(b'/')
            # Begin workaround of yorick bug that sometimes forgot to write
            # string (and maybe pointer?) to Primitive-Types extra.
            if not haspointers:
                if typ == b'string':
                    maybe_pointers |= 1
                if typ == b'pointer':
                    maybe_pointers |= 2
            symtab.append((name, (addr, typ, tuple(shape))))
        else:
            errors.append("{} has count mismatch"
                          "".format(name.decode('latin1')))
    # Adjust symbol names to avoid conflicts between yorick record variables
    # and static variables.
    if recordsym:
        names = set(dundertype)
        names.update(s[0] for s in symtab)
        if len(names) < len(dundertype) + len(symtab):
            for i, (name, desc) in enumerate(symtab):
                if name not in dundertype:
                    continue
                # Append number to end of name, incrementing until unused.
                prefix, j = name, 0
                while True:
                    name = prefix + str(j).encode()
                    if name not in names:
                        break
                    j += 1
                symtab[i] = name, desc
                names.add(name)
    # PDBlib (but not yorick) outputs symtab in random hash order.  Sort it
    # into order of increasing address, which is probably declaration order.
    symtab.sort(key=lambda x: x[1][0])
    symtab = OrderedDict(symtab)

    # Continue parsing symtab_contents to extract extras metadata.
    eob = (b'', b'\002')
    extras = {}
    name = body = datebug = None
    nerrs = 0
    for line in symtab_contents:
        line = line.group(1)
        if (not line) and (not name):
            if datebug:
                datebug = False
                # Some files were written with a bug that included the
                # '\n' at the end of the date from the POSIX ctime call.
                # This causes a \n\n at the end of the Version: extra.
                continue
            break
        if line in eob:
            if name:
                extras[name.strip()] = body
            name = body = None
        elif not name:
            parts = line.split(b':', 1)
            if len(parts) > 1:
                name, body = parts
                if body:
                    name = name.strip()
                    extras[name] = [body]
                    datebug = name == 'Version'
                    name = body = None
                else:
                    body = []
            else:
                nerrs += 1
        else:
            body.append(line)
    if nerrs:
        errors.append("{} missing extra block name(s)".format(nerrs))

    # Handle Alignment extra.
    # The alignment extra may be overridden later by the Primitive-Types
    # extra, so interpret it first.
    palign = extras.pop(b'Alignment', [b''])[0]
    if len(palign) == 7:
        # These override values in types argument (from file header).
        palign = map(ord, palign) if PY2 else [p for p in palign]
        # reorder to match primitives in primtypes
        palign = palign[:1] + palign[2:] + palign[1:2]
        ptypes = []
        for (name, soad), a in zip(primtypes, palign):
            soad = list(soad)
            soad[2] = a
            ptypes.append((name, tuple(soad)))
        primtypes = tuple(ptypes)
    elif version == 2:
        errors.append("bad or missing Alignment extra")

    # Major-Order   (MajorOrder in PDB-3)
    # Struct-Alignment   (StructAlignment in PDB-3)
    # Offset   (DefaultIndexOffset in PDB-3)
    if extras.pop(b'Major-Order', [b'101'])[0].strip() in (b'102', b'column'):
        # Little endian shapes, reverse all shape lists.
        _flip_shapes(structs, symtab)
    extras.pop(b'Offset', None)  # ignore all minimum index values
    structal = extras.pop(b'Struct-Align', [b'0'])[0]
    try:
        structal = int(structal)
    except ValueError:
        structal = -1
    if structal < 0:
        errors.append("bad Struct-Alignment extra")
        structal = 0
    if not structal:
        structal = 1
    # Other potentially interesting extras:
    # Dynamic Spaces: <n>\n      (ia_<n> is next pointer?)
    # Use Itags: 1 or 0\n
    # Previous-File: <name>\n

    # Parse Blocks extra.
    garbled = False
    block_lines = iter(extras.pop(b'Blocks', ()))
    name = None
    # Some block addresses may not be computable until after all data type
    # item sizes are known.
    deferred = {}
    for block in block_lines:  # block_lines may increment in loop body
        # Blocks:\n
        # <name>\001<nblocks>\x<addr> <nitems>...\n    where \x = \n or space
        # \002\n
        #   -- nitems is multiple of declared dimensions except slowest
        name = block.split(b'\001', 1)
        if len(name) != 2:
            if not garbled:
                errors.append("garbled line(s) in Blocks extra")
                garbled = True
            continue
        name, n = name
        try:
            n = list(map(int64, n.split()))
        except ValueError:
            errors.append("bad block count for {}"
                          "".format(name.decode('latin1')))
            continue
        n, addcnt = n[0], n[1:]
        n += n
        if n > len(addcnt):
            for line in block_lines:  # increment block_lines iterator
                try:
                    addcnt += list(map(int64, line.split()))
                except ValueError:
                    errors.append("bad block address for {}"
                                  "".format(name.decode('latin1')))
                    break
                if len(addcnt) >= n:
                    break
            else:
                break
        if n != len(addcnt):
            errors.append("block count address mismatch for {}"
                          "".format(name.decode('latin1')))
            continue
        addcnt = array(addcnt)
        addr, count = addcnt[0::2], addcnt[1::2]
        name = basismap.get(name, name)
        basisrecs.discard(name)
        sym = symtab.get(name)
        if not sym:
            errors.append("block for non-existent variable {}"
                          "".format(name.decode('latin1')))
            continue
        shape = sym[2]
        if not shape:
            shape = (1,)
        # Without if len test, chunk and later count become floats.
        chunk = prod(shape[1:]) if len(shape) > 1 else 1
        if ((count % chunk).any() or (count[0] < shape[0]*chunk)
                or (addr[0] != sym[0])):
            errors.append("block shape disagreement for {}"
                          "".format(name.decode('latin1')))
            continue
        count //= chunk  # number of slowest index positions
        # QnD interface wants simple list of block addresses.
        if (count == 1).all():
            # Address of each block listed separately, so addr list
            # can be used as-is.
            # OrderedDict guarantees that replacing item value does not
            # change its position in the sequence.
            symtab[name] = addr.tolist(), sym[1], shape[1:]
        else:
            # Otherwise, we need to defer converting symtab[name] to
            # blocks until we know the number of bytes per chunk.
            deferred[name] = addr, chunk, count
        name = None
    if name is not None:
        errors.append("block count address mismatch for {}"
                      "".format(name.decode('latin1')))
    # Treat basis variables marked with @history as record variables
    # even if they have no blocks.
    for name in basisrecs:
        sym = symtab.get(name)
        if sym is None:
            continue  # Impossible to get here?
        addr, typ, shape = sym
        if not isinstance(addr, list):  # Impossible to fail this test?
            addr = [addr]
        shape = shape[1:] if shape and shape[0] == 1 else shape
        symtab[name] = addr, typ, shape

    # Parse Primitive-Types extra to an OrderedDict:
    # typename --> size, order, align
    #          or  size, order, align, fpformat
    ptypes = OrderedDict()
    garbled = False
    for line in extras.pop(b'Primitive-Types', []):
        line = line.split(b'\001')[:-1]   # must end with \001
        name = line[0].strip()
        try:
            if not name or len(line) < 5:
                if not garbled:
                    errors.append("garbled line(s) in Primitive-Types extra")
                garbled = True
                raise ValueError
            size, align, order = list(map(int, line[1:4]))
            flag = line[4]
            i = 5
            if flag == b'ORDER':
                i += size
                if len(line) < i:
                    errors.append("bad ORDER in primitive {}"
                                  "".format(name.decode('latin1')))
                    raise ValueError
                order = tuple(map(int, line[5:i]))
            elif flag == b'DEFORDER':
                order = tuple(range(1, size+1) if order == 1 else
                              range(size, 0, -1))
            else:
                errors.append("unknown ORDER in primitive {}"
                              "".format(name.decode('latin1')))
                raise ValueError
            flag, line = line[i], line[i+1:]
            if flag == b'FLOAT':
                if len(line) < 8:
                    errors.append("bad FLOAT in primitive {}"
                                  "".format(name.decode('latin1')))
                    raise ValueError
                fbits = tuple(map(int, line[:8]))
                ptypes[name] = size, order, align, fbits
                line = line[8:]
            elif (flag in (b'FIX', b'NO-CONV')):
                if flag == b'NO-CONV':
                    order = ()
                ptypes[name] = size, order, align
            else:
                errors.append("unrecognized primitive type {}"
                              "".format(name.decode('latin1')))
                raise ValueError
            # line may still contain UNSGNED an ONESCMP items
            # also things like (TUPLE, float_complex, 2, -1) ?!
        except (ValueError, IndexError):
            errors.append("garbled primitive type {}"
                          "".format(name.decode('latin1')))

    # Yorick writes idiosyncratic chart and PrimitiveTypes, which is
    # important for how it interprets its own brand of pointers.
    # Yorick writes first writes seven primitives in order:
    #   * short integer long float double char
    # Yorick writes no other primitives in the chart; rather they are
    # written to the Primitive-Types extra.
    # This differs from the order of the standard type indices it uses
    # when writing pointer values:
    #   char short int long float double string pointer
    # Two additional yorick internal virtual primitives come after these:
    #   char * (PDB member pointer) and  char* (PDB symbol pointer)
    # These 10 types have predefined type indices 0-9 in yorick pointer
    # data.  The remaining type indices come first from the Primitive-Types
    # extra, if and only if referenced (Directory first if used), then
    # from the order of non-primitive compound types in the chart.
    # (Note that yorick's string and pointer primitives appear only if
    # used in the file.  They are aliases for the long type primitive.)
    # Yorick does not write the standard primitives in the chart, PDBLib does.
    # Yorick writes complex data as struct complex {double re, im;}.

    # We have primtypes (char, short, int long, float, double, *),
    # chartprims, and ptypes (from Primitive-Types extra).  The chartprims
    # are nearly useless, since they have no order or alignment information.

    # Yorick assumes that all chartprims are mentioned either as one of
    # the standard primtypes or in the ptypes from the Primitive-Types extra,
    # so we simply ignore them here (in yorick they can only generate errors).
    # Convert primtypes to an OrderedDict and append ptypes to it, keeping
    # the original primtypes at the beginning, and adding any additional
    # primitives from the extras in order afterwards.  The Primitive-Types
    # extra supersedes the original primtypes.
    primtypes = OrderedDict(primtypes)
    primtypes.update(ptypes)
    if not haspointers:
        string, pointer = ptypes.get(b'string'), ptypes.get(b'pointer')
        if not string and maybe_pointers & 1 and b'string' not in structs:
            primtypes[b'string'] = string = primtypes[b'long']
        if not pointer and maybe_pointers & 2 and b'pointer' not in structs:
            primtypes[b'pointer'] = pointer = primtypes[b'long']
        desc = string or pointer
        if string and pointer and string != pointer:
            desc = None
        if desc:
            # At least one of string or pointer primitive present and
            # both the same if both present.  Check that this common type
            # is identical to long primitive.
            ldesc = primtypes[b'long']
            if ldesc[0::2] == desc[0::2]:
                # Orders hard to compare, may be either permuation or flag.
                lord, pord = ldesc[1], desc[1]
                if lord == pord:
                    haspointers = 2
                elif isinstance(lord, tuple) and not isinstance(pord, tuple):
                    sord = tuple(range(1, len(lord)+1))
                    if pord == 2:
                        sord = sord[::-1]
                    elif pord != 1:
                        sord = ()
                    if sord == lord:
                        haspointers = 2
                elif isinstance(pord, tuple) and not isinstance(lord, tuple):
                    sord = tuple(range(1, len(pord)+1))
                    if lord == 2:
                        sord = sord[::-1]
                    elif lord != 1:
                        sord = ()
                    if sord == pord:
                        haspointers = 2

    _endparse(root, structal, haspointers, primtypes, structs, symtab, errors,
              deferred, recordsym)


def _endparse(root, structal, haspointers, primtypes, structs, symtab, errors,
              deferred, recordsym=None):
    chart = root.chart
    chart.structal = structal
    chart.haspointers = haspointers
    # TODO: In PDB-3 it is possible that primtypes does not contain a
    # definition for int, or indeed for any primitive not used in symtab.
    chart.byteorder = primtypes[b'long'][1]  # of long

    # Make sure that Directory primitive goes first if present.
    # This ensures yorick pointer type index will correspond to
    # position in chart.primitives OrderedDict.
    dirdesc = primtypes.get(b'Directory')
    if dirdesc is not None:
        chart.add_primitive(b'Directory', dirdesc)
    type_mismatch = False
    for name, desc in itemsof(primtypes):
        if chart.add_primitive(name, desc) is None:
            type_mismatch = True
    if recordsym:
        dundertype = structs.pop(b'__')  # make sure this is last
    for name, desc in itemsof(structs):
        if chart.add_struct(name, desc) is None:
            type_mismatch = True
    if recordsym:
        if chart.add_struct(b'__', dundertype) is None:
            type_mismatch = True
        else:
            dundertype = chart.structs.pop(b'__')
    if type_mismatch:
        raise IOError("data type mismatch building PDB structure chart")

    if recordsym:
        # Convert yorick default __@history style records to less
        # efficient native PDB blocks style record variables.
        raddr = symtab.pop(recordsym)[0]
        stype, members = dundertype[0], dundertype[3]
        # Names in symtab were adjusted if necessary to not conflict with
        # field names in stype.
        fields = stype.fields
        offsets = [fields[nm][1] for nm in stype.names]
        raddr = array(raddr)
        for off, (sname, (tname, shape)) in zip(offsets, itemsof(members)):
            symtab[sname] = (raddr + off).tolist(), tname, shape

    primitives, structs = chart.primitives, chart.structs
    undefined = set()
    groups = {b'': root}
    addr0 = root.handle.zero_address()
    attributes = OrderedDict()
    attr_address = None;
    for name, (addr, tname, shape) in itemsof(symtab):
        defblock = deferred.get(name)
        if name.startswith(b'/'):
            dirname, name = name.rsplit(b'/', 1)
            grp = groups.get(dirname)
            if grp is None:
                grp = _declare_group(groups, dirname)
        else:
            grp = root
        if tname == b'Directory':
            continue
        if not name:
            errors.append("{}/ not type Directory"
                          "".format(dirname.decode('latin1')))
            name = b'?'  # soldier on with bogus name
        dtype = primitives.get(tname)
        if dtype is None:
            dtype = structs.get(tname)
            if dtype is None:
                if tname not in (b'string', b'pointer'):
                    undefined.add(tname)
                    continue
                # Work around yorick bug that sometimes forgets to write
                # Primitive-Types extras.
        if defblock:
            # Complicated Blocks: extra could not be processed until now
            # that we know the number of bytes per item.
            addr = _make_simple_list(dtype[0].itemsize, *defblock)
        unlim = hasattr(addr, '__iter__')  # match test in PDBLeaf.__init__
        if unlim:
            addr = (array(addr) + addr0).tolist()
        else:
            addr = addr + addr0
        if not PY2:
            name = name.decode('latin1')
        stype, dtype, align, _ = dtype  # from primitives[] or structs[]
        dtype = dtype, stype, align, tname
        if name.startswith(":") :
            if attr_address is None or addr < attr_address:
                attr_address = addr
            attributes[name[1:]] = dtype, shape, addr
        else:
            grp._leaf_declare(name, dtype, shape, addr)
    if undefined:
        errors.append("undefined types: {}".format(undefined))
    if len(attributes):
        root.has_attributes = True
        varname, vname, var = "", "", root
        for aname, (dtype, shape, addr) in itemsof(attributes):
            prevname = varname
            varname = aname.rsplit(":", 1)
            bad = len(varname) < 2
            if not bad:
                varname, name = varname  # name is attribute name
                if varname != prevname:  # otherwise grp and vname unchanged
                    path = varname.split(":")
                    vname = path.pop()
                    grp = root
                    for gname in path:
                        grp = grp.lookup(gname)
                        if grp.islist() == 2:
                            grp = var.parent()
                        elif not grp.isgroup():
                            bad = True
                            break
            if bad:
                errors.append(
                    "cannot find attribute: {}".format(aname))
                continue
            if name != "0":
                grp._read_attr(vname, name, dtype[0], shape, addr)
            elif grp._read_shape(vname, dtype[0], shape, addr):
                errors.append(
                    "bad zero-length array indicator: {}".format(aname))
        # Ensure that attributes are overwritten if file is extended:
        handle = root.handle
        handle.declared(handle.zero_address() | int64(attr_address), None, 0)
    if errors:
        # Can filter this by module name.
        print(errors)
        warn("{} errors parsing PDB metadata".format(len(errors)))


def _declare_group(groups, dirname):
    dname, name = dirname.rsplit(b'/', 1)
    grp = groups.get(dname)
    if grp is None:  # recurse to define all undefined ancestors
        grp = _declare_group(groups, dname)
    if not PY2:
        name = name.decode('latin1')
    grp = grp.declare(name, dict, None)
    groups[dirname] = grp
    return grp


_magic1 = re.compile(br'!<><PDB><>!\s*([0-9]+)\s*[\r\n\037]'
                     br'\s*([0-9]+)\s*[\r\n\037]\s*([0-9]+)\s*[\r\n\037]')

_magic2 = re.compile(br'!<<PDB:II>>![\r\n\037](.)')
_magic2a = re.compile(br'\s*([0-9]+)\001\s*([0-9]+)\001[\r\n\037]'
                      br'\s*([0-9]+)\001\s*([0-9]+)\001[\r\n\037]')

_magic3 = re.compile(br'StructureChartAddress:\s*([0-9]+)\s*[\r\n\037]+'
                     br'SymbolTableAddress:\s*([0-9]+)\s*[\r\n\037]+'
                     br'!<<PDB:3>>![\r\n\037]+')

# IEEE 754-2008 float, double, and quad precision formats
_binary32 =  ( 32,  8,  23, 0, 1,  9, 0,   127)  # noqa IEEE 754 float
_binary64 =  ( 64, 11,  52, 0, 1, 12, 0,  1023)  # noqa IEEE 754 double
_binary128 = (128, 15, 112, 0, 1, 16, 0, 16383)  # noqa IEEE 754 quadruple
# Intel 80 bit format is an example of IEEE 754-2008 binary64ext
_intel80 =   ( 80, 15,  64, 0, 1, 16, 1, 16383)  # noqa Intel 80 bit registers
_oldmac96 =  ( 96, 15,  64, 0, 1, 32, 1, 16382)  # noqa ??
# Remaining formats have different rules for Inf, Nan, and denormal
_cray64 =    ( 64, 15,  48, 0, 1, 16, 1, 16384)  # noqa Cray 1, XMP, YMP
_vaxf32 =    ( 32,  8,  23, 0, 1,  9, 0,   129)  # noqa VAX F float
_vaxd64 =    ( 64,  8,  55, 0, 1,  9, 0,   129)  # noqa VAX D double
_vaxg64 =    ( 64, 11,  52, 0, 1, 12, 0,  1025)  # noqa VAX G double
_vaxh128 =   (128, 15, 112, 0, 1, 16, 0, 16385)  # noqa VAX H quad

# legacy primitive types for PDB version 1
_1layouts = {1: ((b'char',    (1, 0, 0)), (b'short', (2, 1, 2)),   # sun3
                 (b'integer', (4, 1, 2)), (b'long',  (4, 1, 2)),
                 (b'float',  (4, 1, 2, _binary32)),
                 (b'double', (8, 1, 2, _binary64)), (b'*', (4, 1, 4))),
             2: ((b'char',    (1, 0, 0)), (b'short', (2, 2, 2)),  # ibmpc
                 (b'integer', (2, 2, 2)), (b'long',  (4, 2, 2)),
                 (b'float',  (4, 2, 2, _binary32)),
                 (b'double', (8, 2, 2, _binary64)), (b'*', (4, 2, 4))),
             3: ((b'char',    (1, 0, 0)), (b'short', (8, 1, 2)),  # cray
                 (b'integer', (8, 1, 8)), (b'long',  (8, 1, 8)),  # hybrid cc
                 (b'float',  (8, 1, 8, _cray64)),
                 (b'double', (8, 1, 8, _cray64)), (b'*', (8, 1, 8))),
             4: ((b'char',    (1, 0, 0)), (b'short', (2, 2, 1)),  # vax
                 (b'integer', (4, 2, 1)), (b'long',  (4, 2, 1)),
                 (b'float',  (4, (2, 1, 4, 3), 1, _vaxf32)),
                 (b'double', (8, (2, 1, 4, 3, 6, 5, 8, 7), 1, _vaxg64)),
                 (b'*', (4, 2, 1))),
             5: ((b'char',    (1, 0, 0)), (b'short', (2, 1, 4)),  # PDB default
                 (b'integer', (4, 1, 4)), (b'long',  (4, 1, 4)),
                 (b'float',  (4, 1, 4, _binary32)),
                 (b'double', (8, 1, 4, _binary64)), (b'*', (4, 1, 4))),
             6: ((b'char',    (1, 0, 0)), (b'short', (2, 1, 2)),  # mac2
                 (b'integer', (2, 1, 2)), (b'long',  (4, 1, 2)),
                 (b'float',  (4, 1, 2, _binary32)),
                 (b'double', (12, 1, 2, _oldmac96)), (b'*', (4, 1, 2)))}

# Though probably never actually used, PDBlib source code will accept
# ASCII 037 = US = unit separator as a newline character.
# Curiously, PDB treats DOS \r\n as a blank line following the line,
# which would break the PDB metadata.
_line_re = re.compile(br'([^\r\n\037]*)[\r\n\037]')

dim = br'\s*(\d+(?::\d+)?(?:\s*,\s*\d+(?::\d+)?)*)\s*'  # n[:n][,n[:n]]*
dim = br'\s*(?:[\[(]' + dim + br'[])]\s*)?'
_memberdef = re.compile(br'\s*([^[(* \t]+)\s((?:\s*\*)*)'
                        br'([^[(* \t]+)' + dim + br'\s*$')
_ptrheader = re.compile(br'(\d+)\001([^\001]+)\001'
                        br'(?:([-+]?\d+)\001)?(?:(\d+)\001)?[\r\n\037]')

_pdb3_tagline = re.compile(br'^\s*(\S+)\s*:\s*(\S+)?\s*$')
_pdb3_prim = re.compile(br'^\s*(\S+)\s+(\d+)\s+(\d+)\s+([^;]+);\s*$')
_pdb3_pattr = re.compile(br'\s*(FIX|NO-CONV|UNSGNED|ONESCMP|'
                         br'TYPEDEF\s*\(\s*(\S+)\s*\)|'
                         br'ORDER\s*\(\s*(big|little|\d+(?:\s*,\s*\d+)*)\s*\)|'
                         br'FLOAT\s*\(\d+(?:\s*,\s*\d+){7}\s*\))')
_pdb3_struct = re.compile(br'^\s*(\S+)\s*\(\s*(\d+)\s*\)\s*$')
_pdb3_memb = re.compile(br'\s*(\{)?\s*([^ \t*]+)\s*(\*(?:\s*\*)*)?\s*'
                        br'([^[( \t<;]+)' + dim + br'(?:<-\s*\S+\s*)?'
                        br';\s*(\}\s*;\s*)?$')
_pdb3_symb = re.compile(br'\s*([^ \t*]+)\s*(\*(?:\s*\*)*)?\s*([^ \t<;]+)' +
                        dim + br'@\s*(\d+)\s*\(\s*(\d+)\s*\)\s*;\s*$')
del dim
_pdb3_block = re.compile(br'^\s*(\S+)\s*(\d+)\s*$')


def _parse_shape(shape):
    if not shape:
        return ()
    dims = []
    for d in (s.split(b':') for s in shape.split(b',')):
        v, d = int(d[0]), d[1:]
        dims.append(max(int(d[0]) - v + 1, 0) if d else v)
    return tuple(dims)


def _flip_shapes(chart, symtab):
    for name in chart:
        members = chart[name][1]
        if members:
            for mname in members:
                typ, shape = members[mname]
                members[mname] = typ, shape[::-1]
    for name in symtab:
        addr, typ, shape = symtab[name]
        symtab[name] = addr, typ, shape[::-1]


def _make_simple_list(itemsize, addrs, chunk, count):
    # Convert counted list of chunk addresses into a simple address list
    # with one chunk per element.
    # The crucial feature of this algorithm is that there is no explicit
    # interpreted loop; it is entirely numpy array operations, so the
    # number of chunks can be quite large without a significant performance
    # impact.  Equivalent to:
    #   output = []
    #   for a, c in zip(addrs, count):
    #       output.extend(a + arange(c)*chunk)
    chunk *= itemsize
    n, na = count.sum(), addrs.size
    if na < n:
        indx = concatenate(([0], count))[:-1].cumsum()
        ibcast = zeros(n, int64)
        ibcast[indx] = 1
        ibcast = ibcast.cumsum() - 1
        icount = arange(n) - indx[ibcast]
        addrs = addrs[ibcast] + chunk*icount
    return addrs.tolist()


def _parse3(f, root, chart_contents, symtab_contents):
    chart_contents = _line_re.finditer(chart_contents)
    symtab_contents = _line_re.finditer(symtab_contents)
    errors = []

    # Parse the version 3 structure chart.
    ptypes = OrderedDict()
    structs = OrderedDict()
    pointertypes = set()
    structal = skippedprim = skippedcomp = 0
    fortran = False
    openstruct = tag = members = None
    for line in chart_contents:
        match = _pdb3_tagline.match(line)
        if match:
            if openstruct:
                errors.append("struct {} definition incomplete"
                              "".format(openstruct.decode('latin1')))
                openstruct = None
            tag = match.group(1)
            # Each tag is a chart section.  They should appear in this order,
            # although we do not enforce the order here.
            if tag == b'PrimitiveTypes':
                if ptypes:
                    errors.append("multiple PrimitiveTypes sections in chart")
                tag = b':'
            elif tag == b'StructAlignment':
                structal = match.group(2)
                if structal.isdigit():
                    structal = int(structal)
                else:
                    errors.append("bad StructAlignment value")
                    structal = 0
                tag = None
            elif tag == b'DefaultIndexOffset':
                # Ignore all index origin suggestions.
                tag = None
            elif tag == b'MajorOrder':
                fortran = match.group(2) == b'column'
                tag = None
            elif tag == b'CompoundTypes':
                if structs:
                    errors.append("multiple CompoundTypes sections in chart")
                tag = b''
            else:
                tag = None
            continue
        if tag is None:
            continue

        if tag == b':':
            # Parse PrimitiveTypes section.
            # typename --> size, order, align
            #          or  size, order, align, fpformat
            match = _pdb3_prim.match(line)
            if not match:
                if line.strip():
                    skippedprim += 1
                continue
            name, size, align, attribs = match.group(1, 2, 3, 4)
            size, align = int(size), int(align)
            fpbits = None
            order = 0
            isfix, typedefs = False, []
            for attr in attribs.split(b'|'):
                match = _pdb3_pattr.match(attr)
                if not match:
                    continue  # just ignore malformed attributes
                a = match.group(1)
                if a.startswith(b'FL'):
                    fpbits = tuple(int(p) for p in match.group(4).split(b','))
                elif a.startswith(b'O'):
                    order = match.group(3)
                    if order.startswith(b'b'):
                        order = 2
                    elif order.startswith(b'l'):
                        order = 1
                    else:
                        order = tuple(int(p) for p in order.split(b','))
                        if 0 in order:  # PDB-II permutations are 1-origin
                            order = tuple(p+1 for p in order)
                elif a.startswith(b'T'):
                    typedefs.append(match.group(2))
                elif a.startswith(b'FI'):
                    isfix = True
            if isfix and fpbits:
                fpbits = None
                errors.append("{} declared both fix and float in chart"
                              "".format(name.decode('latin1')))
            if (isfix or fpbits) and size > 1 and not order:
                errors.append("{} does not specify byte order in chart"
                              "".format(name.decode('latin1')))
            if fpbits:
                ptypes[name] = size, order, align, fpbits
            else:
                ptypes[name] = size, order, align
            continue

        # Parse CompoundTypes section.
        # typename --> size, members
        #   members is an OrderedDict, membername --> typename, shape
        if not openstruct:
            match = _pdb3_struct.match(line)
            if not match:
                if line.strip():
                    skippedcomp += 1
                continue
            openstruct, size = match.group(1, 2)
            size = int64(size)
            members = OrderedDict()
            continue
        match = _pdb3_memb.match(line)
        if match.group(1):  # has opening brace
            if members:
                errors.append("{} has multiple opening braces"
                              "".format(openstruct.decode('latin1')))
        elif not members:
            errors.append("{} missing open brace"
                          "".format(openstruct.decode('latin1')))
        mtype, ind, mname, shape = match.groups(2, 3, 4, 5)
        if ind:
            mtype += b'*' * ind.count(b'*')
            pointertypes.add(openstruct)
        elif mtype in pointertypes:
            pointertypes.add(openstruct)
        shape = _parse_shape(shape)
        members[mname] = mtype, shape
        if match.groups(6):  # has closing brace
            structs[openstruct] = size, members
            openstruct = members = None
    if openstruct:
        errors.append("struct {} definition incomplete"
                      "".format(openstruct.decode('latin1')))
    if skippedprim:
        errors.append("skipped {} lines parsing PrimitiveTypes"
                      "".format(skippedprim))
    if skippedcomp:
        errors.append("skipped {} lines parsing CompoundTypes"
                      "".format(skippedcomp))

    # Parse the symtab_contents into an OrderedDict,
    # symbolname --> address, typename, shape
    symtab = OrderedDict()
    has_dirs = False
    haspointers = garbage = 0
    slowest = -1 if fortran else 0
    tag = None
    openblock = blocks = nbspecs = None
    deferred = {}
    for line in symtab_contents:
        line = line.group(1)
        match = _pdb3_tagline.match(line)
        if match:
            if openblock:
                errors.append("blocks for {} incomplete"
                              "".format(openblock.decode('latin1')))
                openblock = blocks = None
            tag = match.group(1)
            if tag == b'StructureChartAddress':
                break
            if tag == b'SymbolTable':
                if symtab:
                    errors.append("multiple SymbolTable sections in symtab")
                tag = b':'
            elif tag == b'Blocks':
                if symtab:
                    tag = b'b:'
                else:
                    errors.append("Blocks before SymbolTable in symtab")
                    tag = None
            elif tag == b'Checksums':
                tag = None  # Ignore checksums for now.
            # DynamicSpaces: <n>\n      (ia_<n> is next pointer?)
            # UseItags: 1 or 0\n
            # PreviousFile: <name>\n
            else:
                tag = None
            continue
        if tag is None:
            continue

        if tag == b':':
            # Parse SymbolTable entry
            match = _pdb3_symb.match(line)
            if not match:
                if line.strip():
                    garbage += 1
                continue
            typ, ind, sname, shape, addr, size = match.group(1, 2, 3, 4, 5, 6)
            addr, size = int64(addr), int64(size)
            if not shape and size > 1:
                shape = (size,)
            else:
                shape = _parse_shape(shape)
            if prod(shape) == size:
                if ind:
                    typ += b'*' * ind.count(b'*')
                    haspointers |= 1
                elif typ in pointertypes:
                    haspointers |= 1
                symtab[name] = addr, typ, shape
                if not has_dirs:
                    has_dirs |= typ == 'Directory' or name.startswith(b'/')
            else:
                errors.append("{} has count mismatch"
                              "".format(name.decode('latin1')))

        elif tag == b'b:':
            # Parse next line of blocks table.
            match = _pdb3_block.match(line)
            if not match:
                if line.strip():
                    garbage += 1
                continue
            addr, count = match.group(1, 2)
            count = int64(count)
            if not openblock:
                if count:
                    openblock, nbspecs = addr, count
                    blocks = []
                continue
            if not addr.isdigit():
                errors.append("garbled blocks for {}"
                              "".format(openblock.decode('latin1')))
                # Attempt to resynchronize, but probably hopelessly lost.
                openblock, nbspecs = addr, count
                continue
            addr = int64(addr)
            blocks.append((addr, count))
            nbspecs -= 1
            if not nbspecs:
                name, openblock = openblock, None
                addr, count = array(blocks).T
                sym = symtab.get(name)
                shape = sym[2] if sym else ()
                if not shape:
                    errors.append("block shape mismatch for {}"
                                  "".format(name.decode('latin1')))
                    continue
                # Without if len test, chunk and later count become floats.
                chunk = prod(shape[1:]) if len(shape) > 1 else 1
                if ((count % chunk).any() or (count[0] < shape[slowest]*chunk)
                        or (addr[0] != sym[0])):
                    errors.append("block shape disagreement for {}"
                                  "".format(name.decode('latin1')))
                    continue
                count //= chunk  # number of slowest index positions
                # QnD interface wants simple list of block addresses.
                if (count == 1).all():
                    # All chunk addresses given, convert symtab entry
                    addr, typ, shape = symtab[name]
                    # OrderedDict guarantees that replacing item value does
                    # not change its position in the sequence.
                    symtab[name] = addr.tolist(), typ, shape[1:]
                else:
                    deferred[name] = addr, chunk, count
    if openblock:
        errors.append("blocks for {} incomplete"
                      "".format(openblock.decode('latin1')))
    if garbage:
        errors.append("skipped {} lines parsing symtab".format(garbage))

    # PDBlib outputs symtab in random hash order.  Sort it into order of
    # increasing address, which is probably declaration order.
    symtab = list(itemsof(symtab))
    symtab.sort(key=_addrkey)
    symtab = OrderedDict(symtab)

    if fortran:
        # Shapes are in little-endian order, reverse them.
        _flip_shapes(structs, symtab)

    _endparse(root, structal, haspointers, ptypes, structs, symtab, errors,
              deferred)


def _addrkey(item):
    addr = item[1][0]
    return addr[0] if isinstance(addr, list) else addr


# ----------------------------------------------------------------------------
# standard numpy types, present on all modern machines:

#   Array protocol (array interface) type strings:
#   (> or < depending on platform, except |i1, |u1)
# i1   i2   i4   i8
# u1   u2   u4   u8
# f4   f8
# c8   c16
#   all i and u are two's complement
#     f4 = IEEE 754-2008 binary32 format
#     f8 = IEEE 754-2008 binary64 format

#   Additionally, an unspecified binary64ext extended precision format
#   will be present, either f12 or f16, usually the 80-bit intel format
#   with 2 or 6 unused bytes.  On some non-intel platforms, the f16
#   format will be the specific IEEE 754-2008 binary128 format.
#   The only portable way to discover this type is as
# numpy.longdouble
#   On 32-bit intel machines, this will be f12, on 64 bit machines f16.


# ------------------------------------------------------- PDB metadata summary
# In PDB metadata, in all versions, newline \n may be \n, \r, or \037 (ESC).

# PDB-I header:
# Begin with 11 bytes:
#   header[:11] = '!<><PDB><>!'
# Following this are three lines containing one decimal number each:
#   platform \n chart_address \n symtab_address \n
# The platform numbers are detailed in the _1layouts global variable above;
# all are obsolete machines unlikely to match any modern machines.

# PDB-II header:
# Header format, newline \n may be \n, \r, or \037 (ESC):
# Begin with 13 bytes
#   header[:13] = '!<<PDB:II>>!\n'
# Followed by single byte values as follows:
#   N = ord(header[13]) = count of single byte values
#   sP, sS, sI, sL, sF, sD = ord(header[14:20]) = type byte sizes
#   oS, oI, oL = ord(header[20:23]) = integer type byte orders
#   pF = ord(header[23:23+sF]) = float byte permutation
#   pD = ord(header[23+sF:23+sF+sD]) = double byte permutation
#   fF = ord(header[23+sF+sD:30+sF+sD])
#   fD = ord(header[30+sF+sD:37+sF+sD])
# Notes:
#   N = 24 + sF + sD
#   oS, oI, oL are 1 for big-endian, 2 for little-endian order
#   pF, pD are 1-origin big-endian byte positions
#     e.g.- pF = [4, 3, 2, 1] for little-endian, [1, 2, 3, 4] for big
#     The type alignments are not specified in this header;
#     the extras later override this header primitive specification.
#   fF, fD are floating point format [ N  e#  s#  -&  e&  s&  1? ]
# Bytes header[37+sF+sD:] are four ASCII decimal numbers organized
# in two lines (with same options for newline as header[12]):
#   biasF\001 biasD\001\n chart_address\001 symtab_address\001\n
# The biases complete the floating point formats.
# The addresses are byte addresses; the symtab always immediately
# follows the chart, and the symtab and extras table extend to the
# end of the file.
#
# PDB-II chart:
# <type>\001<nbytes>\001<type> <name> [<dimlist>]\001 ...\001\n
# ... note: primitives first, including *, with no members
# \002\n
#
# PDB-II symtab:
# <name>\001<type>\001<nitems>\001<address>\001[<origin>\001<length>\001]*\n
# \n    so \n\n marks end of symtab, beginning of extras
#
# PDB-II extras related to type interpretation
#   (primitives also listed at front of version II chart):
# Primitive-Types:\n
# <type>\001<nbytes>\001<alignment>\001<order>\001<moreorder>\001<more>\001\n
#   <order> = 1 MSB first, 2 LSB first, -1 NO-CONV or FLOAT
#   <moreorder> = DEFORDER  or  ORDER\0010\0011\0012\0013
#   <more> = NO-CONV or FIX or FLOAT\001nbits\001...bias
# \002\n
# Alignment:cpsilfd\n    alignments as 7 single bytes
# StructAlignment: <n>\n
# Offset:<n>\n    default index origin (ignore)
# MajorOrder: 101 or 102\n   (row or column = first slowest, first fastest)
# Longlong-Format-Alignment:soa\n    nbytes, order, alignment as bytes
#    order = 1 big-endian, 2 little-endian, 3 text(?), 4 external, 0 none
#  When directories initialized, creates / and /&ptrs/
#  and NO-CONV 1 byte 0 align no order Directory primitive.
#
# more PDB-II extras
# Has-Directories: 1 or 0 \n
# Version:<11 for yorick>|<date>\n     bug in some files- <date> contains \n
# Blocks:\n
# <name>\001<nblocks>\x<addr> <nitems>...\n    where \x either \n or space
# \002\n
#   -- note that nitems is multiple of declared dimensions except slowest
# Checksums:
#   ...
# Dynamic Spaces: <n>\n      count of itags (ia_<n> is next pointer?)
# Use Itags: 1 or 0\n
# Previous-File: <name>\n
#
# in PDB-II, \002 is optional at end of block
# extra-name: blah-blah \n
# extra-name:\n
# blah-blah\n
# \n
# extra-name:\n
# blah-blah\n
# \002 blah-blah \n
#
# PDB-II pointer data:
# Symbols which are pointers take no space; instead the pointee data begins
# at the symbol address.  For symbols which are arrays of structs containing
# pointers, the pointer data begins at the address immediately following the
# declared struct array, proceeding in depth-first order (that is, marching
# through the top level symbol, but interrupting with the pointees of any
# pointers encountered in the top-level pointees, and so on).  Each pointee
# is preceeded by a bytes header:
# "%d\001%s\001%d\001%d\001\n", nitems, full_type, address, dataHere
#     -- a NULL pointer is represented by nitems==0, address==-1, dataHere==0
#        -- if address is -1, nitems is ignored (NULL pointer)
#        -- address is the address of the pointee header with dataHere!=0
#           not the address of the data itself (this nitems and full_type
#           are ignored, strangely)
#     -- dataHere may be missing, treated as if dataHere==1
#     -- address and dataHere may be missing, which causes header to be
#        treated as if NULL, again ignoring nitems
# If dataHere non-zero, the full_type[nitems] data begins at the byte after \n
# followed by any pointee data it may refer to.
#
# PDB-II yorick pointers:
# Pointer itself (string or general pointer) is a long containing address
# of the header.  An address <0 (-1 preferred) is a NULL pointer.
# string: Header is long byte count n, followed immediately by n data bytes.
# pointer: Header is 2+ndims longs; first two are typenumber, ndims (<=10).
#          Data follows at next data-aligned address after header.
# Typenumber is tricky; always begins with standard types char (0), short, int
# long, float, double, followed by string, pointer, followed by char *, char*,
# followed by PrimitiveTypes actually used (beginning with Directory if used),
# followed by compound types in chart order.

# PDB-3 has no header, marked near EOF with
# StructureChartAddress: <address>\n
# SymbolTableAddress: <address>\n
# !<<PDB:3>>!\n
#
# PDB-3 chart
# PrimitiveTypes:                primitive types section first, not optional
# <type> <nbytes> <alignment> <attr1>[|<attr2>]*;\n
#   ORDER(big little or 0,1,2,3)  or  NO-CONV
#   FIX  or  FLOAT(nbits,nexp,nmant,osign,oexp,omant,leadbit,bias)
#   UNSGNED   ONESCMP  TYPEDEF(<type>)
#
# Directory 1 0 NO-CONV|FIX;        in particular
#
# StructAlignment: <n>\n
# DefaultIndexOffset: <n>\n
# MajorOrder: row or column\n   (C or Fortran, respectively)
#
# CompoundTypes:                 compound types section last, not optional
# <type> (<nbytes>)\n
# {<type> <name> [<dimlist>] [<- <cast>];\n
#  <type> <name> [<dimlist>] [<- <cast>];\n
#  <type> <name> [<dimlist>] [<- <cast>];};\n
#
# PDB-3 symtab
# SymbolTable:
# <type> <name>[<dimensions>] @ <address> (<nitems>);\n
# Directory / @ 0 (1);
# Directory /&ptrs/ @ 1 (1);
# double /&ptrs/ia_1 @ 10 (10);
#
# PDB-3 extras
# Version: <n> [(date)]\n       date is file creation date, not version date
# Version: 19 (Fri Apr 28 14:04:36 2006)
#    version 24 appeared dated 02/03/2010
# Blocks:
# <name> <nblocks>\n
#   <addr> <nitems>\n
#   ...
# Checksums:
#   ...
# DynamicSpaces: <n>\n      count of itags (ia_<n> is next pointer?)
# UseItags: 1 or 0\n
# PreviousFile: <name>\n
