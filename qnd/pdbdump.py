"""Low level functions to create and write a PDB file."""
from __future__ import absolute_import

from datetime import datetime
import sys

from numpy import array, prod

from .pdbparse import _binary32, _binary64

PY2 = sys.version_info < (3,)
if PY2:

    def itemsof(d): return d.iteritems()  # noqa
else:

    def itemsof(d): return d.items()  # noqa


def flusher_for(root):

    def _flusher(f):
        return flusher(f, root)

    return _flusher


def initializer_for(root):

    def _initializer(f):
        return initializer(f, root)

    return _initializer


# header[:13] = b'!<<PDB:II>>!\n'
# header[13]  = N    count of single byte values (= 24 + sz_float + sz_double)
# header[14:20]) = sz_ptr, sz_short, sz_int, sz_long, sz_float, sz_double
# header[20:23]) = ord_shot, ord_int, ord_long  (1 = >, 2 = < order)
# header[23:23+sz_float] = (4, 3, 2, 1)  for <f4 byte permutation
# header[... +sz_double] = (8, 7, 6, 5, 4, 3, 2, 1)  for <f8
# header[14+N:21+N] = (nbits, e#, s#, -&, e&, s&, 1?) float format
# header[21+N:28+N] = (nbits, e#, s#, -&, e&, s&, 1?) double format
# header[28+N:...] = biasF\001biasD\001\n  completes float and double formats
#   ( 32,  8,  23, 0, 1,  9, 0,  biasF=  127)  for f4
#   ( 64, 11,  52, 0, 1, 12, 0,  biasD= 1023)  for f8
#   b'127\0011023\001\n'  (10 bytes)
# Since numpy assumes these f4 and f8 formats, N = 21, and header has
# a total of 59 bytes to this point.
# header[...] = chart_addr\001symtab_addr\001\n
# Yorick allows an additional 128 unused bytes here to allow for longer
# chart and symtab addresses, but that is clearly overkill.  A full 64-bit
# integer has a range of 2**64 = 1.6e19, and can thus be represented in 20
# decimal digits or fewer.  Hence, the address fields can occupy at most 43
# bytes, for a 102 byte header.  If both float and double were 16 bytes for
# some unfathomable reason, that would add 20 bytes to the header, bringing
# it to 122.  Hence, 128 bytes is essentially guaranteed to be enough for
# all possible PDB headers.


def initializer(f, root):
    handle, chart = root.handle, root.chart  # root is a PDBGroup
    ifile = handle.current_file()
    if ifile:
        # Copy header and non-record vars from file 0 of the family.
        f = handle.open(0)
        header = f.read(128)
        # We will also need to copy all non-record variables, so read
        # them all now while we have file 0 open.
        nrvars = [(item, item.read()) for item in _iter_nonrec(root)]
        f = handle.open(ifile)
        # Write wrong chart and symtab addresses, fixed by flusher later.
        f.write(header)
        handle.declared(0, None, len(header))
        # Then write out all the non-record variables we read from file 0.
        addr0 = handle.zero_address()
        for item, value in nrvars:
            item = item.shifted_copy(addr0)
            stype, shape, addr = item.tsa
            handle.declared(addr, stype[1], prod(shape) if shape else 1)
            item.write(value)
        return
    # Initializing first file of family, so we need to initialize chart.
    order, structal = chart.byteorder, chart.structal
    if order is None:
        # PDB byte order is 1 for >, 2 for <
        # Here is a literal way to compute native byte order:
        order = chart.byteorder = int(array([1]).view('u1')[0] + 1)
    if structal is None:
        chart.structal = structal = 0
    primitives = chart.primitives
    if not primitives:
        descs = ((b'char', (1, 0, 0)), (b'short', (2, order, 2)),
                 (b'integer', (4, order, 4)), (b'long', (8, order, 8)),
                 (b'float', (4, order, 4, _binary32)),
                 (b'double', (8, order, 8, _binary64)),
                 (b'text', (1, 0, 0)))  # QnD-specific type for strings
        for name, desc in descs:
            chart.add_primitive(name, desc)
    names = b'short', b'integer', b'long', b'float', b'double'
    prims = [primitives.get(name) for name in names]
    iords = ((1 if p[0].str[0] == '>' else 2) for p in prims[:3])
    ford, dord = [list(range(1, p[0].itemsize+1)) for p in prims[3:]]
    if prims[3][0].str[0] == '<':
        ford = ford[::-1]
    if prims[4][0].str[0] == '<':
        dord = dord[::-1]
    sizes = [p[0].itemsize for p in [prims[2]] + prims]
    header = (b'!<<PDB:II>>!\n' + tobytes((24 + sizes[4] + sizes[5],)) +
              tobytes(sizes) + tobytes(iords) + tobytes(ford) + tobytes(dord))
    fbits, dbits = [p[3] for p in prims[3:]]
    if fbits is None:
        fbits = _binary32
    if dbits is None:
        dbits = _binary64
    header += tobytes(fbits[:7]) + tobytes(dbits[:7])
    header += _byt(fbits[7]) + b'\x01' + _byt(dbits[7]) + b'\x01\n'
    # Record location of chart and symtab addresses, insert dummy values.
    chart.csaddr = len(header)
    header += b'128\x01128\x01\n'
    header += b'\x00' * (128 - len(header))
    f.write(header)
    # Set nextaddr to 128, which is where first data should begin.
    handle.declared(0, None, len(header))


def _iter_nonrec(root):
    for name in root:
        item = root.items[name]
        if item.isleaf():
            yield item
        elif item.isgroup() and '__class__' not in item:
            # Recurse into groups, but not into lists or objects of any type.
            _iter_nonrec(item)


def flusher(f, root):
    handle, chart = root.handle, root.chart
    nextaddr, chart_addr = handle.next_address(both=1)
    blockadds = handle.zero_address(), chart_addr
    attributes = None
    if root.has_attributes:
        attributes = root._detached_subgroup()  # fake group holds attributes
        _dump_attributes(root, attributes, "")
        _, chart_addr = handle.next_address(both=1)
        handle.nextaddr = nextaddr  # reset nextaddr to clobber attributes
    f.seek(chart_addr)
    # Begin by writing just * short integer long float double char to chart,
    # reserving any additional primitives to the PrimitiveTypes extra.
    # This is what yorick does.
    primitives, aligns = chart.primitives, []
    for name in (b'*', b'short', b'integer', b'long', b'float', b'double',
                 b'char'):
        if name == b'*':
            prim = primitives.get(b'char *', primitives[b'long'])
        else:
            prim = primitives[name]
        aligns.append(prim[2])  # prim = stype, dtype, align, desc or None
        size = prim[0].itemsize
        f.write(name + b'\x01' + _byt(size) + b'\x01\n')
    structs = chart.structs
    for name, struct in itemsof(structs):
        size = struct[0].itemsize
        f.write(name + b'\x01' + _byt(size) + b'\x01')
        for mname, (tname, shape) in itemsof(struct[3]):
            if shape:
                shape = b'[' + b','.join(_byt(s) for s in shape) + b']'
            else:
                shape = b''
            f.write(tname + b' ' + mname + shape + b'\x01')
        f.write(b'\n')
    f.write(b'\x02\n')
    # Next comes the symbol table.
    symtab_addr = f.tell()
    blocks = []
    prefix = b''
    for name, item in itemsof(root.items):
        if item.isgroup() or item.islist() == 2:
            prefix = b'/'
            break
    _dump_group(f, prefix, 0, root, blockadds, blocks)
    if attributes is not None:
        _dump_group(f, prefix + b':', 0, attributes, None, None)
        del attributes
    f.write(b'\n')
    # Finally comes the extras section.
    f.write(b'Offset:1\n')  # Default index origin always 1 to match yorick.
    # Alignment: char, *, short, integer, long, float, double, \n
    f.write(b'Alignment:' + tobytes([1] + aligns[:6]) + b'\n')
    f.write(b'Struct-Alignment:' + _byt(chart.structal) + b'\n')
    # Yorick also writes synonym Struct-Align for old yorick bug workaround?
    # Date format "Sun Dec  7 06:00:00 1941" exactly 24 characters.
    date = datetime.today().ctime()
    if not PY2:
        date = date.encode()
    f.write(b'Version:11|' + date + b'\n')  # Yorick writes version 11.
    # PDBLib requires Major-Order: extra before Blocks: extra.
    # Write array dimensions in Fortran order to match yorick and basis.
    f.write(b'Major-Order:102\n')
    hasdirs = bool(prefix)
    f.write(b'Has-Directories:' + _byt(int(hasdirs)) + b'\n')
    f.write(b'Blocks:\n')
    for name, addr, nitems, nblocks in blocks:
        nitems = b' ' + _byt(nitems)
        f.write(name + b'\x01' + _byt(nblocks) + b' ' +
                b' '.join(_byt(a) + nitems for a in addr) + b'\n')
    f.write(b'\x02\n'
            b'Casts:\n\x02\n'
            b'Primitive-Types:\n')
    if hasdirs:
        f.write(b'Directory\x011\x010\x01-1\x01DEFORDER\x01NO-CONV\x01\n')
    for name, prim in itemsof(chart.primitives):
        if name in (b'*', b'char', b'short', b'integer', b'long',
                    b'float', b'double', b'Directory'):
            continue  # skip standard data types, as yorick does
        stype, _, align, desc = prim
        if desc:  # desc = size, order, align, fpbits
            size, order, align = desc[:3]
            fpbits = desc[3] if len(desc) > 3 else None
            if fpbits:
                kind = b'FLOAT\x01' + b'\x01'.join(_byt(i) for i in fpbits)
            else:
                kind = b'FIX'
            if not order:
                order, kind = b'-1\x01DEFORDER\x01', b'NO-CONV'
            elif not hasattr(order, '__getitem__'):
                order = _byt(int(order)) + b'\x01DEFORDER\x01'
            else:
                order = _prim_order(order, size)
        else:
            size, kind, order = stype.itemsize, stype.kind, stype.str[0]
            if kind == 'V':
                order, kind = b'-1\x01DEFORDER\x01', b'NO-CONV'
            else:
                if kind == b'f':
                    order = _prim_order(range(1, size+1) if order == '>' else
                                        range(size+1, 0, -1), size)
                    kind = _FLOAT_KINDS[size]
                else:
                    order = ((b'1' if order == '>' else b'2') +
                             b'\x01DEFORDER\x01')
                    kind = b'FIX'
        f.write(name + b'\x01' + _byt(size) + b'\x01' + _byt(align) +
                b'\x01' + order + kind + b'\x01\n')
    f.write(b'\x02\n'
            b'\n\n')  # PDB extras and file ends with two newlines.
    # Finally, poke chart and symtab addresses into file header.
    csaddr = chart.csaddr
    if csaddr:
        f.seek(csaddr)
        f.write(_byt(chart_addr) + b'\x01' + _byt(symtab_addr) + b'\x01\n')


def _dump_attributes(group, attributes, prefix):
    ignore_first = False
    if group.islist() == 2:
        group = group.parent()
        ignore_first = True
    if group.attrs:
        for aname, value in itemsof(group.attrs):
            aname = prefix + ":" + aname
            attributes._write_attr(aname, value)
    subgroups = []
    for name, item in itemsof(group.items):
        if ignore_first:
            ignore_first = False;
            continue
        name = prefix + name
        if item.isgroup() or item.islist() == 2:
            subgroups.append((item, name))
        else:
            if item.islist():
                item = item.parent()
            shape = item.tsa[1]  # "0" pseudo-attribute for zero length arrays
            if shape and not all(shape):
                attributes._write_attr(name + ":0", array(shape))
            if item.attrs:
                for aname, value in itemsof(item.attrs):
                    aname = name + ":" + aname
                    attributes._write_attr(aname, value)
    for item, name in subgroups:
        _dump_attributes(item, attributes, name + ":")
    return prefix, attributes


def _dump_group(f, prefix, islist, group, blockadds, blocks):
    # First dump the group itself as a bogus Directory object.
    ignore_first = False
    if islist:
        group = group.parent()
        ignore_first = True
    if prefix and (blocks is not None):
        f.write(prefix + b'\x01Directory\x011\x01127\x01\n')
    for name in group:
        if ignore_first:
            # Write a bogus '_' symtab entry to indicate a QList.
            f.write(prefix + b'_\x01char\x011\x01127\x01\n')
            ignore_first = False
            continue
        item = group.lookup(name)
        islist = item.islist()
        name = prefix + (name if PY2 else name.encode('utf8'))
        if item.isgroup() or islist == 2:
            # dump subdirectory
            _dump_group(f, name + b'/', islist, item, blockadds, blocks)
            continue
        # dump leaf (including block variables)
        typename, shape, addr = (item.parent() if islist else item).tsa
        typename = typename[3]  # dtype, stype, align, typename
        if islist:
            shape = (1,) + (shape or ())
        if shape and all(shape):
            size = prod(shape)
            # Set all index origins to 1 to match yorick.
            shape = b'\x01'.join(b'1\x01' + _byt(s)
                                 for s in reversed(shape)) + b'\x01'
        else:
            size = 1
            shape = b''
        if islist:
            a = []
            amin, amax = blockadds
            for ad in addr:
                ad -= amin
                if ad >= 0 and ad < amax:
                    a.append(ad)
            blocks.append((name, a, size, len(a)))
            addr = a[0]
        f.write(name + b'\x01' + typename + b'\x01' + _byt(size) +
                b'\x01' + _byt(addr) + b'\x01' + shape + b'\n')


if PY2:
    def _byt(number): return str(number)  # noqa
    def tobytes(seq): return ''.join(chr(i) for i in seq)  # noqa
else:
    def _byt(number): return str(number).encode()  # noqa
    tobytes = bytes


def _prim_order(order, size):
    return ((b'1' if order[0] <= (size >> 1) else b'2') + b'\x01ORDER\x01' +
            b'\x01'.join(tobytes(order) + b'\x01'))


# Primitive-Types representations of IEEE 754 binary32 and binary64 formats:
_FLOAT_KINDS = {4: b'FLOAT\x0132\x018\x0123\x010\x011\x019\x010\x01127\x01',
                8: b'FLOAT\x0164\x0111\x0152\x010\x011\x0112\x010\x011023\x01'}


# To support qnd attributes, here is the strategy in the flusher:
# 1. Create a separate symbol table for attributes, global attributes
#    having names ':aname' and variable attributes 'vname:aname'.
# 2. Write this entire tree of fake symbols beginning at the first address
#    after the data.
# 3. The chart address becomes the first address after the attributes.
# 4. Write the fake symbol metadata to the symtab after the real symbols.
# 5. Be sure the nextaddr for the file is set to the beginning address
#    of the attribute data, rather than to the chart address.
# On read, we can recognize the first symbol with a : in its name as
# the first attribute, and make sure to reset nextaddr to its address.
# All of the attributes should be read when the file is opened; they are
# part of the metadata like the chart, symtab, and extras.
