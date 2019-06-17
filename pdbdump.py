"""Low level functions to create and write a PDB file."""
from __future__ import absolute_import

from datetime import datetime
import sys

from numpy import array, prod

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
# header[13]  = N    count of single byte values (= 9 + sz_float + sz_double)
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
    chart = root.chart  # root is a PDBGroup
    byteorder = chart.byteorder
    if byteorder is None:
        # PDB byte order is 1 for >, 2 for <
        # Here is a literal way to compute native byte order:
        byteorder = chart.byteorder = int(array([1]).view('u1')[0] + 1)
        # The fact that we're here means we should take the opportunity
        # to initialize the chart to the native set of primitive types.
    header = _be_header if byteorder == 1 else _le_header
    chart.csaddr = len(header)
    header += b'128\001128\001\n'  # bogus chart_addr = symtab_addr = 128
    header += b'\x00' * (128 - len(header))  # null fill to 128 bytes
    f.seek(0)
    f.write(header)
    # Set nextaddr to 128, which is where first data should begin.
    handle = root.handle
    handle.declared(0, None, len(header))


# Given that the byte order is the only difference among its
# predefined i1, i2, i4, i8, f4, and f8 datatypes, there are
# only two possible PDB headers - the big-endian and little-endian
# variants.
_be_header = (b'!<<PDB:II>>!\n\x15\x08\x02\x04\x08\x04\x08'
              b'\x01\x01\x01\x01\x02\x03\x04\x01\x02\x03\x04\x05\x06\x07\x08'
              b'\x20\x08\x17\x00\x01\x09\x00\x40\x0b\x44\x00\x01\x0c\x00'
              b'127\0011023\001\n')
_le_header = (b'!<<PDB:II>>!\n\x15\x08\x02\x04\x08\x04\x08'
              b'\x02\x02\x02\x04\x03\x02\x01\x08\x07\x06\x05\x04\x03\x02\x01'
              b'\x20\x08\x17\x00\x01\x09\x00\x40\x0b\x44\x00\x01\x0c\x00'
              b'127\0011023\001\n')


def flusher(f, root):
    handle, chart = root.handle, root.chart
    _, chart_addr = handle.next_address(both=1)
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
        for tname, shape in itemsof(struct[3]):
            if shape:
                shape = b'[' + b','.join(_byt(s) for s in shape) + b']'
            else:
                shape = b''
            f.write(tname + b' ' + name + shape + b'\x01')
        f.write(b'\n')
    f.write(b'\x02\n')
    # Next comes the symbol table.
    symtab_addr = f.tell()
    prefix = b''
    for name, item in itemsof(root.items):
        if item.isgroup() or item.islist() == 2:
            prefix = b'/'
            break
    blocks = []
    _dump_group(f, prefix, False, root, blocks)
    f.write(b'\n')
    # Finally comes the extras section.
    f.write(b'Offset:0\n')  # Default index origin always 0.
    # Alignment: char, *, short, integer, long, float, double, \n
    f.write(b'Alignment:' + bytes([1] + aligns[:6]) + b'\n')
    f.write(b'Struct-Alignment:' + _byt(chart.structal) + b'\n')
    # Yorick also writes synonym Struct-Align for old yorick bug workaround?
    # Date format "Sun Dec  7 06:00:00 1941" exactly 24 characters.
    date = datetime.today().ctime()
    if not PY2:
        date = date.encode()
    f.write(b'Version:11|' + date + b'\n')  # Yorick writes version 11.
    # PDBLib requires Major-Order: extra before Blocks: extra.
    f.write(b'Major-Order:102\n')  # Array dimensions always in C order.
    hasdirs = bool(prefix)
    f.write(b'Has-Directories:' + _byt(int(hasdirs)) + b'\n')
    f.write(b'Blocks:\n')
    for name, addr, nitems in blocks:
        f.write(name + b'\x01' + _byt(nitems) + b' ' +
                b' '.join(_byt(a) + b' 1' for a in addr) + b'\n')
    f.write(b'\x02\n'
            b'Casts:\n\x02\n'
            b'Primitive-Types:\n')
    if hasdirs:
        f.write(b'Directory\x011\x010\x01-1\x01DEFORDER\x01NO-CONV\x01\n')
    fbias, dbias, fsize, dsize = 127, 1023, 4, 8
    for name, prim in itemsof(chart.primitives):
        if name in (b'float', b'double'):
            fpbits = prim[3]
            if fpbits:
                size, fpbits = fpbits[0], fpbits[3]
                if name.startswith(b'f'):
                    fbias, fsize = fpbits[7], size
                else:
                    dbias, dsize = fpbits[7], size
            continue
        if name in (b'*', b'char', b'short', b'integer', b'long',
                    b'Directory'):
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


def _dump_group(f, prefix, islist, group, blocks):
    # First dump the group itself as a bogus Directory object.
    if islist:
        group = group.parent()
        ignore_first = True
    f.write(prefix + b'\x01Directory\x011\x01127\x01\n')
    for name in group:
        if ignore_first:
            # Write a bogus '_' symtab entry to indicate a QList.
            f.write(prefix + b'_\x01char\x011\x01127\x01\n')
            ignore_first = False
            continue
        item = group[name]
        islist = item.islist()
        name = prefix + name
        if item.isgroup() or islist == 2:
            # dump subdirectory
            _dump_group(f, name + b'/', islist, group, blocks)
            continue
        # dump leaf (including block variables)
        typename, shape, addr = (item.parent() if islist else item).tsa
        typename = typename[3]  # dtype, stype, align, typename
        if islist:
            shape = (1,) + (shape or ())
        if shape:
            size = prod(shape)
            shape = b'\0x01'.join(b'0\x01' + _byt(s)
                                  for s in shape) + b'\0x01'
        else:
            size = 1
            shape = b''
        if islist:
            blocks.append((name, addr, size))
            addr = addr[0]
        f.write(name + b'_\x01' + typename + b'\x01' + _byt(size) +
                b'\x01' + _byt(addr) + b'\x01' + shape + b'\n')


if PY2:
    def _byt(number): return str(number)  # noqa
else:
    def _byt(number): return str(number).encode()  # noqa


def _prim_order(order, size):
    return ((b'1' if order[0] <= (size >> 1) else b'2') + b'\x01ORDER\x01' +
            b'\x01'.join(bytes((i,) for i in order) + b'\x01'))


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
