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
#
# Given that the byte order is the only difference among its
# predefined i1, i2, i4, i8, f4, and f8 datatypes, there are
# only two possible PDB headers - the big-endian and little-endian
# variants.


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
    header += b'128\001128\001\n'  # bogus chart_addr = symtab_addr = 128
    header += b'\x00' * (128 - len(header))  # null fill to 128 bytes
    f.seek(0)
    f.write(header)
    # Set nextaddr to 128, which is where first data should begin.
    handle = root.handle
    handle.declared(0, _itemsize1, len(header))


_be_header = (b'!<<PDB:II>>!\n\x15\x08\x02\x04\x08\x04\x08'
              b'\x01\x01\x01\x01\x02\x03\x04\x01\x02\x03\x04\x05\x06\x07\x08'
              b'\x20\x08\x17\x00\x01\x09\x00\x40\x0b\x44\x00\x01\x0c\x00'
              b'127\0011023\001\n')
_le_header = (b'!<<PDB:II>>!\n\x15\x08\x02\x04\x08\x04\x08'
              b'\x02\x02\x02\x04\x03\x02\x01\x08\x07\x06\x05\x04\x03\x02\x01'
              b'\x20\x08\x17\x00\x01\x09\x00\x40\x0b\x44\x00\x01\x0c\x00'
              b'127\0011023\001\n')
_itemsize1 = dtype('u1')


def flusher(f, root):
    handle, chart = root.handle, root.chart
    _, chart_addr = handle.next_address(both=1)


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
