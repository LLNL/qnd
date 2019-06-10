"""Generic file or file family open.

"""
from __future__ import absolute_import

from os.path import expanduser, expandvars, abspath, exists, commonprefix
from itertools import count as itertools_count
from string import ascii_lowercase, ascii_uppercase
import glob
import re
import sys

PY2 = sys.version_info < (3,)
if PY2:
    range = xrange
    chr = unichr
else:
    basestring = str
_glob_group = re.compile(r'(\*|\?|\[.+?\])+')
_glob_digit = re.compile(r'(?:\*|\?|\[0-9\])')
_glob_ranges = re.compile(r'[^-[](?:-[^-])?')
_digits_format = re.compile(r'%(?:0(\d+))?d')


def opener(filename, mode, **kwargs):
    """Generic file or file family opener.

    Parameters
    ----------
    filename : str
       Name of file to open.  See notes below for family conventions.
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
       Other keywords.

    Returns
    -------
    handle
       A file handle implementing the generic interface, consisting of::

          handle.callback(flusher)  # call flusher(handle) to flush
          addr = handle.next_address()  # next unused address
          f = handle.seek(addr)  # return ordinary file handle at addr
          handle.flush()  # make file readable, calling flusher(handle)
          # flush() restores next_address to value on entry
          handle.close()  # flush if necessary, then close
    nexisting : int
       Number of existing paths matching `filename`.

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

    A single filename may also contain a ``%d`` or ``%0nd`` print format
    directive, which will be converted to the corresponding number of
    ``[0-9]`` glob patterns.

    """
    isstr = isinstance(filename, basestring)
    if isstr:
        filename = expanduser(expandvars(filename))
        match = _digits_format.search(filename)
        if match:  # %0nd --> [0-9][0-9]... n times
            n = int(match.group(1) or '1')  # %d --> [0-9] once
            filename = filename.replace(match.group(0), ''.join(['[0-9]']*n))
        match = _glob_group.search(filename)
        if match:
            existing = glob.glob(filename)
        else:
            existing = [f for f in [filename] if exists(f)]
    else:
        match = None
        filename = [expanduser(expandvars(f)) for f in filename]
        existing = [f for f in filename if exists(f)]
        if len(existing) < len(filename):
            for g, f in zip(existing, filename):
                if f == g:
                    continue
                raise IOError("intermediate file {} missing".format(f))
        if not filename:
            raise IOError("sequence of filenames is empty")
    mode = mode.lower()
    if existing:
        existing.sort(key=lambda path: (len(path), path))
        if mode.startswith('w-'):
            raise IOError("protecting existing file {}".format(existing[0]))
    elif mode.startswith('r'):
        raise IOError("no file matches {}".format(filename))
    if mode.startswith('r'):
        mode = 'r+b' if mode.startswith('r+') else 'rb'
    elif mode.startswith('a'):
        mode = 'r+b' if existing else 'w+b'
    elif mode.startswith('w'):
        existing = []  # ignore any existing files
        mode = 'w+b'
    else:
        raise IOError("open mode {} not understood".format(mode))
    # Compute:
    # pattern = pattern containing {...} if either existing or future
    # existing = list of existing {...} items (int or str) for pattern
    # future = iterable yielding future {...} items (int or str) for pattern
    future = None
    if match:
        if '+' in mode:
            prefix = filename[:match.start()]
            predictable = 2
            while match:
                predictable <<= 1
                suffix = suffix[match.end():]
                match = _glob_group.search(suffix)
            p, s = len(prefix), -len(suffix) if suffix else None
            existing = [f[p:s] for f in existing]
            if predictable:
                # With a single wildcard group, we may be able to predict
                # future names in the family.
                # We handle two cases:
                # 1. a sequence of *, ?, [0-9] we guess means decimal numbers
                # 2. a single range like [a-z] we take as sequence of chars
                fmt = '{}'
                pat = filename[p:s]
                chunks = _glob_digit.findall(pat)
                if ''.join(chunks) == pat:
                    nast = chunks.count('*')
                    ndig = len(chunks) - nast
                    fmt = '{' + ':0{}d'.format(3*nast + ndig) + '}'
                    future = itertools_count(0)
                else:
                    future = iter(_expand_ranges(pat))
                if existing:
                    if all(f.isdigit() for f in existing):
                        # existing matches are all decimal numbers
                        nums = map(int, existing)
                        fmt = '{' + ':0{}d'.format(len(existing[0])) + '}'
                        if all(f==fmt.format(n)
                               for f, n in zip(existing, nums)):
                            existing = nums
                            future = itertools_count(existing[-1] + 1)
                        else:
                            fmt, future = '{}', None
                    elif fmt != '{}':
                        # pattern looked numerical, but matched non-digits
                        fmt, future = '{}', None
                    elif all(len(f)==1 for f in existing):
                        # existing matches all non-digit single characters
                        final = existing[-1]
                        for f in future:
                            if f == final:
                                break
            pattern = prefix + fmt + suffix
        else:
            filename = existing
            isstr = False
    elif isstr:
        pattern = '{}'
        if not existing:
            future = iter([filename])
    if not isstr:
        prefix = commonprefix(filename)
        if len(filename) > 1:
            suffix = commonprefix([f[::-1] for f in filename])[::-1]
            if suffix == prefix:
                suffix = ''  # all filenames identical (not an error??)
        else:
            suffix = ''
        n = len(prefix)
        m = -len(suffix) if suffix else None
        existing = [f[n:m] for f in existing]
        future = iter([f[n:m] for f in filename[len(existing):]])
        pattern = prefix + '{}' + suffix
    return MultiFile(pattern, existing, future, mode), len(existing)


def _expand_ranges(pat):
    chunks = _glob_ranges.findall(pat)
    if '[{}]'.format(''.join(chunks)) != pat:
        return None
    for i, c in enumerate(chunks):
        if len(c) != 3:
            continue
        c2 = ord(c[2])
        if c2 >= 256:
            return None
        chunks[i] = ''.join(map(chr, range(ord(c[0]), c2+1)))
    chunks = ''.join(chunks)
    return chunks if len(chunks) == len(set(chunks)) else None


# MultiFile implements the basic file handle methods used by pdbf.
class MultiFile(object):
    """A binary file or family of binary files."""
    # Main purpose is to combine address spaces in a file family, so
    # "address" is (file_index | byte_address) packed into 64 bit int.
    # Default number of address bits is 42, with top 22 bits reserved
    # for file index.  This allows for 4 million 4 terabyte files, which
    # should be adequate for most purposes as of 2019.  You can always
    # override abits in an instance if necessary.
    # The address -1 is reserved to mean "no address", as in NULL pointer.
    abits = 42

    def __init__(self, pattern, existing, future, mode):
        if not existing:
            try:
                existing = [next(future)]
            except StopIteration:
                raise IOError("filename specified no filename?")
        elif not isinstance(existing, list):
            existing = list(existing)
        current = 0
        self.f = open(pattern.format(existing[current]), mode)
        self.state = [mode, pattern, existing, current, future]
        self._callbacks = None, None
        self.nextaddr = 0

    def callbacks(self, flusher, initializer):
        """set callback function that flushes file metadata"""
        self._callbacks = flusher, initializer

    def filename(self, n=None):
        """current or n-th existing filename in family"""
        _, pattern, existing, current, _ = self.state
        if n is None:
            n = current
        return pattern.format(existing[n])

    def zero_address(self, n=None):
        """multifile address of first byte in current or n-th file"""
        return self.state[3] << self.abits

    def open(self, n):
        """open n-th file of family"""
        mode, pattern, existing, current, future = self.state
        writeable = mode.startswith('w') or '+' in mode
        isnew = n == len(existing)
        if n < len(existing):
            if writeable:
                mode = 'r+b'
            member = existing[n]
        elif not isnew:
            raise IOError("cannot open file {} in family of {}"
                          "".format(n, len(existing)))
        elif not writeable:
            raise IOError("cannot create new file in read-only family")
        else:
            try:
                member = next(future)
            except StopIteration:
                raise IOError("no rule to compute family name following {}"
                              "".format(pattern.format(existing[-1])))
            if not mode.startswith('w'):
                mode = 'w+b'
        f = open(pattern.format(member), mode)
        if isnew:
            existing.append(member)
        self.state[3] = n
        self.f = f
        self.nextaddr = n << self.abits
        initializer = self._callbacks[1]
        if initializer is not None:
            initializer(f)
        return f

    def flush(self):
        """flush metadata and ordinary file buffers"""
        mode = self.state[0]
        if mode.startswith('r') and '+' not in mode:
            return  # quick no-op for read-only files
        flusher = self._callbacks[0]
        if flusher is not None:
            addr = self.tell()
            flusher(self.seek(self.next_address()))
            self.seek(addr)
        self.f.flush()

    def close(self):
        """flush and close the current file"""
        self.flush()
        self.f.close()

    # This is not used by pdbf module, but provide it anyway.
    def tell(self):
        """return current multi-file address"""
        current = self.state[3]
        abits = self.abits
        mask = (1 << abits) - 1
        addr = self.f.tell()
        if addr & ~mask:
            raise IOError("file too large for {} bit address".format(abits))
        return (current << abits) | addr

    def seek(self, addr):
        """seek to multi-file address, opening alternate file if needed"""
        current = self.state[3]
        abits = self.abits
        mask = (1 << (64 - abits)) - 1
        c = (addr >> abits) & mask
        if c != current:
            self.close()
            f = self.open(c)
        else:
            f = self.f
        mask = (1 << abits) - 1
        f.seek(addr & mask)
        return f

    def next_address(both=False, newfile=False):
        """next unused multi-file address, or None if newfile cannot create"""
        if newfile:
            try:
                self.open(len(self.state[2]))
            except StopIteration:
                return None
        nextaddr = self.nextaddr
        if both:
            mask = (1 << self.abits) - 1
            return nextaddr, nextaddr & mask
        return nextaddr

    def declared(self, addr, dtype, nitems):
        """declare that array has been declared, maybe update next_address"""
        addr = dtype.itemsize * nitems + addr
        nextaddr = self.nextaddr
        if addr > nextaddr:
            self.nextaddr = nextaddr
