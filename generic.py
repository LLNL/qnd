"""Generic file or file family open.

"""
from __future__ import absolute_import

from os.path import expanduser, expandvars, abspath, exists, commonprefix
from os.path import getsize
from itertools import count as itertools_count
import glob
import re
import sys

from numpy import int64, uint64  # Needed to work around Windows misfeatures.

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
    **kwargs
       Other keywords.  This opener consumes one item from kwargs:
    nextaddr_mode : int
       Affects setting of nextaddr for families opened with 'a' or 'r+'
       mode.  0 (default) sets nextaddr to the end of the final existing file,
       1 sets nextaddr to 0 (beginning of first file), and 2 sets nextaddr
       to the beginning of the next file after all existing files.

    Returns
    -------
    handle
       A file handle implementing the generic interface, consisting of::

          handle.callbacks(flusher, initializer)
          addr = handle.next_address()  # next unused address
          f = handle.seek(addr)  # return ordinary file handle at addr
          f = handle.open(n)  # open nth file, calling initializer(f)
          handle.flush()  # make file readable, calling flusher(f)
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
        filename = abspath(expanduser(expandvars(filename)))
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
            prefix, suffix = filename[:match.start()], filename
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
                        nums = list(map(int, existing))
                        fmt = '{' + ':0{}d'.format(len(existing[0])) + '}'
                        if all(f == fmt.format(n)
                               for f, n in zip(existing, nums)):
                            existing = nums
                            future = itertools_count(existing[-1] + 1)
                        else:
                            fmt, future = '{}', None
                    elif fmt != '{}':
                        # pattern looked numerical, but matched non-digits
                        fmt, future = '{}', None
                    elif all(len(f) == 1 for f in existing):
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
    return MultiFile(pattern, existing, future, mode, **kwargs), len(existing)


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
    # Force abits to be int64 to work around Windows numpy misfeature in
    # which 42 becomes int32.  We also explicitly convert
    abits = int64(42)

    def __init__(self, pattern, existing, future, mode, **kwargs):
        newfile = not existing
        nextaddr = int64(0)  # Definite type important on Windows.
        nextaddr_mode = (kwargs.pop('nextaddr_mode', 0)
                         if mode.startswith('r+') else 0)
        if newfile:
            try:
                existing = [next(future)]
            except (StopIteration, TypeError):  # TypeError if future is None
                raise IOError("filename specified no filename?")
        else:
            if not isinstance(existing, list):
                existing = list(existing)
            if nextaddr_mode == 0:
                # In a or r+ mode, should begin with a guess at nextaddr,
                # which is the end of the last existing file.
                # This may be modified by the non-generic caller (see pdbf).
                i = len(existing) - 1
                nextaddr = int64(getsize(pattern.format(existing[i])))
                nextaddr |= int64(i) << self.abits
            elif nextaddr_mode == 2:
                nextaddr = int64(i) << self.abits
        current = 0
        self.f = open(pattern.format(existing[current]), mode)
        self.state = [mode, pattern, existing, current, future]
        self._callbacks = None, Ellipsis if newfile else None
        self.nextaddr = nextaddr

    def callbacks(self, flusher, initializer):
        """set callback function that flushes file metadata"""
        newfile = self._callbacks[1] is Ellipsis  # set in __init__ only
        self._callbacks = flusher, initializer
        if newfile:
            initializer(self.f)

    def filename(self, n=None):
        """current or n-th existing filename in family"""
        _, pattern, existing, current, _ = self.state
        if n is None:
            n = current
        return pattern.format(existing[n])

    def filemode(self):
        return self.state[0]

    def open(self, n):
        """open n-th file of family"""
        mode, pattern, existing, current, future = self.state
        if n == current:
            return self.f
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
            # Do not catch StopIteration here.  Instead handle in caller
            # of next_address(newfile=1).
            if future is None:
                raise StopIteration
            member = next(future)
            if not mode.startswith('w'):
                mode = 'w+b'
        if isnew:
            self.flush()
        f = open(pattern.format(member), mode)
        self.f.close()
        self.state[3] = n
        self.f = f
        # Starting nextaddr at 0 is correct only when isnew.
        # In other cases we leave it set to its most recent known value
        # in a different file.  This assures nextaddr is in last file.
        if isnew:
            existing.append(member)
            self.nextaddr = int64(n) << self.abits
            initializer = self._callbacks[1]
            if initializer not in (None, Ellipsis):
                # Note that initializer may call open recursively, but if
                # so caller must take care to ensure that such a recursive
                # call does not reach this point.
                initializer(f)
        return f

    def flush(self):
        """flush metadata and ordinary file buffers"""
        mode = self.state[0]
        if mode.startswith('r') and '+' not in mode:
            return  # quick no-op for read-only files
        existing, current = self.state[2:4]
        if current + 1 == len(existing):
            flusher = self._callbacks[0]
            if flusher is not None:
                i, nextaddr = self.split_address(self.nextaddr)
                if nextaddr == 0 and i == current + 1:
                    return  # handle special case nextaddr_mode==2
                if i != current:
                    raise AssertionError("(BUG) impossible current file value")
                f = self.f
                addr = f.tell()
                f.seek(nextaddr)
                flusher(f)
                f.seek(addr)
        self.f.flush()

    def close(self):
        """flush and close the current file"""
        self.flush()
        self.f.close()

    def current_file(self):
        """Index of current file in family, argument to open method."""
        return self.state[3]

    def zero_address(self, n=None):
        """multifile address of first byte in current or n-th file"""
        return int64(self.state[3] if n is None else n) << self.abits

    # This is not used by pdbf module, but provide it anyway.
    def tell(self):
        """return current multi-file address"""
        current = int64(self.state[3])
        abits = self.abits
        one = int64(1)
        mask = (one << abits) - one
        addr = int64(self.f.tell())
        if addr & ~mask:
            raise IOError("file too large for {} bit address".format(abits))
        return (current << abits) | addr

    def seek(self, addr):
        """seek to multi-file address, opening alternate file if needed"""
        i, addr = self.split_address(addr)
        f = self.open(i)
        f.seek(addr)
        return f

    def split_address(self, addr):
        """return file index, address for a multifile address"""
        # There is a serious long standing bug in numpy type promotion rules
        # which prevents uint64 from being useful when combined with any other
        # integer type -- numpy will promote uint64 to float64 in a stupid
        # attempt to find a signed type that can hold the result.
        # This is completely different from C type promotion rules, in
        # which signed gets promoted to unsigned (also a bad idea).
        addr = int64(addr)  # even on 64 bit Windows, addr can be int32
        one = int64(1)
        abits = self.abits
        mask = (one << (int64(64) - abits)) - one
        i = (addr >> abits) & mask
        mask = (one << abits) - one
        return i, addr & mask

    def next_address(self, both=False, newfile=False):
        """next unused multi-file address, or None if newfile cannot create"""
        # Use special value of nextaddr as implicit newfile flag.
        nfiles = len(self.state[2])
        if newfile or self.nextaddr == int64(nfiles) << self.abits:
            try:
                self.open(nfiles)
            except StopIteration:
                # Signal caller that we have run out of filenames.
                return (None, None) if both else None
        nextaddr = int64(self.nextaddr)
        if both:
            one = int64(1)
            mask = (one << self.abits) - one
            return nextaddr, nextaddr & mask
        return nextaddr

    def declared(self, addr, dtype, nitems):
        """declare that array has been declared, maybe update next_address"""
        addr = int64((nitems if dtype is None else nitems * dtype.itemsize)
                     + addr)
        nextaddr = int64(self.nextaddr)
        if uint64(addr) > uint64(nextaddr):
            self.nextaddr = addr
