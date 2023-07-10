from __future__ import absolute_import

from numbers import Integral

from numpy import cumprod


def leading_args(args, shape):
    if not args:
        return args, shape, 0
    stride = cumprod((1,) + shape[::-1])[-2::-1]
    offset = 0
    args, shape = list(args), list(shape)
    # First trim any fixed integer values off the leading args.
    # This requires shifting the address by the corresponding offset,
    # and reducing the effective shape of the stored data.  If the
    # leading remaining dimension is a slice, we shrink the shape around
    # the slice, and modify the slice correspondingly.
    for i, a in enumerate(args):
        n = shape[i]  # failure here means more simple args than dimensions
        if isinstance(a, Integral):
            if a < 0:
                a += n
            if a < 0 or a >= n:
                raise IndexError("scalar index out of range")
            offset += a * stride[i]
            continue
        if isinstance(a, slice):
            j0, j1, di = _slice_endpoints(a, n)
            # now x[i0:i1:di] == x[j0:j1][::di]
            if j1 == j0:
                j0 = j1 = 0
            if j0:
                offset += j0 * stride[i]
            shape[i] = j1 - j0
            args[i] = slice(None, None, di) if di != 1 else slice(None)
        if i:
            args, shape = args[i-1:], shape[i-1:]
        break
    else:
        # All args were scalar int.
        args, shape = [], shape[len(args):]
    # Next, trim any full length trailing dimensions.  This is used to
    # detect whether or not a write operation requires a read-modify-write
    # cycle (presuming yes if the remaining args list is not empty).
    # This program may be too difficult to carry out if any argument
    # corresponds to multiple dimensions or is otherwise complex, so
    # we need to work out the args <--> shape correspondence first.
    nargs, ndim = len(args), len(shape)
    iellipsis = None
    for i, a in enumerate(args):
        if a is Ellipsis:
            iellipsis = i
        elif not isinstance(a, slice):
            break
    else:
        # We may be able to remove some trailing dimensions.
        soffset = 0 if iellipsis is None else ndim - nargs
        if soffset >= -1:
            nremove = 0
            for i, a in reversed(list(enumerate(args))):
                if i == iellipsis:
                    soffset = 0
                    # Note that if ellipsis is present and not removed,
                    # we cannot remove anything.
                    iellipsis = None
                else:
                    n = shape[i+soffset]
                    if _slice_endpoints(a, n) != (0, n, 1):
                        break
                nremove += 1
            if nremove and iellipsis is None:
                args = args[:-nremove]
    return tuple(args), tuple(shape), offset


def _slice_endpoints(s, n):
    i0, i1, di = s.indices(n)
    # Although i0 is the first point, i1 may not be one beyond the
    # nominal final point.  Furthermore, if di < 0, i1 will be less
    # than i0.  Unscramble all of this to find j0 < j1 so that:
    #   x[i0:i1:di] == x[j0:j1][::di]
    if di > 1:
        i1 = ((i1 - i0 - 1)//di)*di + i0 + 1
    elif di < 0:
        i1, i0 = i0 + 1, ((i1 - i0 + 1)//di)*di + i0
    return i0, i1, di
