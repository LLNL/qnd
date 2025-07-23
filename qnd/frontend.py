"""Quick and Dirty, a high level file access wrapper.

The QnD interface looks like this::

   f = format_specific_open(filename, mode)  # e.g. openh5
   var = f.varname  # read var
   var = f['varname']
   var = f.get('varname', default)

   f.var = var_value  # declare and write var
   f['var'] = something
   f.var = dtype, shape  # declare var without writing
   f.update({...vars...}, another_var=value, ...)
   f.grpname = {...vars...}  # declare a subgroup and some members

   if name in f: do_something
   varnames = list(f)
   for name in f: do_something
   for name, var in f.items(): do_something

   g = f.grpname
   f = g.root()  # Use the root method to get the top-level QGroup.
   f.close()  # important if you have written to f

Generally, a QnD QGroup like `f` in the example behaves like a dict.
However, you may also reference variables or subgroups as if they were
attributes.  Use attributes to access variables when you know the
variable name.  In short, use square brackets when the variable name
is the value of an expression.  (QnD will remove a single trailing
underscore from any attribute reference, so you can use ``f.yield_``
for ``f['yield']`` or ``f.items_`` for ``f['items']``.)  The adict
module has an `ADict` class and a `redict` function to produce ordinary
in-memory dict objects with their items accessible as attributes with
the same rules.  You can read a whole file (or a whole subgroup) like
this::

   ff = f(2)

The optional `2` argument is the auto-read mode flag.  By default, the
auto-read mode flag is set to 1, which causes ``f.varname`` to read an
array variable and return its value, but to simply return a QGroup
object (like `f`) if the name refers to a subgroup.  When the `auto`
flag equals `2`, any subgroups are read recursively, and their values
become ADict instances.  (QnD also supports QList variables, and
``auto=2`` mode returns those as python list instances.)

The ``items()`` method also accepts an optional `auto` argument to
temporarily change auto-read mode used for the iteration.

You can turn auto-read mode off by setting the `auto` flag to `0`.  In
this mode, referencing a variable returns a QLeaf instance without
reading it.  This enables you to query a variable without reading it.
You can also do that by retrieving the attributes object::

   with f.push():  # Use f as a context manager to temporarily change modes.
       f.auto(0)  # Turn off auto-read mode.
       v = f.varname
   value = v()  # Read a QLeaf by calling it...
   value = v[:]  # ...or by indexing it.
   v(value)  # Write a QLeaf by calling it with an argument...
   v[:] = value  # ...or by setting a slice.
   v.dtype, v.shape, v.size, v.ndim  # properties of the QLeaf v
   # An alternate method which pays no attention to auto mode:
   va = f.attrs.varname  # Get attributes of varname.
   va.dtype, va.shape, va.size, va.ndim  # Built-in pseudo-attributes.
   # You can use va to get or set real attributes of varname as well:
   units = va.units  # retrieve units attribute
   va.centering = 1  # set centering attribute

When you call a QGroup like `f` as a function, you may also pass it a
list of variable names to read only that subset of variables.  With
auto-read mode turned off, this results in a sort of "casual subgroup"::

   g = f(0, 'vname1', 'vname2', ...)
   h = f(1, 'vname1', 'vname2', ...)
   ff = f(2, 'vname1', 'vname2', ...)

Here, g is an ADict containing QLeaf and QGroup objects, with nothing at
all read from the file, while h is and ADict containing ndarray and QGroup
objects, while ff is an ADict containing ndarray and ADict objects, with
no references at all to `f`.

If you want to use `f` as a context manager in the manner of other
python file handles, so that the file is closed when you exit the with
statement, just do it::

   with openh5(filename, "a") as f:
       do_something(f)
   # f has been properly flushed and closed on exit from the with.

------

QnD also supports old netCDF style UNLIMITED dimensions, and their
equivalents in HDF5.  Unlike the netCDF or HDF5 interface, in QnD the
first (slowest varying) dimension of these arrays maps to a python
list, so we regard the entire collected variable as a list of
ndarrays.  The netCDF record number is the index into the list, while
any faster varying dimensions are real ndarray dimensions.  This
subtle difference in approach is more consistent with the way these
variables are stored, and also generalizes to the fairly common case
that the array dimensions -- often mesh dimensions -- change from one
record to the next.

To write records using QnD, turn on "recording mode"::

   f.recording(1)  # 0 for off, 2 for generalized records
   f.time = 0.
   f.x = x = arange(10)
   f.time = 0.5
   f.x = x**2

Ordinarily, when you set the value of ``f.time`` or ``f.x``, any
previous value will be overwritten.  But in recording mode, each time
you write a variable, you create a new record, saving the new value
without overwriting the previous value.  If you want all record
variables to have the same number of records, you need to be sure
you write them each the same number of times.  One way to do that is
to use the update function rather than setting them one at a time::

   record = ADict()
   record.time, record.x = 0., arange(10)
   f.recording(1)
   f.update(record)
   record.time, record.x = 0.5, record.x**2
   f.update(record)

You cannot change a variable from not having records to having records
(or from recording mode 1 to recording mode 2); the recording mode in
force when a variable was first declared determines if and how all
future write operations behave.

Reading back record variables introduces "goto mode".  Initially, goto
mode is off or None, so that reading a record variable gets the whole
collection of values as a QList, or as an ordinary python list if
auto mode is on::

   f.goto(None)  # explicitly turn off goto mode
   f.auto(2)
   times = f.time  # python list of f.time values
   xs = f.x  # python list of f.x arrays
   f.auto(0)
   time = f.time  # QList for the collection of time values
   nrecords = len(time)

On the other hand, with goto mode turned on, the fact that `time` and `x`
are record variables disappears, so that your view of ``f.time`` and
``f.x`` matches what it was when you recorded them.  You use the goto
function to set the record::

   f.goto(0)  # first record is 0, like any python list
   t = f.time  # == 0.
   f.goto(1)  # set to second record
   t = f.time  # == 0.5
   x = f.x  # == arange(10)**2
   f.goto(-1)  # final record, negative index works like any python list
   # You can also pass a keyword to goto, which can be the name of any
   # scalar record variable, to go to the record nearest that value.
   f.goto(time=0.1)  # will select record 0 here
   current_record = f.goto()  # goto() returns current record number

   for r in f.gotoit(): do_something  # f.goto(r) is set automatically

Note the ``gotoit()`` method returns an iterator over all records,
yielding the record number for each pass, and setting the goto record
for each pass automatically.  You can use ``f.push()`` in a with
statement to temporarily move to a different record.

If you set the recording mode to `2`, the record variables need not
have the same shape or same type from one record to the next (indeed,
they can be a subgroup on one record and an array on another).  This
cannot be represented as an UNLIMITED array dimension in an HDF5 or
netCDF file, so the QList variable in QnD will become an HDF5 group in
this case, where variable names in the group are _0, _1, _2, and so on
for QList element 0, 1, 2, and so on (plus a hidden element _ which
identifies this group as a list when it is empty).  You can create a
QList of this general type without using recording or goto mode at
all::

   f.recording(0)  # Turn off recording and goto modes.
   f.goto(None)
   f.varname = list  # Make an empty QList
   ql = f.varname
   ql.append(value0)
   ql.extend([value1, value2, ...])
   var = ql[1]  # retrieves value1
   nelements = len(ql)  # current number of elements (also works for QGroup)
   ql.auto(0)  # a QList has auto mode just like a QGroup
   for var in ql: do_something  # var depends on ql auto mode setting

------

"""
from __future__ import absolute_import

# The three kinds of objects we support here are:
# 1. QGroup --> dict with str keys
# 2. QList --> list
# 3. QLeaf --> ndarray with dtype.kind in buifcS, with U encoded as S
#              and V handled as recarray.
# We attempt to store arbitrary python objects as a QGroup with member
#   __class__ = 'modulename.classname' (or just classname for builtin)
#   remainder of QGroup is the instance __dict__, unless __class__ has a
#   __setstate__, in which case argument to that method stored in the
#   __setstate__ variable.
#   If __class__ has a __getnewargs__, result is written sub-QGroup with
#   _0, _1, ..., which will be passed to the class constructor -- otherwise
#   the class starts empty and neither __new__ nor __init__ is called.
#   List or tuple objects not distinguished, becoming QList items.
#   Dict objects with non-text keys stored with __class__ = 'dict' and
#   members _0, _1, _2, etc., where even item is key and following odd item
#   is corresponding value.
# We ask the backend to support the variable value None as a QLeaf in
# addition to the arrays, if possible.  Arrays with zero length dimensions
# should also be supported if possible.
#
# This qnd module also provides a low level QnDList implementation of the
# QList in terms of the backend QGroup (recording=2) and QLeaf (recording=1)
# implementations, for backends which do not support a native list type.
# The convention is that a generic QList is a QGroup with a blank member _
# and members _0, _1, _2, etc.  If the backend supports QLeaf arrays with an
# UNLIMITED leading dimension, these can also be presented as QList
# variables by the QnD API.

# Backend object methods assumed here:
# qnd_group methods: close(), flush(), root()
#   isgroup() -> 1, islist() -> 0, isleaf() -> 0
#   __len__, __iter__ returns names
#   lookup(name) -> None if not found
#   declare(name, dtype, shape, unlim)  dtype can be dict, list, or None
#   attget(vname) --> variable attributes, vname='' for group attributes
#   attset(vname, aname, dtype, shape, value) --> variable attributes
# qnd_list methods: root()
#   isgroup() -> 0, isleaf() -> 0
#   islist() -> 1 if this is UNLIMITED dimension, -> 2 if anonymous group
#   __len__, __iter__ returns unread elements
#   index(i) -> None if i out of range
#   declare(dtype, shape)  dtype can be dict, list, or None
# qnd_leaf methods: root()
#   isgroup() -> 0, islist() -> 0, isleaf() -> 1
#   query() -> dtype, shape, sshape  (None, (), () for None)
#   read(args)
#   write(value, args)

import sys
from weakref import proxy, ProxyTypes
from importlib import import_module
import re
# Major change in array(x) function semantics when x is a list of ragged
# arrays:  This now generates a DeprecationWarning as of 1.19, and
# presumably an exception for some future numpy (beyond 1.21).  See the
# _categorize function below for the workaround.
from warnings import catch_warnings, simplefilter

try:
    from numpy import VisibleDeprecationWarning
except ImportError:
    VisibleDeprecationWarning = DeprecationWarning
from numpy import (dtype, asarray, asanyarray, arange, interp, where, prod,
                   ndarray)
from numpy.core.defchararray import encode as npencode, decode as npdecode

from .adict import ItemsAreAttrs, ADict

PY2 = sys.version_info < (3,)
if PY2:
    range = xrange
else:
    basestring = str
_NOT_PRESENT_ = object()
_us_digits = re.compile(r"^_\d*$")


class QGroup(ItemsAreAttrs):
    """Group of subgroups, lists, and ndarrays.

    You reference QGroup items by name, either as ``qg['name']`` like a
    dict item, or equivalently as ``qg.name`` like an object attribute.
    Use ``[]`` when the item name is an expression or the contents of
    a variable; use ``.`` when you know the name of the item.  You can
    use ``[]`` or ``.`` to both get and set items in the QGroup.  To
    read the entire group into a ADict, call it like a function, ``qg()``;
    you may supply a list of names to read only a subset of items.
    A QGroup acts like a dict in many ways::

       if 'name' in qg: do_something
       for name in qg: do_something
       item_names = list(qg)  # qg.keys() exists but is never necessary
       for name, item in qg.items(): do_something
       qg.update({name0: val0, ...}, [(name1, val1), ...], name2=val2, ...)
       value = qg.get('name', default)

    A QGroup has several possible states or modes:

    1. Recording mode, turned on by ``qg.recording(1)`` and off by
       ``qg.recording(0)``, affects what happens when you set group items.
       With recording mode off, setting an item to an array creates the
       item as an array if its name has not been used, or otherwise writes
       its new value, requiring it be compatible with the dtype and shape
       of the previous declaration.  With recording mode on, setting an
       item for the first time creates a QList and sets its first element
       to the given value, and subsequently setting that item appends the
       given value to the existing QList.  There is also a recording mode
       ``qg.recording(2)`` in which subsequent values need not match the
       dtype or shape of the first item.  You may not switch recording
       modes for a given item; the mode in effect when an item is first
       created governs the behavior of that item.
    2. Goto mode, in which you set a current record with ``qg.goto(rec)``.
       Any item you retrieve or query which is a QList retrieves or queries
       the element with 0-origin index ``rec`` instead of the whole QList.
       You turn off goto mode with ``qg.goto(None)``.  There is also a
       ``qg.gotoit()`` function which returns an iterator over all the
       records (generally the longest QList in ``qg``).
    3. Auto mode, turned on by ``qg.auto(1)`` and off by ``qg.auto(0)``,
       in which getting any item reads and returns its value, rather than
       a QLeaf object.  There is also a ``qg.auto(2)`` mode in which
       the auto-read feature applies to any QGroup or QList (if goto mode
       is off) items recursively.

    A QGroup has `push` and `drop` methods which can be used to save and
    restore all its modes.  The `drop` method is called implicitly upon
    exit from a with statement, so you can use the QGroup as a context
    manager::

       with openh5('myfile.h5', 'a') as qg:
           do_something(qg)
           with qg.push():
               qg.goto(rec)
               do_something_else(qg)
           # qg restored to goto mode state before with.
           do_even_more(qg)
       # qg flushed and closed upon exit from with clause that has no
       # no corresponding push

    Attributes
    ----------
    islist
    isleaf
       Always 0.
    isgroup
       Always 1.
    dtype
       Always ``dict``, the builtin python type.
    shape
    ndim
    size
    sshape
       Always None.

    """
    __slots__ = "_qnd_group", "_qnd_state", "_qnd_cache", "__weakref__"
    isgroup = 1
    islist = isleaf = 0
    dtype, shape, ndim, size, sshape = dict, None, None, None, None

    def __init__(self, item=None, state=None, auto=None, recording=None,
                 goto=None):
        object.__setattr__(self, "_qnd_group", item)
        object.__setattr__(self, "_qnd_state",
                           QState() if state is None else QState(state))
        object.__setattr__(self, "_qnd_cache", None)
        state = self._qnd_state
        if auto is not None:
            state.auto = int(auto)
        if recording is not None:
            state.recording = int(recording)
        if goto is not None:
            state.goto = int(goto)

    def recording(self, flag):
        """Change recording mode for this QGroup.

        With recording mode off, writing to a variable overwrites that
        variable.  With recording mode on, new variables are declared as
        a QList and subsequent write operations append a new element to
        this QList instead of overwriting any previously stored values.
        In netCDF parlance, variables declared in recording mode are
        record variables.  Writing to a variable declared when recording
        mode was off will always overwrite it; once declared, you cannot
        convert a variable to a QList simply by turning on recording mode.

        See goto mode for handling record variable read operations.

        A `flag` value of 0 turns off recording mode.  A `flag` of 1 turns
        on recording mode, utilizing a trailing UNLIMITED array dimension
        in netCDF or HDF5 parlance, which promises that all values written
        will have the same dtype and shape.  A `flag` of 2 places no
        restrictions on the dtype or shape of the QList elements; such
        an unrestricted QList resembles an anonymous QGroup.

        """
        self._qnd_state.recording = int(flag)

    def goto(self, record=_NOT_PRESENT_, **kwargs):
        """Set the current record for this QGroup, or turn off goto mode.

        Pass `record` of None to turn off goto mode, so that QList
        variables appear as the whole QList.  Setting an integer `record`
        makes any QList variable appear to be the specified single
        element.  A `record` value may be negative, with the usual python
        interpretation for a negative sequence index.  If different QList
        variables have different lengths, the current `record` may be
        out of range for some variables but not for others.  (Hence using
        goto mode may be confusing in such situations.)

        Note that you can temporarily set goto mode using a with clause.

        This `goto` method also accepts a keyword argument instead of a
        `record` number.  The keyword name must match the name of a
        QList variable in this QGroup, whose vaules are scalars.  This
        will set `record` to the record where that variable is nearest
        the keyword value.  Thus, ``goto(time=t)`` selects the record
        nearest `time` t.

        As a special case, you can get the current record number by calling
        `goto` with neither a `record` nor a keyword::

            current_record = qg.goto()

        """
        if kwargs:
            if record is not _NOT_PRESENT_:
                raise TypeError("either use keyword or record index")
            if len(kwargs) != 1:
                raise TypeError("only one keyword argument accepted")
            name, val = list(kwargs.items())[0]
            records, values = self._qnd_goto_recs(name)
            val = float(val)
            n = values.size
            if n > 1:
                # result of interp is float scalar in older numpy versions
                # rather than numpy.float64, cannot use astype
                i = int(interp(val, values, arange(n) + 0.5))
            else:
                i = 0
            record = records[min(i, n-1)]
        elif record is _NOT_PRESENT_:
            return self._qnd_state.goto
        elif record is not None:
            record = int(record)
        self._qnd_state.goto = record

    def _qnd_goto_recs(self, name):
        cache = self._qnd_cache
        values = cache.get(name) if cache else None
        if values is None:
            item = self._qnd_group.lookup(name)
            if item is not None and item.islist():
                with self.push():
                    self.goto(None)
                    self.auto(2)
                    values = self[name]
                values = asarray(values, float)
                if values.ndim != 1 or values.size < 1:
                    values = None
            if values is None:
                raise TypeError("{} is not scalar record variable"
                                "".format(name))
            values = _monotonize(values)
            if not cache:
                cache = {}
                object.__setattr__(self, "_qnd_cache", cache)
            cache[name] = values
        return values  # returned by _monotonize

    def gotoit(self, name=None):
        """Iterate over goto records, yielding current record.

        Optional `name` argument is the name of a `goto` method keyword,
        which may implicitly remove records corresponding to non-monotonic
        changes of that variable.  If `name` is a decreasing variable,
        the record order will be reversed.

        As a side effect, the current record of this QGroup will be set
        during each pass.  If the loop completes, the original goto state
        will be restored, but breaking out of the loop will leave the
        goto record set.

        """
        if name is not None:
            records, _ = self._qnd_goto_recs(name)
        else:
            # scan through all variables to find largest recrod count
            nrecords = 0
            for name in self._qnd_group:
                item = self._qnd_group.lookup(name)
                if item.islist():
                    n = len(item)
                    if n > nrecords:
                        nrecords = n
            records = arange(nrecords)
        r0 = self._qnd_state.goto
        for r in records:
            self._qnd_state.goto = r
            yield r
        self._qnd_state.goto = r0

    def auto(self, recurse):
        """Set the auto-read mode for this QGroup.

        In auto-read mode, getting an item returns its value, rather than a
        QLeaf.  If the item is a QGroup or QList, that is returned if
        the `recurse` value is 1, whereas if `recurse` is 2, the QGroup
        or QList variables will be read recursively.  Setting `recurse` to
        0 turns off auto-read mode entirely.

        Note that you can temporarily set auto mode using a with clause.

        """
        self._qnd_state.auto = int(recurse)

    def push(self):
        """Push current recording, goto, and auto mode onto state stack."""
        self._qnd_state.push()
        return self

    def drop(self, nlevels=None, close=False):
        """Restore previous recording, goto, and auto mode settings.

        Default ``drop()`` drops one pushed state, ``drop(n)`` drops n,
        ``drop('all')`` drops all pushed states.  By default, `drop` is
        a no-op if no pushed states to drop, ``drop(close=1)`` closes
        the file if no pushed states to drop, which is called implicitly
        on exit from a with suite.

        """
        if nlevels is None:
            nlevels = 1
        elif nlevels == "all":
            nlevels = len(self._qnd_state) - 3
        while nlevels >= 0:
            if self._qnd_state.drop() and close:
                self.close()
            nlevels -= 1

    def close(self):
        """Close associated file."""
        this = self._qnd_group
        if this is not None:
            for nm in ["_qnd_group", "_qnd_state", "_qnd_cache"]:
                object.__setattr__(self, nm, None)
            this.close()

    def flush(self):
        """Flush associated file."""
        this = self._qnd_group
        if this is not None:
            this.flush()

    def root(self):
        """Return root QGroup for this item."""
        qgroup = self._qnd_group
        root = qgroup.root()
        if root is qgroup:
            return self
        state = QState(self._qnd_state)  # copy
        return QGroup(root, state)

    def attrs(self):
        """Return attribute tree for variables in this group."""
        return QAttributes(self._qnd_group)

    def get(self, key, default=None):
        """like dict.get method"""
        try:
            return self[key]
        except KeyError:
            return default

    def items(self, auto=None):
        """like dict.items method (iteritems in python2)"""
        if auto == self._qnd_state.auto:
            auto = None
        for name in self._qnd_group:
            if auto is None:
                value = self[name]
            else:
                with self.push():
                    self.auto(auto)
                    value = self[name]
            yield name, value

    def __repr__(self):
        this = self._qnd_group
        if this is not None:
            return "<QGroup with {} items>".format(len(this))
        else:
            return "<closed QGroup>"

    def __len__(self):
        return len(self._qnd_group)

    def __contains__(self, name):
        return self._qnd_group.lookup(name) is not None

    def __iter__(self):
        return iter(self._qnd_group)

    keys = __iter__

    __enter__ = push

    def __exit__(self, etype, evalue, etrace):
        self.drop(close=1)

    def __call__(self, auto=None, *args):
        # Make qg() shorthand for qg[()], returning whole group.
        if auto == self._qnd_state.auto:
            auto = None
        if auto is None:
            value = self[args]
        else:
            with self.push():
                self.auto(auto)
                value = self[args]
        return value

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if not key:  # qg[()] retrieves entire group
            key = (list(self._qnd_group),)
        name, args = key[0], key[1:]
        if isinstance(name, basestring):
            if "/" in name:
                if name.startswith("/"):
                    return self.root()[(name[1:],) + args]
                name = name.split("/")
                name, args = name[0], tuple(name[1:]) + args
        else:
            # qg[["name1", "name2", ...], slice0, ...]
            # returns [qg.name1[slice0, ...], qg.name2[slice0, ...], ...]
            items = []
            for key in name:
                if not isinstance(key, basestring):
                    # Prevent recursive name lists inside name lists.
                    raise KeyError("expecting item name or list of item names")
                items.append((key, self[(key,)+args]))
            return ADict(items)
        item = self._qnd_group.lookup(name)
        if item is None:
            raise KeyError("no such item in QGroup as {}".format(name))
        state = self._qnd_state
        auto, recording = state.auto, state.recording
        record = state.goto
        if item.islist():
            item = QList(item, auto)
            if record is None:
                if not args and auto <= 1:
                    return item
            else:
                args = (record,) + args
            return item[args] if args else item[:]
        if item.isleaf():
            return _reader(item, args) if args or auto else QLeaf(item)
        # Item must be a group, set up inherited part of state.
        # Note that goto record was not used, so subgroup inherits it.
        cls = item.lookup("__class__")
        if cls is not None and auto:
            return _load_object(item, cls)
        item = QGroup(item, auto=auto, recording=recording, goto=record)
        return item() if auto > 1 else item

    def __setitem__(self, key, value):
        name, args = (key[0], key[1:]) if isinstance(key, tuple) else (key, ())
        if not isinstance(name, basestring):
            name = "/".join(name)
        if "/" in name:
            if name.startswith("/"):
                item = self.root()
                name = name[1:]
            else:
                path, name = name.rsplit("/", 1)
                with self.push():
                    self.auto = 0
                    item = self[path]
            item[(name,) + args] = value
            return
        dtype, shape, value = _categorize(value)
        state = self._qnd_state
        recording, record = state.recording, state.goto
        this = self._qnd_group
        item = this.lookup(name)
        if item is None:
            # Declare item now.
            if args:
                raise KeyError("partial write during declaration of {}"
                               "".format(name))
            if recording:
                # numpy (1.16.4) misfeature dtype('f8') tests == None
                # (other dtypes are != None as expected), so cannot
                # ask if (list, dict, object, None) contains dtype
                if recording != 1 or (dtype is None or
                                      dtype in (list, dict, object)):
                    # Declare an anonymous-group-style list.
                    item = this.declare(name, list, None)
                else:  # Declare item with UNLIMITED dimension.
                    item = this.declare(name, dtype, shape, 1)
                # item now an empty list
            elif dtype == dict:
                item = this.declare(name, dict, None)
                if value:
                    QGroup(item).update(value)
                return
            elif dtype == list:
                item = this.declare(name, list, None)
                if value:
                    QList(item).extend(value)
                return
            elif dtype == object:
                item = this.declare(name, dict, None)
                _dump_object(item, value)
                return
            else:
                item = this.declare(name, dtype, shape)
            if value is None:
                return
        while item.islist():
            if recording:
                if args:
                    raise KeyError("partial write while recording {}"
                                   "".format(name))
                record = len(item)  # index of next record
                item = item.declare(dict if dtype == object else dtype,
                                    shape)  # declare the next item
                state.goto = record
                if dtype is None:
                    return
                if dtype == list:
                    if value:
                        QList(item).extend(value)
                    return
                if dtype == dict:
                    if value:
                        QGroup(item).update(value)
                    return
                if dtype == object:
                    _dump_object(item, value)
                    return
                break
            if record is None:
                if args:
                    record, args = args[0], args[1:]
                elif dtype == list and not value:
                    # qg.lst = list is no-op for existing QList
                    return
                else:
                    raise ValueError("cannot set existing QList {}, use "
                                     "append, goto, or recording".format(name))
            item = item.index(record)
            if item is None:
                raise KeyError("no such item in QList as {}".format(record))
            record = None  # makes no sense to use record recursively
            recording = 0
        if item.isgroup():
            if args:
                QGroup(item)[args] = value
                return
            if dtype == dict and not value:
                # qg.grp = {} is no-op for existing QGroup
                return
            raise ValueError("cannot set existing QGroup {}, use update"
                             "".format(name))
        # item is a leaf (neither a list nor a group)
        if dtype in (dict, list, object):
            raise TypeError("type mismatch in QLeaf {}".format(name))
        elif item.query()[0] is None:
            # None QLeaf objects need not support write() method.
            if dtype is None and not args:
                return
            raise TypeError("QLeaf {} declared as None".format(name))
        item.write(value, args)


def _monotonize(values):
    # This function ensures values are monotonically increasing,
    # searching backwards for decreasing sequences.
    decreasing = values[-1] < values[0]
    if decreasing:
        values = -values
    mask = values == values
    vnext = values[-1]
    for i in range(-2, -values.size-1, -1):
        v = values[i]
        if v >= vnext:
            mask[i] = False
        else:
            vnext = v
    records, values = where(mask)[0], values[mask]
    if decreasing:
        # Reverse both records and values so latter is strictly increasing.
        # In this way, values can always be used as x in the interp function.
        records = records[::-1]
        values = -values[::-1]
    return records, values


def _categorize(value, attrib=False):
    # This function defines the various sorts of values QnD recognizes:
    # 1. None           dtype = shape = value = None
    # 2. list [, seq]   dtype = list, shape = None, value = [] or seq
    # 3. {...}          dtype = dict, shape = None, value = {...}
    # 4. type|dtype [, shape]  dtype = dtype(type), shape, value = None
    # 5. array_like     value.dtype, value.shape, value = asanyarray(...)
    # 6. object or dtype('O')  dtype = object, shape = None, value
    if value is None:
        dtype = shape = None
    elif isinstance(value, (type, _dtype)):
        if value == list:
            dtype, shape, value = list, None, []
        else:
            dtype, shape, value = _dtype(value), (), None
    elif isinstance(value, dict):
        if all(isinstance(key, basestring) for key in value):
            dtype = dict
        else:
            dtype = object
        shape = None
    elif (isinstance(value, tuple) and len(value) == 2+bool(attrib) and
          isinstance(value[0], (type, _dtype))):
        dtype = value[0]
        if dtype is not None and dtype not in (list, dict, object):
            dtype = _dtype(dtype)  # no-op if already a dtype
        if not attrib:
            if dtype == list:
                value, shape = value[1], None
            else:
                value, shape = None, tuple(value[1])
        else:
            shape, value = value[1:]
    else:
        # The array(a) constructor used to accept essentially any argument a.
        # At numpy 1.19 it began issues a VisibleDeprecationWarning when a
        # was a list whose items were of differing lengths (or shapes).
        # Prior to that, it simply produced an ndarray of dtype object whose
        # items were the python entities in the original list.  This is the
        # behavior we want in QnD, so we do not want to print a warning.
        # Moreover, when the feature is eventually removed, this case will
        # throw a (currently unknown) exception, which we need to avoid.
        # Passing the dtype=object keyword to the array() constructor
        # produces the pre-1.19 behavior (as far as I can tell), but of
        # course we cannot do that here.
        # The following code must work in three cases: (1) pre-1.19 numpy,
        # (2) numpy 1.19-1.21 (at least) which print unwanted warnings without
        # special treatment, and (3) future numpy which throws an error
        # without the dtype=object keyword.  Since QnD must always run in
        # all three cases, there is no way to remove the protection against
        # the deprecation wawrning, even when numpy move past it.
        # As a further complication, VisibleDeprecationWaring has been removed
        # in numpy 2.0.
        with catch_warnings():
            # Make case 2 (numpy 1.19) behave like case 3 (future numpy)
            simplefilter("error", VisibleDeprecationWarning)
            try:
                v = asanyarray(value)
            except Exception:
                # As far as I can tell, the original numpy array() constructor
                # would accept any argument whatsoever, returning either a
                # scalar or 1D array of type object if its argument could not
                # be interpreted.  Therefore I believe only a ragged array
                # argument reaches this point, and we can return to the
                # original behavior by specifying dtype explicitly.
                # Nevertheless, we protect against a possible exception.
                simplefilter("ignore", VisibleDeprecationWarning)
                try:
                    v = asanyarray(value, dtype=object)
                except Exception:
                    return object, None, value
        dtype, shape = v.dtype, v.shape
        if dtype.kind == "O":
            if not shape:
                dtype, shape = object, None
            else:
                # Note that this does not work as expected when the contents
                # of the list were themselves lists (not ndarrays) of numbers
                # of varying lengths, since the asanyarray function will not
                # convert those inner lists to ndarrays.  Hence v.tolist() is
                # really the same as the original value here.
                # A QnD user must ensure that the inner lists are ndarrays if
                # that is what they intended.
                dtype, shape, value = list, None, v.tolist()
        else:
            value = v
    if isinstance(dtype, _dtype):
        kind = dtype.kind
        if kind == "U":
            if value is not None:
                value = npencode(value, "utf8")  # convert to 'S'
                dtype = value.dtype
        elif kind == "O":
            raise ValueError("numpy dtype.kind 'O' not supported")
    return dtype, shape, value


def _reader(item, args):
    value = item.read(args)
    dtyp = getattr(value, "dtype", None)
    if dtyp is not None:
        kind = dtyp.kind
        if kind == "V":
            if dtyp.names:
                # The recarray has some significant misfeatures.  The worst
                # is that it will not print (repr or str) if it is aligned,
                # or simply if the itemsize does not match what it expects.
                # value = value.view(recarray)
                pass
        elif kind in "SU":
            value = _totext(value)
    if isinstance(value, ndarray) and not value.shape:
        value = value[()]
    return value


def _totext(value):
    if not PY2 and value.dtype.kind == "S":
        try:
            value = npdecode(value, "utf8")
        except UnicodeDecodeError:
            value = npdecode(value, "latin1")
    if not value.ndim:
        value = value.tolist()  # make scalar text a real string instance
    return value


_dtype = dtype  # to allow access in methods using local name dtype
_builtin_module = str.__class__.__module__


def _dump_object(item, value):
    # item.isgroup() is true, as yet empty, value is an object
    item = QGroup(item)
    if isinstance(value, dict):
        # special case for dict with non-text keys
        item["__class__"] = "dict"
        items = value.iteritems if PY2 else value.items
        for i, (k, v) in enumerate(items()):
            item["_" + str(2*i)] = k
            item["_" + str(2*i+1)] = v
    else:
        cls = value.__class__
        cname, module = cls.__name__, cls.__module__
        if module is not None and module != _builtin_module:
            cname = ".".join((module, cname))
        item["__class__"] = cname
        # We do not support the python2-only __getinitargs__.
        mydict = getattr(value, "__dict__", None)
        getter = getattr(value, "__getstate__", None)
        newargs, newkwargs = (), {}
        if hasattr(value, "__getnewargs_ex__"):
            newargs, newkwargs = value.__getnewargs_ex__()
        elif hasattr(value, "__getnewargs__"):
            newargs = value.__getnewargs__()
        if newargs:
            item["__getnewargs__"] = list, newargs
        if newkwargs:
            item["__getnewargs_ex__"] = newkwargs
        state = getter() if getter else None
        if state:
            slots = None
            if hasattr(value, "__slots__"):
                if mydict is not None:
                    mydict, slots = state
                else:
                    slots = state
            elif mydict is not None:
                mydict = state
            if slots is not None:
                item["__slots__"] = slots
        else:
            # No high level pickle interface, check for simple low level:
            reduced = value.__reduce__()
            if reduced[0] == cls and not any(v is not None
                                             for v in reduced[2:]):
                item["__reduce__"] = list, reduced[1]
                mydict = None
        if mydict is not None:
            setter = hasattr(value, "__setstate__")
            item["__setstate__" if setter else "__dict__"] = mydict


def _load_object(qgroup, cls):
    # If you fail here, you can still read the group with ADict(qgroup)
    # which avoids this special treatment.
    cls = cls.read()  # assume QLeaf yields a text string
    dtyp = getattr(cls, "dtype", None)
    if dtyp is None or not dtyp.kind in "SU" or cls.ndim:
        raise TypeError("Expecting __class__ member of QGroup to be text.")
    cls = _totext(cls)
    qgroup = QGroup(qgroup, auto=2)
    if cls == "dict":
        obj = {}
        names = list(name for name in qgroup if name != "__class__")
        if len(names) & 1:
            names[0] = ""  # die in first pass
        key = None
        for i, n in enumerate(sorted(names)):
            if "_{}".format(i) != n:
                raise TypeError("QGroup with __class__ dict error")
            value = qgroup[n]
            if i & 1:
                obj[key] = value
            else:
                key = value
        return obj
    elif cls == "ellipsis":
        return Ellipsis
    cls = cls.rsplit(".", 1)
    try:
        module = (import_module(cls[0]) if len(cls) > 1 else
                  sys.modules[_builtin_module])
        cls = getattr(module, cls[-1])
    except (ImportError, AttributeError):
        # If the named class does not exist or does not have
        # the specified class, just return an ADict.
        return ADict(qgroup)
    reduced = qgroup.get("__reduce__")
    if reduced is not None:
        return cls(*reduced)
    newargs = qgroup.get("__getnewargs__", ())
    newkwargs = qgroup.get("__getnewargs_ex__", {})
    obj = cls.__new__(cls, *newargs, **newkwargs)
    if "__setstate__" in qgroup:
        obj.__setstate__(qgroup["__setstate__"])
    else:
        if "__dict__" in qgroup:
            obj.__dict__.update(qgroup["__dict__"])
        if "__slots__" in qgroup:
            obj.__slots__.update(qgroup["__slots__"])
    return obj


class QState(list):
    """State information for a QGroup."""
    __slots__ = ()

    def __init__(self, recording=0, goto=None, auto=0):
        if hasattr(recording, "__iter__"):
            seq = tuple(recording)[:3]
        else:
            if goto is not None:
                goto = int(goto)
            recording, auto = int(recording), int(auto)
            seq = recording, goto, auto
        super(QState, self).__init__(seq)

    @property
    def recording(self):
        return self[0]

    @recording.setter
    def recording(self, value):
        self[0] = int(value)

    @property
    def goto(self):
        return self[1]

    @goto.setter
    def goto(self, value):
        self[1] = None if value is None else int(value)

    @property
    def auto(self):
        return self[2]

    @auto.setter
    def auto(self, value):
        self[2] = int(value)

    def push(self):
        state = self[:3]
        self.append(state)

    def drop(self):
        if len(self) < 4:
            return 1
        self[:3] = super(QState, self).pop()
        return 0


class QList(object):
    """List of subgroups, lists, and ndarrays.

    You reference QList elements by index or slice, like ordinary list
    elements, including the python convention for negative index values.
    To read the entire list, call it like a function, ``ql()``, which is
    equivalent to ``ql[:]``.  A QList has __iter__, append, and extend::

       for element in ql: do_something
       ql.append(value)
       ql.extend(iterable)

    In general, the elements of a QList are unrelated to one another;
    it's like an anonymous QGroup.  However, a common use case is to
    represent a so-called UNLIMITED dimension in netCDF or HDF5.  In
    this case, every element will have the same dtype and shape.  The
    `islist` method returns 1 for this special restricted case, while
    it returns 2 for an unrestricted QList.  Whether this makes any
    difference depends on the underlying file format.  The QGroup
    `recording` and `goto` methods allow you to access QList items in
    the group transparently, as if they were individual elements at
    a current record or index.

    Attributes
    ----------
    isgroup
    isleaf
       Always 0.
    islist
       This is 1 if this QList is a record array declared in recording
       mode 1, and 2 if it was declared in any other way (including as a
       record array in recording mode 2).
    dtype
       Always ``list``, the builtin python type.
    shape
    ndim
    size
    sshape
       Always None.

    """
    __slots__ = "_qnd_list", "_qnd_auto"
    isgroup = isleaf = 0
    dtype, shape, ndim, size, sshape = list, None, None, None, None

    def __init__(self, item=None, auto=0):
        object.__setattr__(self, "_qnd_list", item)
        self.auto(auto)

    def auto(self, recurse):
        """Set auto read mode, analogous to QGroup.auto method."""
        object.__setattr__(self, "_qnd_auto", int(recurse))

    def root(self):
        """Return root QGroup for this item."""
        return QGroup(self._qnd_list.root(), QState(auto=self._qnd_auto))

    @property
    def islist(self):
        return self._qnd_list.islist()

    def extend(self, iterable):
        """append multiple new elements to this QList"""
        for value in iterable:
            self.append(value)

    def append(self, value):
        """append a new element to this QList"""
        dtype, shape, value = _categorize(value)
        item = self._qnd_list.declare(dtype, shape)
        if dtype is None:
            return
        if dtype == list:
            if value:
                QList(item).extend(value)
            return
        if dtype == dict:
            if value:
                QGroup(item).update(value)
            return
        if dtype == object:
            _dump_object(item, value)
            return
        if value is not None:
            item.write(value, ())
            # Being unable to do partial write on declaration is consistent
            # with behavior of QGroup __setitem__.  The way to get it is to
            # make a declaration with value = (type, shape) instead of an
            # actual value in both cases.

    def __repr__(self):
        return "<QList with {} items>".format(len(self))

    def __len__(self):
        return len(self._qnd_list)

    def __iter__(self):
        auto = self._qnd_auto
        recurse = auto > 1
        for item in self._qnd_list:
            if item.isgroup():
                cls = item.lookup("__class__") if auto else None
                if cls is None:
                    item, readit = QGroup(item), recurse
                else:
                    item, readit = _load_object(item, cls), 0
            elif item.islist():
                item, readit = QList(item), recurse
            else:
                item, readit = QLeaf(item), auto
            yield item() if readit else item

    def __call__(self):
        return self[:]

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        index, args = key[0], key[1:]
        this = self._qnd_list
        if isinstance(index, slice):
            index = range(*index.indices(len(this)))
        if hasattr(index, "__iter__"):
            return [self[(i,) + args] for i in index]
        item = this.index(index)
        if item is None:
            raise IndexError("QList index {} out of range".format(index))
        auto = self._qnd_auto
        if item.islist():
            item = QList(item, auto)
            if args:
                return item[args]
            return item[:] if auto > 1 else item
        if item.isleaf():
            return _reader(item, args) if args or auto else QLeaf(item)
        # Item must be a group, set up inherited part of state.
        # Note that goto record was not used, so subgroup inherits it.
        cls = item.lookup("__class__")
        if cls is not None and auto:
            return _load_object(item, cls)
        item = QGroup(item, auto=auto)
        return item() if auto > 1 else item

    def __setitem__(self, key, value):
        if not isinstance(key, tuple):
            key = (key,)
        index, args = key[0], key[1:]
        if isinstance(index, slice) or hasattr(index, "__iter__"):
            raise TypeError("QList does not support multi-element setitem")
        dtype, shape, value = _categorize(value)
        item = self._qnd_list.index(index)
        if item is None:
            raise IndexError("QList index {} out of range".format(index))
        if item.islist() or item.isgroup():
            idtype = list if item.islist() else dict
            if idtype == dtype and not value:
                return
            raise TypeError("cannot set existing QGroup or QList")
        # item is a QLeaf
        if item.query()[0] is None:
            if dtype is None and not args:
                return
            raise TypeError("QLeaf {} declared as None".format(index))
        # Work around numpy (1.16.4) misfeature dtype('f8') tests == None:
        if dtype is None or dtype in (list, dict, object):
            raise TypeError("type mismatch setting QLeaf {}".format(index))
        item.write(value, args)


class QLeaf(object):
    """An ndarray or None stored in a file.

    You can read the data by calling the leaf instance ``ql()``, or by
    indexing it ``ql[:]``, which also provides a means for partial reads.
    A QLeaf has `dtype`, `shape`, `ndim`, and `size` properties with the
    same meanings as an ndarray (except None has all these properties
    equal None).  Additionally, the `sshape` property may return a symbolic
    shape with optional strings in the tuple representing dimension names.

    You can write data by calling ``ql(value)``, or by setting a slice,
    which provides a means for partial writes.

    Attributes
    ----------
    isgroup
    islist
       Always 0.
    isleaf
       Always 1.
    dtype
       The numpy dtype of this ndarray, or None if this leaf is None.
       This is the dtype in memory, not necessarily as stored.
    shape
    ndim
    size
       The numpy ndarray properties, or None if this leaf is None.
    sshape
       A symbolic shape tuple, like shape except dimension lengths may be
       type str instead of int.

    """
    __slots__ = "_qnd_leaf",
    isgroup = islist = 0
    isleaf = 1

    def __init__(self, item):
        object.__setattr__(self, "_qnd_leaf", item)

    def root(self):
        """Return root QGroup for this item."""
        return QGroup(self._qnd_leaf.root())

    def __call__(self, value=_NOT_PRESENT_):
        if value is _NOT_PRESENT_:
            return self[()]
        else:
            self[()] = value

    def __getitem__(self, key):
        return _reader(self._qnd_leaf,
                       key if isinstance(key, tuple) else (key,))

    def __setitem__(self, key, value):
        self._qnd_leaf.write(value, key if isinstance(key, tuple) else (key,))

    @property
    def dtype(self):
        return self._qnd_leaf.query()[0]

    @property
    def shape(self):
        return self._qnd_leaf.query()[1]

    @property
    def ndim(self):
        return len(self._qnd_leaf.query()[1])

    @property
    def size(self):
        shape = self._qnd_leaf.query()[1]
        return prod(shape) if shape else 1

    @property
    def sshape(self):
        _, s, ss = self._qnd_leaf.query()
        return ss if ss else s


class QAttributes(ItemsAreAttrs):
    """Attributes for a QGroup and its members.

    Usage::

       qa = qgroup.attrs()
       qa0 = qa.vname  # for variables in this group, or qa['vname']
       qa1 = qa._  # or qa[''] for attributes of this group
       value = qa0.aname  # or qa0['aname'], None if no such attribute
       qa0.aname = value  # or qa0['aname'] = value
       qa0.aname = dtype, shape, value
       if 'aname' in qa0: do_something
       for aname in qa0: do_something
       for aname, value in qa0.items(): do_something

    """
    __slots__ = "_qnd_parent", "_qnd_vname", "__weakref__"

    def __init__(self, parent, vname=None):
        if not isinstance(parent, ProxyTypes):
            parent = proxy(parent)
        object.__setattr__(self, "_qnd_parent", parent)
        object.__setattr__(self, "_qnd_vname", vname)

    def __repr__(self):
        vname = self._qnd_vname
        if vname is None:
            return "<QAttributes accessor for QGroup items>"
        elif not vname:
            return "<QAttributes for whole QGroup>"
        return "<QAttributes for item {}>".format(vname)

    def get(self, key, default=None):
        parent, vname = self._qnd_parent, self._qnd_vname
        if vname is None:
            # Get group attribute, even though that is inconsistent...
            # Should we implement matching set() or just let it go?
            vname = ""
        else:
            parent = parent._qnd_parent
        return parent.attget(vname).get(key, default)

    def keys(self):
        group, vname = self._qnd_group_vname()
        return iter(group.attget(vname))

    def items(self):
        group, vname = self._qnd_group_vname()
        return group.attget(vname).items()

    def _qnd_group_vname(self):
        parent, vname = self._qnd_parent, self._qnd_vname
        if vname is None:
            raise TypeError("need to specify QGroup item name")
        return parent._qnd_parent, vname

    def __getattr__(self, name):
        vname = self._qnd_vname
        if vname is None or name not in self._qnd_builtins_:
            return super(QAttributes, self).__getattr__(name)
        # Handle builtin pseudo-attributes here; they do not show up
        # in the actual attribute dict referenced by [key].
        # Can use dtype_, shape_, etc. attributes if real attributes
        # have these names.
        item = self._qnd_parent._qnd_parent.lookup(vname)
        if item.isgroup():
            return dict if name == "dtype" else None
        if item.islist():
            return list if name == "dtype" else None
        dsss = item.query()
        if dsss[0] is None:
            return None
        if name == "ndim":
            return len(dsss[1])
        if name == "size":
            return prod(dsss[1])
        return dsss[self._qnd_builtins_.index(name)]

    _qnd_builtins_ = ["dtype", "shape", "sshape", "size", "ndim"]

    def __getitem__(self, key):
        parent, vname = self._qnd_parent, self._qnd_vname
        if vname is None:
            # key is vname
            item = parent.lookup(key) if key else True
            if item is None:
                raise KeyError("no such item in QGroup as {}".format(key))
            return QAttributes(self, key)
        return parent._qnd_parent.attget(vname).get(key)

    def __setitem__(self, key, value):
        group, vname = self._qnd_group_vname()
        # Note that value can be (dtype, shape, value) to be explicit.
        dtype, shape, value = _categorize(value, 1)
        if dtype in (list, dict, object):
            raise TypeError("an attribute cannot be a dict or list")
        group.attset(vname, key, dtype, shape, value)

    def __iter__(self):
        group, vname = self._qnd_group_vname()
        return iter(group.attget(vname))

    def __contains__(self, name):
        group, vname = self._qnd_group_vname()
        return name in group.attget(vname)

    def __len__(self):
        group, vname = self._qnd_group_vname()
        return len(group.attget(vname))


class QnDList(object):
    """Implmentation of a low level QList type using QGroup.

    A backend which has no direct support for QList objects can use
    this to produce a pseudo-list, which is a group with member names
    _ (None or a single signed or unsigned byte, value never read) and
    names _0, _1, _2, etc.

    This implementation will handle both UNLIMITED index-style lists
    made with recording = 1 (that is group.declare with unlim flag)
    and general lists.  If UNLIMITED dimensions are supported, pass the
    QnDLeaf to this constructor::

       item = QnDList(QnDLeaf)  # if at least one record exists
       item = QnDList(QnDLeaf, 1)  # if no records yet exist

    Use the fromgroup constructor to check if a QnDGroup is a pseudo-list::

       item = QnDList.fromgroup(QnDGroup)

    """
    __slots__ = "_qnd_parent", "_qnd_current",

    def __init__(self, parent, empty=None):
        self._qnd_parent = parent
        current = empty
        if empty is not None:
            if parent.isgroup():
                parent.declare("_", None, ())
            elif not isinstance(parent, QnDList):
                current = -1
        self._qnd_current = current

    @staticmethod
    def fromgroup(parent):
        item = parent.lookup("_")
        if item is not None:
            if all(_us_digits.match(name) for name in parent):
                return QnDList(parent)  # parent is a pseudo-list
        return parent

    def parent(self):
        parent = self._qnd_parent
        return parent._qnd_parent if isinstance(parent, QnDList) else parent

    @staticmethod
    def isgroup():
        return 0

    def isleaf(self):
        return int(isinstance(self._qnd_parent, QnDList))

    def islist(self):
        if self._qnd_parent.isgroup():
            return 2
        return int(not isinstance(self._qnd_parent, QnDList))

    def root(self):
        return self._qnd_parent.root()

    # len, iter, index, declare are list methods, assume isleaf() false

    def __len__(self):
        if self._qnd_parent.isgroup():
            return len(self._qnd_parent) - 1  # subtract _ item
        if self._qnd_current is not None and self._qnd_current < 0:
            return 0  # leaf.query() probably returns 1
        return self._qnd_parent.query()[1][0]

    def __iter__(self):
        parent = self._qnd_parent
        if parent.isgroup():
            for i in range(len(self)):
                yield parent.lookup("_" + str(i))
        else:
            for i in range(len(self)):
                yield QnDList(self, i)

    def index(self, ndx):
        nrecs = max(len(self), 1)
        if ndx < 0:
            ndx = ndx + nrecs
        if ndx < 0 or ndx >= nrecs:
            return None  # out of range, let caller raise any exception
        parent = self._qnd_parent
        if parent.isgroup():
            return parent.lookup("_" + str(ndx))
        return QnDList(self, ndx)

    def declare(self, dtype, shape):
        parent = self._qnd_parent
        nrecs = len(self)
        if parent.isgroup():
            return parent.declare("_" + str(nrecs), dtype, shape)
        return QnDList(self, nrecs)

    # query, read, write are leaf methods, assume isleaf() true

    def query(self):
        qndlist = self._qnd_parent
        dtype, shape, sshape = qndlist._qnd_parent.query()
        shape = shape[1:]
        if sshape:
            sshape = sshape[1:]
        return dtype, shape, sshape

    def read(self, args=()):
        current = self._qnd_current
        qndlist = self._qnd_parent
        check = qndlist._qnd_current
        if check is not None and check < 0:
            raise TypeError("attempt to read from empty UNLIMITED array")
        return qndlist._qnd_parent.read((current,) + args)

    def write(self, value, args=()):
        qndlist = self._qnd_parent
        qndlist._qnd_parent.write(value, (self._qnd_current,) + args)
        # Turn off special empty list state (if on):
        qndlist._qnd_current = None
