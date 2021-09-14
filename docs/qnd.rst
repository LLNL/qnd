User Interface for Binary Files
===============================

In terms of the scipy environment, qnd addresses the storage and
retrieval of numpy ndarrays, excepting arrays with the general python
object data type (dtype.kind 'O').  In scipy programs, these arrays
may be attributes of objects surrounded by methods and other non-data,
but we assume that the programmer provides some means of initializing
all such objects given just the (numerical) ndarrays at their heart.

The python language provides two kinds of collections which qnd
directly supports: the dict (a collection of named objects) and the
list (a sequence of heterogeneous objects).  We presume that a
programmer will provide a way to map a nested tree (no loops) of dicts
and lists with ndarray leaves (and dict keys which are strings) to and
from whatever objects their program requires.  The dict and list are
precisely the collections provided by the simple and popular JSON data
interchange format.  Thus, in order to use the qnd storage interface,
we are essentially asking the programmer support a portable
organization of the program data.

Note the contrast to the goal of the python pickle module; we
acknowledge that some extra design and maintenance work may be
required to support such a mapping.  Often the additional effort pays
off in a simplified overall design.

Basic Usage
-----------

The first step is to obtain a file handle, say `f`, by opening a file.
The open function belongs to the particular backend, called ``openXXX``,
and defined in a backend module ``XXXf``.  For example::

  from qnd.h5f import openh5
  f = openh5('filename.h5', 'r+')

The qnd mode choices are the same for all backends (copied from the
excellent h5py module mode sematics):

  * 'r' opens the file read-only
  * 'w' creates a new file, clobbering any existing file
  * 'a' opens the file read-write, creating it if it did not exist
  * 'r+' opens the file read-write, raising an error if it did not exist
  * 'w-' creates a new file, raising an error if it exists beforehand

These mode flags are not semantically identical to the python ``open``
function: The 'w-' is not recognized by ``open`` at all, and the
``open`` 'a' guarantees that any existing file bytes will not be
modified, while the qnd 'a' merely means read-write (but has the same
semantics as 'a' in terms of file existence and creation).
Furthermore, qnd files are always readable, even if opened in one of
the 'w' modes.

Python syntax has two operators for extracting named members from a
compound object: The dot operator extracts an `attribute` from an
object, and square brackets extract an `item` from a dict object (or
other mapping).  Qnd file handles support both.  The dot syntax
``f.var`` is best when you know the name 'var' at the time you write
the expression, while the square bracket syntax ``f[expression]`` is
best when the name is the result of an expression or value of a
variable.

The dot operator is overloaded, since it is also used for method
attributes like ``f.close()``.  Thus the qnd file handle `f` really
behaves like a dict for the most part, with its support for the dot
syntax mere sugar to improve code legibility and, at least as
importantly, ease of typing in interactive usage.  If there were a
variable named 'close' in `f` (who knows where `f` came from), you
could always access it as ``f['close']``.  However, qnd provides a
quick and dirty option for using the dot operator even in these cases:
it will remove a single trailing underscore, so that ``f.close_``
refers to the variable ``'close'``, not ``'close_'``.  (``f.close__``
would refer to ``'close_'``.)  This idiom is suggested by the PEP8
python style guide, and you would also need it to escape python
keywords, like ``f.yield_`` to refer to ``f['yield']``.

The bottom line is, you use a qnd file handle `f` as if it were a
python dict, but you are also free to treat items in `f` as if they
were attributes of this dict object::

  x = f.x           # read variable "x" from f, same as f['x']
  f.x = expression  # declare and write "x" to f
  f.update(x=expr1, y=expr2, ...)  # declare and write several variables
  # update also accepts non-keyword dicts and lists of (name, value)
  x = f.get('x', xdefault)  # same as get from dict
  varnames = list(f)        # preferred over f.keys(), as for any dict
  nvars = len(f)
  if name in f: do_something
  for name in f: do_something
  for name, value in f.items(): do_something

In addition to the dict-like `update`, `get`, `keys`, and `items` methods,
qnd files also have a number of non-dict methods and behaviors::

  f.close()
  f.flush()  # like close then reopen
  with openh5('myfile.h5', 'a') as f:
      write_something(f)  # closing f upon exit from with suite
  f.auto(0)        # turn off (or on) auto-read mode
  f.recording(1)   # turn on (or off) recording mode
  f.goto(time=t)   # set to previously recorded record
  with f.push():
      do_something(f)  # temporarily change auto, recording, goto state

The `recording` and `goto` modes are the subject of the next section;
we conclude this section by discussing `auto` mode.  You may have
noticed that ``f.x`` or ``f['x']`` immediately read the variable from
the file, giving you no opportunity to query its data type or shape,
which you might well want to do without incurring the overhead of the
actual read, especially if you know it is a very large array.  We can
get the names of all stored variables with ``list(f)``, but how do we
find out what each one looks like without reading it?

The answer is that a qnd file `f` can be placed into a mode in which
variable references do not trigger an automatic read operation, by
invoking ``f.auto(0)``.  You can also request this mode using the
``auto=0`` keyword when you open the file.  (The default is
``auto=1``.)  With autoread mode off, getting an item returns a qnd
leaf object, which is like a mini-file handle you can use to query,
read, or write only that specific variable.  It has properties similar
to an ndarray::

  f.auto(0)
  xhandle = f.x  # or f['x']
  dtype, shape = xhandle.dtype, xhandle.shape  # also size and ndim
  xhandle = f(0, 'x')  # return handle to x independent of auto mode
  x = xhandle[:]  # read x if x is not scalar
  x = xhandle[()]  # read x no matter what
  xhandle[()] = expression  # write x no matter what
  x = xhandle()  # shorthand for xhandle[()]
  xhandle(expression)  # shorthand for xhandle[()] = expression
  xpart = xhandle[index_expressions]  # read part of x
  xhandle[index_expressions] = xpart  # write part of x

Notice that `xhandle` inherits the obscure indexing behavior of
ndarray scalars, for which ``x[:]`` raises an error.  However,
`xhandle` provides a non-ndarray operation to compensate -- calling a
qnd handle as a function always reads the whole thing, whether or not
it has any dimensions.

Although the qnd leaf handles can be used for partial read and write
operations, if that is all you want to do, you can simply combine the
partial index expressions into a single square bracket::

  xpart = f['x', index_expressions]
  f['x', index_expressions] = xpart

These work no matter how the autoread mode is set, but there is no
equivalent using the dot syntax: Although ``f.x[index_expressions]``
produces the same final result, it reads all of `x` before applying
`index_expressions` to the resulting large ndarray.

(Note that qnd only reads or writes the largest contiguous block of
leading indices specified by `index_expressions`; it only reduces the
intermediate memory footprint when the leading indices are scalar or
small slices of `x`.)

Finally, sometimes you need to declare a variable without writing it.
To do this in qnd, make its value a dtype or a (dtype, shape) tuple::

  f.x = float  # declare x to be a scalar dtype(float), that is f8
  f.y = yy.dtype, yy.shape  # declare y with type and shape of yy
  f.z = bool, yy.shape  # declare z to be boolean with same shape as yy

Such a declaration reserves space for the array in the file, but it is
your responsibility to fill it with sensible values with one later
write or several partial writes.

Recording History
-----------------

Setting an item with ``f.x = value`` or ``f['x'] = value`` both
declares the variable and writes its value.  If you later write it a
second time with ``f.x = value2``, by default this overwrites the
orginal value you wrote.  Sometimes, however, you need to record the
history of a variable which is changing as a simulation progresses.
The idea behind recording mode is to make the second assignment store
the new `value2` in addition to the original `value`, so by repeatedly
assigning values to `x` you can store as many versions of its changing
values as you like.

The HDF5, netCDF, and PDB file formats all support this capability by
allowing the leading dimension of a variable to be "unlimited".  But
in qnd, you can suppress this fictitious leading dimension by using
the `recording` mode to write such variables, and the `goto` mode to
read them::

  f = openh5('myfile.h5', 'w')
  f.x = xa  # x is not a record variable.
  f.recording(1)  # Put f in recoding mode; new variables are recorded.
  f.time = t0  # Time is a record variable with t0 for its first record.
  f.y = y0  # y is a record variable with y0 for its first value.
  f.x = xb  # x remains a non-record variable, xb overwrites xa
  f.time = t1  # Write a second record of time with value t1.
  f.y = y1  # Write a second record of y with value y1.
  f.close()

  f = openh5('myfile.h5', 'r')
  # Initially, goto mode is off (None), and reading a record variable...
  times = f.time[:]  # ...returns a list (not array) of all of its records.
  # Use goto to set a "current record" index for all record variables:
  f.goto(0)  # first record
  t0 = f.time
  y0 = f.y
  xb = f.x  # non-record variables ignore current record
  with f.push():  # current record restored on exit from with suite
      f.goto(-1)  # go to last record, record<0 acts like any other index
      yN = f.y
  # You may use any scalar record variable as a keyword to jump to the
  # record nearest the specified value of that variable (assuming it is
  # monotonic):
  f.goto(time=1.2)  # set to record where f.time nearest 1.2
  y12 = f.y
  for record in f.gotoit():  # iterate over all records
      # gotoit() causes implicit f.goto(record) before each pass
      do_something(f)
  f.goto(None)  # Turn off goto mode.
  ylist = f.y  # list of y arrays at every record

The qnd interface, unlike the existing backend file formats, also supports
the case of record variables whose shape changes from one record to the
next.  To use this feature, set the recording mode to 2 instead of to 1::

  f.recording(2)
  f.x = zeros((nx, ny))  # First x record has shape (nx, ny).
  f.x = zeros((nx+5, ny-2))  # Second x record has shape (nx+5, ny-2).
  f.goto(None)
  xlist = f.x  # list of x arrays at every record

This possibility explains why ``f.recordvar`` returns a list of values at
every record, rather than an array with an extra leading dimension (as in
the fiction employed for the existing file formats).

Groups and lists of variables
-----------------------------

The qnd file handle class is `QGroup`; specifically it is the "root
group" of the file.  But a QGroup may contain subgroups, just as a
python dict may contain other dicts.  To define a subgroup, simply
assign a dict instead of an array-like value to an item::

  f.g = {}  # declare an empty subgroup g
  f.g.update(x=expr1, y=expr2)  # all the methods of f work with g
  g = f.g  # g is a QGroup, a subgroup of f
  y = g.y  # or g['y']
  g.auto(0)  # initially g inherits autoread and other modes from f
  root = g.root()  # returns root QGroup, root is f here
  if f is f.root(): task_if_f_is_root_group()
  f['g/x']  # same as f.g.x
  f['/g/x']  # same as f.root().g.x

Although a subgroup initially inherits its autoread, recording, and goto
modes from its parent, thereafter the modes of `g` are independent of the
modes of `f`.  In a `gotoit` loop, the record number in the iterator will
be necessary to explicitly keep subgroups synchronized::

  g = f.g
  for record in f.gotoit():
      g.goto(record)
      do_something(f, g)

Because of the the fact that a `QGroup` looks like a dict, ``dict(f)``
will read every variable in `f`.  By analogy with the qnd leaf
handles, ``f()`` also reads every item in `f` into a dict, with one
twist: Instead of an ordinary dict, ``f()`` results in a dict subclass
called an `ADict`, which permits access to the dict items as
attributes according to the same rules as for a `QGroup`.  If you want
to convert your own `dict` objects into `Adict` objects, you can use
the `redict` function in the ``qnd.adict`` module.  That module also
contains a generic mix-in class `ItemsAreAttrs` which you can use as a
base class for your own mapping classes.  (Although be sure you read
the comment in the `__getattr__` method before you attempt this, as it
can make your code difficult to debug.)

Note that ``f()`` respects the autoread and goto modes.  Thus if
``auto=0``, you nothing will be read from the file and the returned
dict will contain qnd leaf handles (`QLeaf` objects) rather than
variable values.  When ``auto=1``, the dict item corresponding to any
subgroup will be a `QGroup` object.  If you want to recursively read
all subgroups, set ``auto=2``, which causes subgroups to be read
automatically.  (Note that since ``g = f.g`` produces an `ADict` in
that case rather than a `QGroup`, ``auto=2`` can never be inherited.)

In addition to `QGroup` (a dict with str keys) and `QLeaf` (an
ndarray), the qnd interface provides a third item type, `QList`, which
stores a python heterogeneous list.  A `QList` is a way to store a
sequence of objects anonymously, so that you can reference them simply
by a sequence number instead of by a name.  If you find yourself
inventing sequences of names like 'var00', 'var01', var02', and so on,
to store in a `QGroup`, you want to use a `QList` instead::

  f.var = list  # (the builtin list type) declares empty list var
  var = f.var  # the QList object, assuming f.goto(None)
  var.append(value0)  # QList has list-like append and extend methods
  var.append(value1)
  var.extend([value0, value1, ...])
  value1 = var[1]  # second item of var, negative index, slices work
  var[1] = newvalue1  # overwrite value1
  nitems = len(var)
  var.auto(0)  # QList initially inherits its parent's autoread mode

Although `QList` has an autoread mode like a `QGroup`, it does not have
either a recording mode or a goto mode.  In fact, a record variable is
implemented as a `QList`, so the recording and goto modes in the parent
group will influence how the list presents itself::

  f.goto(1)
  value1 = f.var  # In goto mode, f.var means f.var[current_record].

The ability to store aribtrary str-keyed dict and list trees whose
leaves are ndarrays (or None) gives qnd the ability to support pretty
much arbitrary python objects.  In particular, anything which can be
reduced to JSON format can be stored.

Other attributes
----------------

The HDF5 and netCDF file formats support variable attributes beyond
name, type, and shape.  These attribute metadata are generally not
useful outside a very narrow software suite for which they were
designed, but may provide helpful documentation when first opening a
category of file.  Therefore, qnd supports variable attributes for
backend formats which support them.  In qnd, all attributes belong to
the `QGroup` of the parent.  Thus, `QList` elements may not have
attributes (which is irrelevant since neither HDF5 nor netCDF has
native support for list objects)::

  fattrs = f.attrs()
  attrs = fattrs.x  # or fattrs['x'], attributes of f.x
  attrs = fattrs._  # or fattrs[''], attributes of f itself
  value = attrs.aname  # or attrs['aname'] value of attribute or None
  attrs.aname = value  # declare and set attribute
  attrs.aname = dtype, shape, value  # convert value to dtype and shape
  anames = list(fattrs.x)  # names of attributes of f.x
  if aname in fattrs.x: do_something
  for aname in fattrs.x: do_something
  for aname, avalue in fattrs.x.items(): do_something

Attribute values may not be dict or non-array-like lists.  Also, the
attribute names 'dtype', 'shape', 'size', 'ndim', and 'sshape' will
always return the corresponding properties of the item, even though
they are not stored as variable attributes and are not actually present
in the `attrs` mapping objects.
