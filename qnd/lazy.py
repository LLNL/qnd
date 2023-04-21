"""Lazy file type identification and permission filters."""

from __future__ import absolute_import

from glob import glob
from os.path import expanduser, expandvars
from os import access, R_OK
from importlib import import_module
import sys

PY2 = sys.version_info < (3,)
if not PY2:
    basestring = str


def permfilter(*names):
    """Given a glob string or list of filenames, exclude unreadable ones.

    Parameters
    ----------
    names : str or list of str
       An individual filename or a list of filenames, which may include
       globbing characters like '*' or '?'.
       Accepts any number of such arguments.

    Returns
    -------
    names : list of str
       The subset of all the specified names for which you have read
       permission (and which exist).  The list may be empty.

    Globbed results are sorted first by length, then alphabetically,
    so that, for example, name100 follows name99.  Note, however, that
    name000 also follows name99.

    """
    args = []
    for name in names:
        if isinstance(name, basestring):
            args.append(name)
        else:
            args.extend(name)
    names = []
    for name in args:
        names.extend(sorted(glob(expanduser(expandvars(name))),
                            key=lambda f: (len(f), f)))
    return filter(lambda f: access(f, R_OK), names)


def openb(filename, mode="r", auto=1, **kwargs):
    """Open PDB or HDF5 or netCDF or some othe binary file known to QnD.

    This imports the required package as necessary.  See
    ``qnd.pdbf.openpdb`` or ``qnd.h5f.openh5`` or some other actual
    open function for argument details.

    """
    name = filename
    if isinstance(name, basestring):
        namex = expanduser(expandvars(name))
        name = list((sorted(glob(namex), key=lambda f: (len(f), f))))
        if not name and mode.startswith("w"):
            name = [namex];
    if not name:
        raise ValueError("no such file(s): {}".format(filename))
    name = name[0]
    pkg = None
    ext = name.rsplit(".", 1)
    ext = ext[1].lower() if len(ext) > 1 else ""
    if ext in ["h5", "hdf", "hdf5"]:
        pkg = "h5"
    elif ext == "pdb":
        pkg = "pdb"
    elif ext == "nc":
        pkg = "nc"
    if not pkg:
        with open(name, "rb") as f:
            h = f.read(1024)
            if h.startswith(b"\x89HDF\r\n\x1a\n"):
                pkg = "h5"
            elif h.startswith(b"!<<PDB:II>>!") or h.startswith(b"!<><PDB><>!"):
                pkg = "pdb"
            elif h.startswith(b"CDF\x01") or h.startswith(b"CDF\x02"):
                pkg = "nc"
    if not pkg:
        raise ValueError("unknown file type: {}".format(name))
    module = import_module("." + pkg + "f", __package__ or "qnd")
    opener = getattr(module, "open" + pkg)
    return opener(filename, mode, auto, **kwargs)
