"""An items-as-attributes dict."""
from __future__ import absolute_import

import sys
# from collections import Mapping

PY2 = sys.version_info < (3,)


class ItemsAreAttrs(object):
    """Mix-in class for QArray, QGroup, or QList, and also ADict."""
    __slots__ = ()

    def __getattr__(self, name):
        # Mixing __getattr__ with properties, as we do here, can lead to
        # very obscure errors, because the way __getattr__ works is to
        # simply try to retrieve the attribute, and call __getattr__ if
        # that operation raises an exception, instead of propagating the
        # exception as usual.  When you have a property or any descriptor,
        # the act of retrieving that attribute implicitly calls a method.
        # If that property method raises an exception, the __getattr__
        # machinery interprets that to mean the attribute does not exist
        # and calls __getattr__, removing the actual faulting method from
        # the call chain and making debugging difficult.
        # Beware of errors with this __getattr__ in their calling chain;
        # they may have originated in a different error in a property
        # method!
        if name.startswith('__') and len(name) > 2 or name == 'getdoc':
            # Do not redirect dunder or ipython getdoc calls, as this
            # confuses many simple attempts at introspection.
            return super(ItemsAreAttrs, self).__getattr__(name)
        # Strip single trailing _ as an interactive convenience for the
        # problem of attribute names that match reserved words or property
        # or method names.  This is inspired by the PEP8 advice for dealing
        # with this issue.
        if name.endswith('_'):
            name = name[:-1]
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(*e.args)

    def __setattr__(self, name, value):
        if name.startswith('__') and len(name) > 2:
            super(ItemsAreAttrs, self).__setattr__(name, value)
        if name.endswith('_'):
            name = name[:-1]
        self[name] = value

    def __delattr__(self, name):
        if name.startswith('__') and len(name) > 2:
            super(ItemsAreAttrs, self).__delattr__(name)
        if name.endswith('_'):
            name = name[:-1]
        del self[name]

    def update(self, *args, **kwargs):
        """Multiple __setitem__ from positional arguments or keywords."""
        for arg in args:
            if hasattr(arg, 'keys'):  # dict-like, not list-like
                for key in arg:
                    self[key] = arg[key]
            else:
                key, value = arg
                self[key] = value
        for key in kwargs:
            self[key] = kwargs[key[:-1] if key.endswith('_') else key]


class ADict(ItemsAreAttrs, dict):
    """Subclass of dict permitting access to items as if they were attributes.

    For a ADict ad, ``ad.x`` is equivalent to ``ad['x']`` for getting,
    setting, or deleting items.  The exceptions are dict method names,
    like `keys` or `items`, syntactically illegal names, like `class`
    or `yield`, and any name beginning with `__`.

    Additionally, as a work around for some of these exceptions, ADict
    will remove a single trailing underscore from an attribute name,
    so ``ad.x_`` is also equivalent to ``ad['x']``, and you need
    ``ad.x__`` to get ``ad['x_']`` (a convention inspired by the
    similar PEP8 recommendation for syntatically illegal variable
    names).  The trailing underscore removal does not apply to names
    beginning with `__`.

    The trailing underscore removal convention applies to keywords
    passed to the constructor or to the `update` method as well.

    Use subscript syntax when a variable or expression holds an item name;
    use attribute syntax when you know the item name at parse time::

        ad[variable] = value  # value of variable is the item name
        ad.fixed = value  # 'fixed' is the item name
        value = ad.get('fixed', default)  # except to avoid KeyError

    See Also
    --------
    redict : recursively toggle between dict and ADict
    ItemsAreAttrs : mixin base class to provide this for any class

    """
    __slots__ = ()

    def __init__(self, *args, **kwargs):
        super(ADict, self).__init__(*args)
        self.update(kwargs)

    def __repr__(self):
        return "ADict(" + super(ADict, self).__repr__() + ")"


def redict(d, cls=None):
    """Recursively convert a nested dict to a nested ADict and vice versa.

    Parameters
    ----------
    d : dict or ADict instance
        A dict, possibly nested, to be converted.
    cls : dict or subclass of dict, optional
        The dict-like cls to recursively convert `d` and any sub-dicts
        into.  By default, if `d` is a `ADict`, `cls` is `dict`,
        otherwise `cls` is `ADict`, so repeated calls to `redict` toggle
        between `dict` and `ADict`.

    Returns
    -------
    dnew : dict or ADict
        A copy of `d` whose class is `cls`.  Any items which are dict
        instances are similarly copied to be `cls` instances.  Non-dict
        items are not copied unless assignment makes copies.

    """
    if cls is None:
        cls = dict if isinstance(d, ADict) else ADict
    dnew = cls(d)
    for key, value in (d.iteritems() if PY2 else d.items()):
        if hasattr(value, '__iter__') and hasattr(value, 'keys'):
            dnew[key] = redict(value, cls)
    return dnew
