# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_Walksat')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_Walksat')
    _Walksat = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Walksat', [dirname(__file__)])
        except ImportError:
            import _Walksat
            return _Walksat
        try:
            _mod = imp.load_module('_Walksat', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _Walksat = swig_import_helper()
    del swig_import_helper
else:
    import _Walksat
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        object.__setattr__(self, name, value)
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_method(set):
    def set_attr(self, name, value):
        if (name == "thisown"):
            return self.this.own(value)
        if hasattr(self, name) or (name == "this"):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr


import Numberjack.solvers.SatWrapper
class WalksatSolver(Numberjack.solvers.SatWrapper.SatWrapperSolver):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self):
        this = _Walksat.new_WalksatSolver()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Walksat.delete_WalksatSolver
    __del__ = lambda self: None

    def truth_value(self, x):
        return _Walksat.WalksatSolver_truth_value(self, x)

    def initialise(self, *args):
        return _Walksat.WalksatSolver_initialise(self, *args)

    def solve(self):
        return _Walksat.WalksatSolver_solve(self)

    def solveAndRestart(self, *args):
        return _Walksat.WalksatSolver_solveAndRestart(self, *args)

    def store_solution(self):
        return _Walksat.WalksatSolver_store_solution(self)

    def setFailureLimit(self, cutoff):
        return _Walksat.WalksatSolver_setFailureLimit(self, cutoff)

    def setNodeLimit(self, cutoff):
        return _Walksat.WalksatSolver_setNodeLimit(self, cutoff)

    def setTimeLimit(self, cutoff):
        return _Walksat.WalksatSolver_setTimeLimit(self, cutoff)

    def setRestartLimit(self, cutoff):
        return _Walksat.WalksatSolver_setRestartLimit(self, cutoff)

    def setVerbosity(self, degree):
        return _Walksat.WalksatSolver_setVerbosity(self, degree)

    def setRandomSeed(self, seed):
        return _Walksat.WalksatSolver_setRandomSeed(self, seed)

    def is_sat(self):
        return _Walksat.WalksatSolver_is_sat(self)

    def printStatistics(self):
        return _Walksat.WalksatSolver_printStatistics(self)

    def getTime(self):
        return _Walksat.WalksatSolver_getTime(self)
WalksatSolver_swigregister = _Walksat.WalksatSolver_swigregister
WalksatSolver_swigregister(WalksatSolver)


import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Walksat", "SatWrapper", model, X, FD, clause_limit, encoding)



