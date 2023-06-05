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
        mname = '.'.join((pkg, '_Toulbar2')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_Toulbar2')
    _Toulbar2 = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Toulbar2', [dirname(__file__)])
        except ImportError:
            import _Toulbar2
            return _Toulbar2
        try:
            _mod = imp.load_module('_Toulbar2', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _Toulbar2 = swig_import_helper()
    del swig_import_helper
else:
    import _Toulbar2
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


class Toulbar2_Expression(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    nbj_ident = _swig_property(_Toulbar2.Toulbar2_Expression_nbj_ident_get, _Toulbar2.Toulbar2_Expression_nbj_ident_set)
    _solver = _swig_property(_Toulbar2.Toulbar2_Expression__solver_get, _Toulbar2.Toulbar2_Expression__solver_set)
    _name = _swig_property(_Toulbar2.Toulbar2_Expression__name_get, _Toulbar2.Toulbar2_Expression__name_set)
    _iinf = _swig_property(_Toulbar2.Toulbar2_Expression__iinf_get, _Toulbar2.Toulbar2_Expression__iinf_set)
    _isup = _swig_property(_Toulbar2.Toulbar2_Expression__isup_get, _Toulbar2.Toulbar2_Expression__isup_set)
    _size = _swig_property(_Toulbar2.Toulbar2_Expression__size_get, _Toulbar2.Toulbar2_Expression__size_set)
    _domain = _swig_property(_Toulbar2.Toulbar2_Expression__domain_get, _Toulbar2.Toulbar2_Expression__domain_set)
    _wcspIndex = _swig_property(_Toulbar2.Toulbar2_Expression__wcspIndex_get, _Toulbar2.Toulbar2_Expression__wcspIndex_set)

    def has_been_added(self):
        return _Toulbar2.Toulbar2_Expression_has_been_added(self)

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_Expression(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_Expression
    __del__ = lambda self: None

    def getVariableId(self):
        return _Toulbar2.Toulbar2_Expression_getVariableId(self)

    def next(self, v):
        return _Toulbar2.Toulbar2_Expression_next(self, v)

    def get_value(self):
        return _Toulbar2.Toulbar2_Expression_get_value(self)

    def get_size(self):
        return _Toulbar2.Toulbar2_Expression_get_size(self)

    def get_min(self):
        return _Toulbar2.Toulbar2_Expression_get_min(self)

    def get_max(self):
        return _Toulbar2.Toulbar2_Expression_get_max(self)

    def contain(self, v):
        return _Toulbar2.Toulbar2_Expression_contain(self, v)

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_Expression_add(self, solver, top_level)
Toulbar2_Expression_swigregister = _Toulbar2.Toulbar2_Expression_swigregister
Toulbar2_Expression_swigregister(Toulbar2_Expression)
cvar = _Toulbar2.cvar
MAXCOST = cvar.MAXCOST

class Toulbar2_IntVar(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_IntVar(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_IntVar
    __del__ = lambda self: None
Toulbar2_IntVar_swigregister = _Toulbar2.Toulbar2_IntVar_swigregister
Toulbar2_IntVar_swigregister(Toulbar2_IntVar)

class Toulbar2_AllDiff(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_AllDiff(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_AllDiff_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_AllDiff
    __del__ = lambda self: None
Toulbar2_AllDiff_swigregister = _Toulbar2.Toulbar2_AllDiff_swigregister
Toulbar2_AllDiff_swigregister(Toulbar2_AllDiff)

class Toulbar2_Gcc(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, vals, lb_card, ub_card):
        this = _Toulbar2.new_Toulbar2_Gcc(vars, vals, lb_card, ub_card)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_Gcc_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_Gcc
    __del__ = lambda self: None
Toulbar2_Gcc_swigregister = _Toulbar2.Toulbar2_Gcc_swigregister
Toulbar2_Gcc_swigregister(Toulbar2_Gcc)

class Toulbar2_PostNullary(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, cost):
        this = _Toulbar2.new_Toulbar2_PostNullary(cost)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostNullary_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostNullary
    __del__ = lambda self: None
Toulbar2_PostNullary_swigregister = _Toulbar2.Toulbar2_PostNullary_swigregister
Toulbar2_PostNullary_swigregister(Toulbar2_PostNullary)

class Toulbar2_PostUnary(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, var, costs):
        this = _Toulbar2.new_Toulbar2_PostUnary(var, costs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostUnary_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostUnary
    __del__ = lambda self: None
Toulbar2_PostUnary_swigregister = _Toulbar2.Toulbar2_PostUnary_swigregister
Toulbar2_PostUnary_swigregister(Toulbar2_PostUnary)

class Toulbar2_PostBinary(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, var1, var2, costs):
        this = _Toulbar2.new_Toulbar2_PostBinary(var1, var2, costs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostBinary_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostBinary
    __del__ = lambda self: None
Toulbar2_PostBinary_swigregister = _Toulbar2.Toulbar2_PostBinary_swigregister
Toulbar2_PostBinary_swigregister(Toulbar2_PostBinary)

class Toulbar2_PostTernary(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, costs):
        this = _Toulbar2.new_Toulbar2_PostTernary(vars, costs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostTernary_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostTernary
    __del__ = lambda self: None
Toulbar2_PostTernary_swigregister = _Toulbar2.Toulbar2_PostTernary_swigregister
Toulbar2_PostTernary_swigregister(Toulbar2_PostTernary)

class Toulbar2_PostNary(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, arity, _defcost, values, costs):
        this = _Toulbar2.new_Toulbar2_PostNary(vars, arity, _defcost, values, costs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostNary_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostNary
    __del__ = lambda self: None
Toulbar2_PostNary_swigregister = _Toulbar2.Toulbar2_PostNary_swigregister
Toulbar2_PostNary_swigregister(Toulbar2_PostNary)

class Toulbar2_PostWSum(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, arity, semantics, baseCost, comparator, rightRes):
        this = _Toulbar2.new_Toulbar2_PostWSum(vars, arity, semantics, baseCost, comparator, rightRes)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostWSum_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostWSum
    __del__ = lambda self: None
Toulbar2_PostWSum_swigregister = _Toulbar2.Toulbar2_PostWSum_swigregister
Toulbar2_PostWSum_swigregister(Toulbar2_PostWSum)

class Toulbar2_PostWVarSum(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, arity, semantics, baseCost, comparator):
        this = _Toulbar2.new_Toulbar2_PostWVarSum(vars, arity, semantics, baseCost, comparator)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostWVarSum_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostWVarSum
    __del__ = lambda self: None
Toulbar2_PostWVarSum_swigregister = _Toulbar2.Toulbar2_PostWVarSum_swigregister
Toulbar2_PostWVarSum_swigregister(Toulbar2_PostWVarSum)

class Toulbar2_PostWAmong(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_PostWAmong(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostWAmong_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostWAmong
    __del__ = lambda self: None
Toulbar2_PostWAmong_swigregister = _Toulbar2.Toulbar2_PostWAmong_swigregister
Toulbar2_PostWAmong_swigregister(Toulbar2_PostWAmong)

class Toulbar2_Regular(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_Regular(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_Regular_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_Regular
    __del__ = lambda self: None
Toulbar2_Regular_swigregister = _Toulbar2.Toulbar2_Regular_swigregister
Toulbar2_Regular_swigregister(Toulbar2_Regular)

class Toulbar2_Same(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_Same(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_Same_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_Same
    __del__ = lambda self: None
Toulbar2_Same_swigregister = _Toulbar2.Toulbar2_Same_swigregister
Toulbar2_Same_swigregister(Toulbar2_Same)

class Toulbar2_PostWSameGcc(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, vals, lb_card, ub_card, type, semantics, baseCost):
        this = _Toulbar2.new_Toulbar2_PostWSameGcc(vars, vals, lb_card, ub_card, type, semantics, baseCost)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostWSameGcc_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostWSameGcc
    __del__ = lambda self: None
Toulbar2_PostWSameGcc_swigregister = _Toulbar2.Toulbar2_PostWSameGcc_swigregister
Toulbar2_PostWSameGcc_swigregister(Toulbar2_PostWSameGcc)

class Toulbar2_PostWOverlap(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, vars, arity, semantics, baseCost, comparator, rightRes):
        this = _Toulbar2.new_Toulbar2_PostWOverlap(vars, arity, semantics, baseCost, comparator, rightRes)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_PostWOverlap_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_PostWOverlap
    __del__ = lambda self: None
Toulbar2_PostWOverlap_swigregister = _Toulbar2.Toulbar2_PostWOverlap_swigregister
Toulbar2_PostWOverlap_swigregister(Toulbar2_PostWOverlap)

class Toulbar2_Table(Toulbar2_Expression):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _Toulbar2.new_Toulbar2_Table(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def add(self, solver, top_level):
        return _Toulbar2.Toulbar2_Table_add(self, solver, top_level)
    __swig_destroy__ = _Toulbar2.delete_Toulbar2_Table
    __del__ = lambda self: None
Toulbar2_Table_swigregister = _Toulbar2.Toulbar2_Table_swigregister
Toulbar2_Table_swigregister(Toulbar2_Table)

class Toulbar2Solver(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    solver = _swig_property(_Toulbar2.Toulbar2Solver_solver_get, _Toulbar2.Toulbar2Solver_solver_set)
    wcsp = _swig_property(_Toulbar2.Toulbar2Solver_wcsp_get, _Toulbar2.Toulbar2Solver_wcsp_set)
    upperbound = _swig_property(_Toulbar2.Toulbar2Solver_upperbound_get, _Toulbar2.Toulbar2Solver_upperbound_set)
    optimum = _swig_property(_Toulbar2.Toulbar2Solver_optimum_get, _Toulbar2.Toulbar2Solver_optimum_set)
    costshift = _swig_property(_Toulbar2.Toulbar2Solver_costshift_get, _Toulbar2.Toulbar2Solver_costshift_set)
    unsatisfiable = _swig_property(_Toulbar2.Toulbar2Solver_unsatisfiable_get, _Toulbar2.Toulbar2Solver_unsatisfiable_set)
    interrupted = _swig_property(_Toulbar2.Toulbar2Solver_interrupted_get, _Toulbar2.Toulbar2Solver_interrupted_set)
    solution = _swig_property(_Toulbar2.Toulbar2Solver_solution_get, _Toulbar2.Toulbar2Solver_solution_set)

    def __init__(self):
        this = _Toulbar2.new_Toulbar2Solver()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Toulbar2.delete_Toulbar2Solver
    __del__ = lambda self: None

    def add(self, arg):
        return _Toulbar2.Toulbar2Solver_add(self, arg)

    def initialise(self, *args):
        return _Toulbar2.Toulbar2Solver_initialise(self, *args)

    def propagate(self):
        return _Toulbar2.Toulbar2Solver_propagate(self)

    def solve(self):
        return _Toulbar2.Toulbar2Solver_solve(self)

    def solveAndRestart(self, policy=0, base=32, factor=1.3333333, decay=0.0, reinit=-1):
        return _Toulbar2.Toulbar2Solver_solveAndRestart(self, policy, base, factor, decay, reinit)

    def store_solution(self):
        return _Toulbar2.Toulbar2Solver_store_solution(self)

    def is_opt(self):
        return _Toulbar2.Toulbar2Solver_is_opt(self)

    def is_sat(self):
        return _Toulbar2.Toulbar2Solver_is_sat(self)

    def is_unsat(self):
        return _Toulbar2.Toulbar2Solver_is_unsat(self)

    def getOptimum(self):
        return _Toulbar2.Toulbar2Solver_getOptimum(self)

    def getBacktracks(self):
        return _Toulbar2.Toulbar2Solver_getBacktracks(self)

    def getNodes(self):
        return _Toulbar2.Toulbar2Solver_getNodes(self)

    def getFailures(self):
        return _Toulbar2.Toulbar2Solver_getFailures(self)

    def getPropags(self):
        return _Toulbar2.Toulbar2Solver_getPropags(self)

    def getTime(self):
        return _Toulbar2.Toulbar2Solver_getTime(self)

    def getChecks(self):
        return _Toulbar2.Toulbar2Solver_getChecks(self)

    def printStatistics(self):
        return _Toulbar2.Toulbar2Solver_printStatistics(self)

    def getNumVariables(self):
        return _Toulbar2.Toulbar2Solver_getNumVariables(self)

    def getNumConstraints(self):
        return _Toulbar2.Toulbar2Solver_getNumConstraints(self)

    def printPython(self):
        return _Toulbar2.Toulbar2Solver_printPython(self)

    def setHeuristic(self, var_heuristic, val_heuristic, rand):
        return _Toulbar2.Toulbar2Solver_setHeuristic(self, var_heuristic, val_heuristic, rand)

    def setFailureLimit(self, cutoff):
        return _Toulbar2.Toulbar2Solver_setFailureLimit(self, cutoff)

    def setNodeLimit(self, cutoff):
        return _Toulbar2.Toulbar2Solver_setNodeLimit(self, cutoff)

    def setTimeLimit(self, cutoff):
        return _Toulbar2.Toulbar2Solver_setTimeLimit(self, cutoff)

    def setRandomized(self, degree):
        return _Toulbar2.Toulbar2Solver_setRandomized(self, degree)

    def setRandomSeed(self, seed):
        return _Toulbar2.Toulbar2Solver_setRandomSeed(self, seed)

    def setVerbosity(self, degree):
        return _Toulbar2.Toulbar2Solver_setVerbosity(self, degree)

    def debug(self, debug):
        return _Toulbar2.Toulbar2Solver_debug(self, debug)

    def writeSolution(self, write):
        return _Toulbar2.Toulbar2Solver_writeSolution(self, write)

    def showSolutions(self, show):
        return _Toulbar2.Toulbar2Solver_showSolutions(self, show)

    def dumpWCSP(self, level, problem):
        return _Toulbar2.Toulbar2Solver_dumpWCSP(self, level, problem)

    def nopre(self):
        return _Toulbar2.Toulbar2Solver_nopre(self)

    def updateUb(self, newUb):
        return _Toulbar2.Toulbar2Solver_updateUb(self, newUb)

    def lds(self, maxlds):
        return _Toulbar2.Toulbar2Solver_lds(self, maxlds)

    def restart(self, maxrestarts):
        return _Toulbar2.Toulbar2Solver_restart(self, maxrestarts)

    def hbfs(self, hbfsgloballimit):
        return _Toulbar2.Toulbar2Solver_hbfs(self, hbfsgloballimit)

    def hbfsAlpha(self, hbfsalpha):
        return _Toulbar2.Toulbar2Solver_hbfsAlpha(self, hbfsalpha)

    def hbfsBeta(self, hbfsbeta):
        return _Toulbar2.Toulbar2Solver_hbfsBeta(self, hbfsbeta)

    def hbfsOpenNodeLimit(self, openlimit):
        return _Toulbar2.Toulbar2Solver_hbfsOpenNodeLimit(self, openlimit)

    def lcLevel(self, level):
        return _Toulbar2.Toulbar2Solver_lcLevel(self, level)

    def QueueComplexity(self, queue):
        return _Toulbar2.Toulbar2Solver_QueueComplexity(self, queue)

    def allSolutions(self, sol):
        return _Toulbar2.Toulbar2Solver_allSolutions(self, sol)

    def approximateCountingBTD(self, aproxim):
        return _Toulbar2.Toulbar2Solver_approximateCountingBTD(self, aproxim)

    def binaryBranching(self, boost):
        return _Toulbar2.Toulbar2Solver_binaryBranching(self, boost)

    def staticVariableOrdering(self, staticOrdering):
        return _Toulbar2.Toulbar2Solver_staticVariableOrdering(self, staticOrdering)

    def lastConflict(self, last):
        return _Toulbar2.Toulbar2Solver_lastConflict(self, last)

    def dichotomicBranching(self, dicho):
        return _Toulbar2.Toulbar2Solver_dichotomicBranching(self, dicho)

    def sortDomains(self, sort):
        return _Toulbar2.Toulbar2Solver_sortDomains(self, sort)

    def weightedDegree(self, wDegree):
        return _Toulbar2.Toulbar2Solver_weightedDegree(self, wDegree)

    def weightedTightness(self, wTight):
        return _Toulbar2.Toulbar2Solver_weightedTightness(self, wTight)

    def variableEliminationOrdering(self, order):
        return _Toulbar2.Toulbar2Solver_variableEliminationOrdering(self, order)

    def nbDecisionVars(self, nbDecision):
        return _Toulbar2.Toulbar2Solver_nbDecisionVars(self, nbDecision)

    def partialAssign(self, certificate):
        return _Toulbar2.Toulbar2Solver_partialAssign(self, certificate)

    def elimDegree(self, degree):
        return _Toulbar2.Toulbar2Solver_elimDegree(self, degree)

    def elimDegree_preprocessing(self, degree_prepoc):
        return _Toulbar2.Toulbar2Solver_elimDegree_preprocessing(self, degree_prepoc)

    def elimSpaceMaxMB(self, size):
        return _Toulbar2.Toulbar2Solver_elimSpaceMaxMB(self, size)

    def costfuncSeparate(self, separate):
        return _Toulbar2.Toulbar2Solver_costfuncSeparate(self, separate)

    def deadEndElimination(self, level):
        return _Toulbar2.Toulbar2Solver_deadEndElimination(self, level)

    def preprocessTernaryRPC(self, size):
        return _Toulbar2.Toulbar2Solver_preprocessTernaryRPC(self, size)

    def preprocessFunctional(self, func):
        return _Toulbar2.Toulbar2Solver_preprocessFunctional(self, func)

    def preprocessNary(self, maxnary):
        return _Toulbar2.Toulbar2Solver_preprocessNary(self, maxnary)

    def btdMode(self, mode):
        return _Toulbar2.Toulbar2Solver_btdMode(self, mode)

    def splitClusterMaxSize(self, size):
        return _Toulbar2.Toulbar2Solver_splitClusterMaxSize(self, size)

    def maxSeparatorSize(self, size):
        return _Toulbar2.Toulbar2Solver_maxSeparatorSize(self, size)

    def minProperVarSize(self, size):
        return _Toulbar2.Toulbar2Solver_minProperVarSize(self, size)

    def boostingBTD(self, boost):
        return _Toulbar2.Toulbar2Solver_boostingBTD(self, boost)

    def btdRootCluster(self, rCluster):
        return _Toulbar2.Toulbar2Solver_btdRootCluster(self, rCluster)

    def btdSubTree(self, sTree):
        return _Toulbar2.Toulbar2Solver_btdSubTree(self, sTree)

    def vac(self, depth):
        return _Toulbar2.Toulbar2Solver_vac(self, depth)

    def vacValueHeuristic(self, vacVal):
        return _Toulbar2.Toulbar2Solver_vacValueHeuristic(self, vacVal)

    def costThreshold(self, cost):
        return _Toulbar2.Toulbar2Solver_costThreshold(self, cost)

    def costThresholdPre(self, cost):
        return _Toulbar2.Toulbar2Solver_costThresholdPre(self, cost)

    def costMultiplier(self, cost):
        return _Toulbar2.Toulbar2Solver_costMultiplier(self, cost)

    def singletonConsistency(self, singleCons):
        return _Toulbar2.Toulbar2Solver_singletonConsistency(self, singleCons)

    def minsumDiffusion(self, min):
        return _Toulbar2.Toulbar2Solver_minsumDiffusion(self, min)

    def incop(self, cmd):
        return _Toulbar2.Toulbar2Solver_incop(self, cmd)
Toulbar2Solver_swigregister = _Toulbar2.Toulbar2Solver_swigregister
Toulbar2Solver_swigregister(Toulbar2Solver)

class Toulbar2ExpArray(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self):
        this = _Toulbar2.new_Toulbar2ExpArray()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Toulbar2.delete_Toulbar2ExpArray
    __del__ = lambda self: None

    def size(self):
        return _Toulbar2.Toulbar2ExpArray_size(self)

    def add(self, arg):
        return _Toulbar2.Toulbar2ExpArray_add(self, arg)

    def erase(self, i):
        return _Toulbar2.Toulbar2ExpArray_erase(self, i)

    def get_item(self, i):
        return _Toulbar2.Toulbar2ExpArray_get_item(self, i)

    def set_item(self, i, item):
        return _Toulbar2.Toulbar2ExpArray_set_item(self, i, item)
Toulbar2ExpArray_swigregister = _Toulbar2.Toulbar2ExpArray_swigregister
Toulbar2ExpArray_swigregister(Toulbar2ExpArray)

class Toulbar2IntArray(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self):
        this = _Toulbar2.new_Toulbar2IntArray()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Toulbar2.delete_Toulbar2IntArray
    __del__ = lambda self: None

    def size(self):
        return _Toulbar2.Toulbar2IntArray_size(self)

    def add(self, arg):
        return _Toulbar2.Toulbar2IntArray_add(self, arg)

    def erase(self, i):
        return _Toulbar2.Toulbar2IntArray_erase(self, i)

    def get_item(self, i):
        return _Toulbar2.Toulbar2IntArray_get_item(self, i)

    def set_item(self, i, item):
        return _Toulbar2.Toulbar2IntArray_set_item(self, i, item)
Toulbar2IntArray_swigregister = _Toulbar2.Toulbar2IntArray_swigregister
Toulbar2IntArray_swigregister(Toulbar2IntArray)

class Toulbar2DoubleArray(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self):
        this = _Toulbar2.new_Toulbar2DoubleArray()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _Toulbar2.delete_Toulbar2DoubleArray
    __del__ = lambda self: None

    def size(self):
        return _Toulbar2.Toulbar2DoubleArray_size(self)

    def add(self, arg):
        return _Toulbar2.Toulbar2DoubleArray_add(self, arg)

    def erase(self, i):
        return _Toulbar2.Toulbar2DoubleArray_erase(self, i)

    def get_item(self, i):
        return _Toulbar2.Toulbar2DoubleArray_get_item(self, i)

    def set_item(self, i, item):
        return _Toulbar2.Toulbar2DoubleArray_set_item(self, i, item)
Toulbar2DoubleArray_swigregister = _Toulbar2.Toulbar2DoubleArray_swigregister
Toulbar2DoubleArray_swigregister(Toulbar2DoubleArray)


import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Toulbar2", "Toulbar2", model, X, FD, clause_limit, encoding)



