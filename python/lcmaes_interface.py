"""usage (yet versatile):

    import lcmaes_interface as lci

    # setup input parameters
    def myfun(x): return sum([xi**2 for xi in x])  # myfun accepts a list of numbers as input
    dimension = 10 
    x0 = [2.1] * dimension
    sigma0 = 0.1

    # run optimization via lci
    res = lci.pcmaes(lci.to_fitfunc(myfun),
                     lci.to_params(x0, sigma0,  # all remaining parameters are optional
                         str_algo=b'aipop',  # b=bytes, because unicode fails
                         lbounds=[-5] * dimension, ubounds=[5] * dimension,
                         restarts=2,  # means 2 runs, i.e. one restart
                         )
                    ) 
    lci.plot()  # plot from file set in lci.to_params

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import lcmaes
import cma_multiplt as cmaplt
fplot_current = b'lcmaes.dat'

def to_params(x0, sigma0, str_algo=b'acmaes', fplot=None, lbounds=None, ubounds=None, scaling=False, vscaling=None, vshift=None, **kwargs):
    """return parameter object instance for `lcmaes.pcmaes`.

    Keys in `kwargs` must correspond to `name` in `set_name` attributes
    of `lcmaes.CMAParameters`, e.g. ``ftarget=1e-7``.

    Details: when `fplot is None` (default), the default output filename
    is used.
    """
    has_bounds = not lbounds==None and not ubounds == None
    p = None
    if has_bounds:
        if scaling==False:
            gp = lcmaes.make_genopheno_pwqb(lbounds,ubounds,len(ubounds))
            p = lcmaes.make_parameters_pwqb(x0,sigma0,gp)
        else:
            gp = lcmaes.make_genopheno_pwqb_ls(lbounds,ubounds,len(ubounds))
            p = lcmaes.make_parameters_pwqb_ls(x0,sigma0,gp,-1,0)
    else:
        if vscaling is None:
            p = lcmaes.make_simple_parameters(x0, sigma0)
        else:
            gp = lcmaes.make_genopheno_ls(vscaling,vshift)
            p = lcmaes.make_parameters_ls(x0,sigma0,gp)
    p.set_str_algo(str_algo)
    if fplot and fplot != True:  # then fplot must be filename
        global fplot_current
        fplot_current = fplot
    if fplot or fplot is None:  # 0 or False or '' or "" prevents writing
        p.set_fplot(fplot_current)
    for key, val in kwargs.items():
        setter = "set_" + key
        if not hasattr(p, setter):
            raise ValueError(setter + " is not known as method of CMAParameters")
        getattr(p, setter)(val)  # call setter with value
    return p

def pcmaes(fitfunc,p):
    has_bounds = isinstance(p,lcmaes.CMAParametersPB) or isinstance(p,lcmaes.CMAParametersPBS)
    has_scaling = isinstance(p,lcmaes.CMAParametersNBS) or isinstance(p,lcmaes.CMAParametersPBS)
    if not has_bounds:
        if not has_scaling:
            return lcmaes.pcmaes(fitfunc,p)
        else:
            return lcmaes.pcmaes_ls(fitfunc,p)
    else:
        if not has_scaling:
            return lcmaes.pcmaes_pwqb(fitfunc,p)
        else:
            return lcmaes.pcmaes_pwqb_ls(fitfunc,p)

def to_fitfunc(f):
    """return function for lcmaes from callable `f`, where `f` accepts a list of numbers as input."""
    return lcmaes.fitfunc_pbf.from_callable(lambda x, n: f(x))

def plot(file=None):
    cmaplt.plot(file if file else fplot_current)
    cmaplt.pylab.ioff()
    cmaplt.pylab.show()

