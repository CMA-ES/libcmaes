"""usage (yet versatile): 

    import lcmaes_interface as lci 
    
    # setup input parameters
    myfun = lambda x: sum([xi**2 for xi in x])  # myfun accepts a list of numbers as input
    x0 = [2.1] * 10
    sigma0 = 0.1
    
    # run optimization via lci
    res = lci.lcmaes.pcmaes(lci.to_fitfunc(myfun), lci.to_params(x0, sigma0))
    lci.plot()
    
Details: for the time being `to_params` is based on `lcmaes.make_simple_parameters`, 
but that might change in future to expose more parameters. 

"""

import lcmaes
import cma_multiplt as cmaplt
outfile_current = 'lcmaes.dat'

def to_params(x0, sigma0, strategy="acmaes", outfile=None):
    """return parameter object instance for ``lcmaes.pcmaes``"""
    if outfile:
        global outfile_current
        outfile_current = outfile
    p = cmaes.make_simple_parameters(x0, sigma0)
    p.set_str_algo("acmaes")
    p.set_fplot(outfile_current)
    return p

def to_fitfunc(f):
    """return function from callable `f`, where `f` accepts a list of numbers as input."""
    return cmaes.fitfunc_pbf.from_callable(lambda x, n: f(x))
    
def plot(outfile=None):
    cmaplt.plot(outfile if outfile else outfile_current)
    
