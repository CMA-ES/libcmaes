import lcmaes
import cma_multiplt as lcmaplt

# input parameters for a 10-D problem
x = [10]*10
sigma = 0.1
p = lcmaes.make_simple_parameters(x,sigma)
p.set_str_algo("acmaes")
p.set_fplot('simple.dat')

# objective function.
def nfitfunc(x,n):
    val = 0.0
    for i in range(0,n):
        val += x[i]*x[i]
    return val

# pass the function and parameters to cmaes, run optimization and collect solution object.
cmasols = lcmaes.pcmaes(lcmaes.fitfunc_pbf.from_callable(nfitfunc),p)

# visualize results
lcmaplt.plot('simple.dat')
