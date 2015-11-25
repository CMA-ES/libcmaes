import lcmaes
import cma_multiplt as lcmaplt

# input parameters for a 10-D problem
x = [10.0] * 10
sigma = 0.1
outfile = 'simple.dat'

p = lcmaes.make_simple_parameters(x, sigma)
p.set_str_algo("acmaes")
p.set_fplot(outfile)

# objective function.
def nfitfunc(x, n):
    assert len(x) == n  # should not be necessary
    return sum([xi**2 for xi in x])
    
# pass the function and parameters to cmaes, run optimization and collect solution object.
cmasols = lcmaes.pcmaes(lcmaes.fitfunc_pbf.from_callable(nfitfunc), p)

# visualize results
lcmaplt.plot(outfile)
lcmaplt.pylab.ioff()
lcmaplt.pylab.show()

if __name__ == "__main__":
    msg = '  --- press return to continue --- '
    try: 
        raw_input(msg) 
    except NameError: 
        input(msg)
