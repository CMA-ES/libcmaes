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
