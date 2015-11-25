import lcmaes

# input parameters for a 10-D problem
x = [10]*10
lambda_ = 10 # lambda is a reserved keyword in python, using lambda_ instead.
seed = 0 # 0 for seed auto-generated within the lib.
sigma = 0.1
p = lcmaes.make_simple_parameters(x,sigma,lambda_,seed)
p.set_str_algo("acmaes")

# objective function.
def nfitfunc(x,n):
    val = 0.0
    for i in range(0,n):
        val += x[i]*x[i]
    return val

# generate a function object
objfunc = lcmaes.fitfunc_pbf.from_callable(nfitfunc);

# pass the function and parameter to cmaes, run optimization and collect solution object.
cmasols = lcmaes.pcmaes(objfunc,p)

# collect and inspect results
bcand = cmasols.best_candidate()
bx = lcmaes.get_candidate_x(bcand)
print("best x=",bx)
print("distribution mean=",lcmaes.get_solution_xmean(cmasols))
cov = lcmaes.get_solution_cov(cmasols) # numpy array
print("cov=",cov)
print("elapsed time=",cmasols.elapsed_time(),"ms")
