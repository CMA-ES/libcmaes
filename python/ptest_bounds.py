import lcmaes

# input parameters for a 10-D problem
x = [3]*10
lambda_ = 10 # lambda is a reserved keyword in python, using lambda_ instead.
seed = 0 # 0 for seed auto-generated within the lib.
sigma = 0.1
lbounds = [-4]*10
ubounds = [4]*10
gp = lcmaes.make_genopheno_pwqb(lbounds,ubounds,10)
p = lcmaes.make_parameters_pwqb(x,sigma,gp,lambda_,seed)

# objective function.
def nfitfunc(x,n):
    val = 0.0
    for i in range(0,n):
        val += x[i]*x[i]
    return val

# generate a function object
objfunc = lcmaes.fitfunc_pbf.from_callable(nfitfunc);

# pass the function and parameter to cmaes, run optimization and collect solution object.
cmasols = lcmaes.pcmaes_pwqb(objfunc,p)

# collect and inspect results
bcand = cmasols.best_candidate()
bcand_pheno = lcmaes.get_best_candidate_pheno(cmasols,gp)
bx = lcmaes.get_candidate_x(bcand)
print("best x=",bx)
print('best x in phenotype space=',bcand_pheno)
print("distribution mean=",lcmaes.get_solution_xmean(cmasols))
cov = lcmaes.get_solution_cov(cmasols) # numpy array
print("cov=",cov)
print("elapsed time=",cmasols.elapsed_time(),"ms")
