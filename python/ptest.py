import lcmaes

x = [2,2]
#p = lcmaes.CMAParametersNB()
olambda = 10
seed = 0
sigma = 0.1
p = lcmaes.make_parameters(x,sigma,olambda,seed)
#p.set_noisy()

def nfitfunc(x,n):
    val = 0.0
    for i in range(0,n):
        val += x[i]*x[i]
    return val
        
objfunc = lcmaes.fitfunc_pbf.from_callable(nfitfunc);
cmasols = lcmaes.pcmaes(objfunc,p)
bcand = cmasols.best_candidate()
bx = lcmaes.get_candidate_x(bcand)
print "best x=",bx
print "distribution mean=",lcmaes.get_solution_xmean(cmasols)
cov = lcmaes.get_solution_cov(cmasols) # numpy array
print "cov=",cov
print "elapsed time=",cmasols.elapsed_time(),"ms"
