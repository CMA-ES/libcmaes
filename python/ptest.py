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
    for i in range(0,n-1):
        val += x[i]*x[i]
    return val
        
objfunc = lcmaes.fitfunc_pbf.from_callable(nfitfunc);
lcmaes.pcmaes(objfunc,p)
