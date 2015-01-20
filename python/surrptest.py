import lcmaes
import numpy as np
import sys
import math

# input parameters for a 10-D problem
x0 = 2
dim = 10
x = [x0]*dim
lambda_ = 10 # lambda is a reserved keyword in python, using lambda_ instead.
seed = 0 # 0 for seed auto-generated within the lib.
sigma = 0.1
p = lcmaes.make_simple_parameters(x,sigma,lambda_,seed)
p.set_str_algo("cmaes")
#p.set_max_iter(20)
p.set_ftarget(1e-3)
p.set_quiet(False)

# surrogate
tree_max_depth = 100000
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.svm import NuSVR

#clf = DecisionTreeRegressor(max_depth=tree_max_depth)
#clf = RandomForestRegressor(n_estimators=100)
#clf = GradientBoostingRegressor()
clf = NuSVR(C=1.0)

def encode(x,Cinv):
    xmean = np.mean(x,axis=0)
    xxmean = []
    for xi in x:
        xxmean.append(list(xi - xmean))
    nx = []
    for xi in xxmean:
        nx.append(list(Cinv.dot(xi)))
    return nx

def traindt(x,y):
    global clf

    #print "training surrogate"
    #clft = DecisionTreeRegressor(max_depth=tree_max_depth,splitter='random')
    #clft = RandomForestRegressor()
    #clft = GradientBoostingRegressor(loss='lad',n_estimators=50,learning_rate=0.3,max_depth=2)
    clft = NuSVR(C=1e6)
    clf = clft.fit(x,y)
    
def predict(x):
    global clf
#    print "predicting from surrogate"
    y = clf.predict(x)
    return y

# objective function.
def fsphere(x,n):
    val = 0.0
    for i in range(0,n):
        val += x[i]*x[i]
    return val

def rosen(x,n):
    val = 0.0
    for i in range(0,n-1):
        val += 100.0*pow(x[i+1]-x[i]*x[i],2) + pow(x[i]-1.0,2)
    return val

# generate a function object
objfunc = lcmaes.fitfunc_pbf.from_callable(rosen)

# surrogate
def ncsurrfunc(c,m):
    x = []
    y = []
    for ci in c:
        x.append(lcmaes.get_candidate_x(ci))
        y.append(ci.get_fvalue())
    #nx = encode(x,m)
    traindt(x,y)
    return 1
trainfunc = lcmaes.csurrfunc_pbf.from_callable(ncsurrfunc)

def nsurrfunc(c,m):
    x = []
    for ci in c:
        x.append(lcmaes.get_candidate_x(ci))
    #nx = encode(x,m)
    y = predict(x)
    i = 0
    for ci in c:
        ci.set_fvalue(y[i])
        i = i + 1
    return 1
predictfunc = lcmaes.surrfunc_pbf.from_callable(nsurrfunc)

# pass the function and parameter to ACM surrogate cmaes, run optimization with 200 training points and collect solution object.
cmasols = lcmaes.surrpcmaes(objfunc,trainfunc,predictfunc,p,True,200)

# collect and inspect results
bcand = cmasols.best_candidate()
bx = lcmaes.get_candidate_x(bcand)
print("best x=",bx," / fevals=",cmasols.fevals())
print("distribution mean=",lcmaes.get_solution_xmean(cmasols))
cov = lcmaes.get_solution_cov(cmasols) # numpy array
#print "cov=",cov
print("elapsed time=",cmasols.elapsed_time(),"ms")
