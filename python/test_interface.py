import lcmaes_interface as lci

# setup input parameters                                                                              
myfun = lambda x: sum([xi**2 for xi in x])  # myfun accepts a list of numbers as input                
x0 = [2.1] * 10                                                                                       
sigma0 = 0.1                                                                                          
mlbounds = [-4]*10
mubounds = [4]*10

# run optimization via lci                                                                            
res = lci.pcmaes(lci.to_fitfunc(myfun),                                                        
                 lci.to_params(x0, sigma0,                                                     
                               str_algo=b'abipop',quiet=True, # b=bytes, unicode fails                               
                               lbounds=mlbounds,ubounds=mubounds,
                               scaling=False,
                               restarts=2))
lci.plot()  # plot from file set in lci.to_params  
