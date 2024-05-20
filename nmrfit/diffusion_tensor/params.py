import numpy as np
from lmfit import Parameters


def make_params_difftens(model, dataset):
    params = Parameters()   
    # params.add('nfield', value=2, vary=False)
    # params.add('field', value=600, vary=False)
    params.add('model', value=int(model), vary=False)
    ## Modelfree
    ## 1 oblate        4 + Rex  7 + S2  10 + Rex + S2 
    ## 2 prolate       5 + Rex  8 + S2  11 + Rex + S2
    ## 3 anisotropic   6 + Rex  9 + S2  12 + Rex + S2

    params.add('S2', value=1, vary=False)
    params.add('ti', value=7.14e-09, vary=False)
    params.add('theta_ddcsa', value=22, vary=False)
    # if theta=="Yes":
    #     params.add('theta_ddcsa', value=22, min=0, max=45, vary=True)

    params.add('Dxx', min=1e5, value=1e7, max=1e9, vary=True)
    params.add('Dyy', min=1e5, value=1e7, max=1e9, vary=False)
    params.add('Dzz_var', min=0, value=1e4, max=1e8, vary=True)
    params.add('Dzz', expr='Dxx + Dzz_var', vary=False)
    
    params.add('Rex', value=0, vary=False)

    params.add('phi', value=0.5, min=1e-3, max=2*np.pi, vary=True)
    params.add('theta', value=0.5, min=1e-3, max=2*np.pi, vary=True)
    params.add('psi', value=0.5, min=1e-3, max=2*np.pi, vary=False)
    
    # if params['model'].value in [1,4,7,10]:
        
    
    if params['model'].value in [3,6,9,12]:
       params.add('Dyy_var', min=0, value=1e4, max=1e8, vary=True)
       params.add('Dyy', expr='Dxx + Dyy_var', vary=False)
       params.add('Dzz', expr='Dyy + Dzz_var', vary=False)
       params.add('psi', value=0.5, min=1e-3, max=2*np.pi, vary=True)
    #    params.add('Diso2', expr='((Dxx+Dyy+Dzz)/3)**2', vary=False)
    #    params.add('L2', expr='(((Dxx*Dyy)+(Dxx*Dzz)+(Dyy+Dzz))/3)**2', vary=False)
       

    if params["model"].value in [4,5,6,7,8,9,10,11,12]:
        for i in range(len(dataset)):
            params.add(f'Rex_{i}',value=3/(600**2), min=0, max=5.5e-5, vary=True)
        
    if params["model"].value in [7,8,9,10,11,12]:
        for i in range(len(dataset)):
            params.add(f'S2_{i}', value=0.8, min=0.15, max=1, vary=True)
            params.add(f'ti_{i}', value=7.14e-09, min=1e-13, max=1e-7, vary=True)
        
    return params


def make_params_difftens_intmol(model, dataset, params_difftens, theta):
    params = Parameters()   
    # params.add('nfield', value=2, vary=False)
    # params.add('field', value=600, vary=False)
    params.add('model', value=int(model), vary=False)
    ## Modelfree
    ## 1 7 + S2  10 + Rex + S2 
    ## 2 8 + S2  11 + Rex + S2
    ## 3 9 + S2  12 + Rex + S2

    params.add('S2', value=0.6, min=0.3, max=1, vary=True)
    params.add('ti', value=5e-09, min=1e-15, max=1e-7, vary=True)
    params.add('Rex',value=3/(600**2), min=0, max=5.5e-5, vary=True)
    
    params.add('theta_ddcsa', value=22, vary=False)
    if theta=="Yes":
        params.add('theta_ddcsa', value=22, min=15, max=30, vary=True)

    params.add('Dxx', value=params_difftens["Dxx"].value, vary=False)
    params.add('Dyy', value=params_difftens["Dyy"].value, vary=False)
    params.add('Dzz', value=params_difftens["Dzz"].value, vary=False)

    params.add('phi', value=params_difftens["phi"].value, vary=False)
    params.add('theta', value=params_difftens["theta"].value, vary=False)
    params.add('psi', value=params_difftens["psi"].value, vary=False)
    
    # if params['model'].value in [1,4,7,10]:
        
    
    # if params['model'].value in [6,9,12]:
    #    params.add('Dyy_var', min=0, value=1e4, max=1e8, vary=True)
    #    params.add('Dyy', expr='Dxx + Dyy_var', vary=False)
    #    params.add('Dzz', expr='Dyy + Dzz_var', vary=False)
    #    params.add('psi', value=0.5, min=1e-3, max=2*np.pi, vary=True)
    #    params.add('Diso2', expr='((Dxx+Dyy+Dzz)/3)**2', vary=False)
    #    params.add('L2', expr='(((Dxx*Dyy)+(Dxx*Dzz)+(Dyy+Dzz))/3)**2', vary=False)
       

    # if params["model"].value in [7,8,9,10,11,12]:
    #     for i in range(len(dataset)):
    #         params.add(f'Rex_{i}',value=3/(np.mean(fields)**2), min=0, max=5.5e-5, vary=True)
        
    # if params["model"].value in [7,8,9,10,11,12]:
    #     for i in range(len(dataset)):
    #         params.add(f'S2_{i}', value=0.8, min=0.15, max=1, vary=True)
    #         params.add(f'ti_{i}', value=7.14e-09, min=1e-13, max=1e-7, vary=True)
        
    return params