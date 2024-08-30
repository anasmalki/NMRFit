from lmfit import Parameters
from numpy.random import rand, uniform

### Models and params
models = [1,2,3,4,5,6,7,8,9,10,11,12]
nmodel = len(models)
name_models = ["Simplified", "Simplifed + Rex",
                "Classic Lipari-Szabo", "Classic Lipari-Szabo + Rex",
                "Extended model free", "Extended model free + Rex",
                "Extended model free 2", "Extended model free 2 + Rex",
                "IDP 2 modes", "IDP 2 modes + Rex", "IDP 3 modes" "IDP 3 modes + Rex"]

def make_params(model, tslow, tfast, theta):
    params = Parameters()
    # params.add('nfield', value=2, vary=False)
    # params.add('field', value=600, vary=False)
    params.add('model', value=int(model), vary=False)
    ## Modelfree
    # 1: simplified with ti<20ps too fast to probe
    # 2: Simplified (1) + Rex
    # 3: ti relaxaiton active, classic Lipari-Szabo
    # 4: Lipari-Szabo (2) + Rex
    # 5: Extended Modelfree, with a very fast and a slower internal motion with different S²
    # 6: Extended Modelfree+ Rex
    # 7: Extended v2
    # 8: Extended v2 + Rex
    # 9 2 modes
    # 10 2 modes + Rex
    # 11 3 modes
    # 12 3 modes + Rex

    if isinstance(tslow, (int, float)) is True:
        params.add('tm', value=tslow, vary=False)
    else:
        if tslow=="Yes":
            params.add('tm', value=8e-9, min=1e-9, max=35e-9, vary=True)
        else:
            params.add('tm', value=8e-9, vary=False)
            
  
    if isinstance(theta, (int, float)) is True:
        params.add('theta', value=theta, vary=False)
    else:
        if theta=="Yes":
            params.add('theta', value=22.5, vary=True)
        else:
            params.add('theta', value=22.5, vary=False)
    
    params.add('ts', value=1e-8, vary=False)
    # params.add('S2s', value=0.8, min=0, max=1, vary=True)

    params.add('tf', value=0, vary=False)
    params.add('S2f', value=0, vary=False)

    params.add('S3', value=0, vary=False)
    params.add('Rex', value=0, vary=False)


    if params["model"].value in [3,4,5,6,7,8] :
        params.add('tf', value=1e-10, min=1e-13, max= 5e-7, vary=True)
        params.add('S2f', value=0.4, min=0.1, max=1, vary=True)

    if params["model"].value in [5,6,7,8] :
        params.add('ts', value=1e-8, min=1e-13, max= 5e-8, vary=True)
        
    if params['model'].value in [9,10]:
        # params.add('tm', value=1e-9, min=1e-10, max=2e-9, vary=True) # mode slow
        # params.add('tm_minus_tf', value=6e-9, min=1e-10, max=1e-8, vary=True) #slow      
        # params.add('tf', value=1e-8, min=2e-9, max= 1e-6,vary=True)
        # params.add('tf', expr='tm + tm_minus_tf', vary=False)
        params.add('tm', value=1e-9, min=1.5e-9, max=35e-9, vary=True) # mode slow
        params.add('S2s', value=0.2, min=0, max=1, vary=True)#slow
        # params.add('tm_minus_tf', value=5e-9, min=1e-10, max=1e-8,vary=True) #slow        
        params.add('tf', value=1e-9, min=1e-12, max= 1.5e-9,vary=True) #median
        # params.add('tm', expr='tf + tm_minus_tf', min=1e-9, max= 35e-9, vary=False) #slow
        # params.add('Amplitude_sum', value=1, vary=False)
        params.add('S2f', expr="1 - S2s", min=0, max=1) #Fast

    if params["model"].value in [11,12]:
        # params.add('S2s', value=0.3, min=0, max=1)#slow
        params.add('tm_minus_tf', value=5e-9, min=1e-10, max=1e-8,vary=True) #slow        
        params.add('tf', value=1e-9, min=1e-11, max= 3e-9,vary=True) #median
        params.add('tm', expr='tf + tm_minus_tf', min=1e-9, max= 35e-9, vary=False) #slow
        # params.add('tm', value=10e-9, min=1e-10, max= 30e-9, vary=True) #slow
        # params.add('tf', expr='tm + tf_minus_tm', max=3e-9, vary=False) #median
        # params.add('S2f', value=0.4, min=0, max=1) #median
        # params.add('S2f_minus_S3', value=0.15, min=0, max=1, vary=True)
        # params.add('S2f', expr='S3 + S2f_minus_S3', min=0, max=1, vary=False) #median                
        # params.add('S3', value=0.3, min=0, max= 1) #fast
        # params.add('tf_minus_ts', value=50e-9, min=1e-12, max=1e-9,vary=True) #slow        
        # # params.add('ts', value=50e-12, min=1e-13, max= 1e-9, vary=False) #fast
        # params.add('ts', expr='tf + tf_minus_ts', min=1e-10, max= 30e-9, vary=False) #fast
        # params.add('Amplitude_check', expr='S2s + S2f', min=0, max=1)        
        # params.add('S3', expr='1 - S2s - S2f', min = 0, max=1) #fast
        # params.add('Amplitude_sum', expr='S2s + S2f + S3 - 1 and S3 == (1 - S2s - S2f)', value=0)
        # params.add('Zero', expr='1 - Amplitude_sum', min=0, max=1e-9
        # params.add('S3', expr='1 - (S2s + S2f)', min = 0, max=1) #fast
        # params.add('Amplitude_sum', expr='S2s + S2f + S3 if S3 != 0 or S2s !=1 else 0', value=1, vary=False)
        # params.add('S23', expr="1 - S2s", min=0, max=1) #Fast
        # params.add('Zero', expr='S2s - S23', value=0, vary=False)
        ###""
        # params.add('S2s', value=0.3, min=0, max=1)#slow
        # params.add('S3', value=0.4, min=0, max=1, vary=True) #median 
        # params.add('S2f', expr='1 - S2s - S3', min=0, max=1) #median                              
        # params.add('Amplitude_sum', expr='1 - (S2s + S2f + S3)', value=0)
        ####
        params.add('S2f', value=0.3, min=0, max=1, vary=True) #median      
        params.add('S3', value=0.4, min=0, max=1, vary=True) # fast
        params.add('S2s', expr='1 - S2f - S3', min=0, max=1)            
        params.add('Amplitude_sum', expr='(S2f + S3 + S2s) - 1', value=0, vary=False)



        if isinstance(tfast, (int, float)) is True:
            params.add('ts', value=tfast, vary=False)
        else:
            if tfast=="Yes":
                params.add('tf_minus_ts', value=50e-12, min=1e-12, max=2e-9, vary=True) #slow      
                params.add('ts', expr='tf - tf_minus_ts', min=10e-12, max=300e-12, vary=False)
            else:
                params.add('ts', value=50e-12, vary=False)


    if params["model"].value in [2,4,6,8,10,12]:
        params.add('Rex',value=3/(600**2), min=0, max=3.3e-5, vary=True)

    return params

def make_params_random(model, tslow, tfast, theta):
    params = Parameters()
    # params.add('nfield', value=2, vary=False)
    # params.add('field', value=600, vary=False)
    params.add('model', value=int(model), vary=False)
    ## Modelfree
    # 1: simplified with ti<20ps too fast to probe
    # 2: Simplified (1) + Rex
    # 3: ti relaxaiton active, classic Lipari-Szabo
    # 4: Lipari-Szabo (2) + Rex
    # 5: Extended Modelfree, with a very fast and a slower internal motion with different S²
    # 6: Extended Modelfree+ Rex
    # 7: Extended v2
    # 8: Extended v2 + Rex
    # 9 2 modes
    # 10 2 modes + Rex
    # 11 3 modes
    # 12 3 modes + Rex
    
    if isinstance(tslow, (int, float)) is True:
        params.add('tm', value=tslow, vary=False)
    else:
        if tslow=="Yes":
            params.add('tm', value=(rand()*1e-9), min=1e-9, max=35e-9, vary=True)
        else:
            params.add('tm', value=8e-9, vary=False)
            

    if isinstance(theta, (int, float)) is True:
        params.add('theta', value=theta, min=15, max=30, vary=False)
    else:
        if theta=="Yes":
            params.add('theta', value=uniform(low=21.8, high=23), min=21.8, max=23, vary=True)
        else:
            params.add('theta', value=22.5, vary=False)
            
    params.add('tm', value=(rand()*1e-9), min=1e-9, max=50e-9, vary=False) #tc

    params.add('ts', value=(rand()*1e-8), vary=False)
    # params.add('S2s', value=rand(), min=0, max=1, vary=True)

    params.add('tf', value=0, vary=False)
    params.add('S2f', value=0, vary=False)

    params.add('S3', value=0, vary=False)
    params.add('Rex', value=0, vary=False)


    if params["model"].value in [3,4,5,6,7,8] :
        params.add('tf', value=(rand()*1e-10), min=1e-13, max= 5e-7, vary=True)
        params.add('S2f', value=(rand()*1e-9), min=0.1, max=1, vary=True)

    if params["model"].value in [5,6,7,8] :
        params.add('ts', value=(rand()*1e-8), min=1e-13, max= 5e-8, vary=True)
        
    if params['model'].value in [9,10]:
        # params.add('tm', value=(rand()*1e-9), min=1e-10, max=2e-9, vary=True) # mode slow
        # params.add('tf', value=(rand()*1e-8), min=2e-9, max= 1e-6,vary=True)
        # params.add('tf', expr='ts + tf_var', vary=False)
        params.add('tm', value=(rand()*1e-9), min=1.5e-9, max=35e-9, vary=True)
        params.add('S2s', value=0.8, min=0, max=1)#slow
        # params.add('tm_minus_tf', value=(rand()*1e-9), min=1e-10, max=1e-8,vary=True) #slow    
        params.add('tf', value=(rand()*1e-8), min=1e-11, max= 1.5e-9,vary=True) #median    
        # params.add('tm', expr='tf + tm_minus_tf', min=1e-9, max= 35e-9, vary=False) #slow
        # params.add('Amplitude_sum', value=1, vary=False)
        params.add('S2f', expr="1 - S2s", min=0, max=1) #Fast

    if params["model"].value in [11,12]:
        params.add('tm_minus_tf', value=(rand()*1e-9), min=1e-10, max=1e-8,vary=True) #slow 
        params.add('tf', value=(rand()*1e-8), min=1e-11, max= 3e-9,vary=True) #median       
        params.add('tm', expr='tf + tm_minus_tf', min=1e-10, max= 35e-9, vary=False) #slow
        # params.add('tm', value=10e-9, min=1e-10, max= 30e-9, vary=True) #slow
        # params.add('tf', expr='tm + tf_minus_tm', max=3e-9, vary=False) #median
        # params.add('S2f', value=rand(), min=0, max=1) #median
        # params.add('S3', value=rand(), expr='1.000000000 - S2s - S2f', min=0, max=1) #fast
        # params.add('S2f_minus_S3', value=0.15, min=0, max=1, vary=True)
        # params.add('S2f', expr='S3 + S2f_minus_S3', min=0, max=1, vary=False) #median                
        # params.add('S3', value=0.3, min=0, max= 1) #fast
        # params.add('tf_minus_ts', value=50e-9, min=1e-12, max=1e-9,vary=True) #slow        
        # # params.add('ts', value=50e-12, min=1e-13, max= 1e-9, vary=False) #fast
        # params.add('ts', expr='tf + tf_minus_ts', min=1e-10, max= 30e-9, vary=False) #fast
        # params.add('Amplitude_check', expr='S2s + S2f', min=0, max=1)        
        # params.add('S3', expr='1 - S2s - S2f and S2s!=1', min = 0, max=1) #fast
        # params.add('Amplitude_sum', expr='S2s + S2f + S3 if S3 != 0 or S2s !=1 else 0', value=1, vary=False)
        # params.add('Amplitude_sum', expr='S2s + S2f + S3 - 1 and S3 == (1 - S2s - S2f)', value=0)
        # params.add('Zero', expr='1 - Amplitude_sum', min=0, max=1e-9)
        # params.add('S23', expr="1 - S2s", min=0, max=1) #Fast
        # params.add('Zero', expr='S2s - S23', value=0, vary=False)
        # #####
        # params.add('S2s', value=rand(), min=0, max=1)#slow
        # params.add('S3', value=rand(), min=0, max=1, vary=True) #median 
        # params.add('S2f', expr='1 - S3 - S3', min=0, max=1) #median                    
        # params.add('Amplitude_sum', expr='S2s + S2f + S3', min=0.999999, max=1.0000001)
        ####                 
        params.add('S2f', value=rand(), min=0, max=1, vary=True) #median      
        params.add('S3', value=rand(), min=0, max=1, vary=True) # fast
        params.add('S2s', value=rand(), expr='1 - S2f - S3', min=0, max=1)     
        params.add('Amplitude_sum', expr='(S2f + S3 + S2s) - 1', value=0, vary=False)
        
        if isinstance(tfast, (int, float)) is True:
            params.add('ts', value=tfast, vary=False)
        else:
            if tfast=="Yes":
                params.add('tf_minus_ts', value=(rand()*1e-12), min=1e-12, max=2e-9, vary=True) #slow      
                params.add('ts', expr='tf - tf_minus_ts', min=10e-12, max=300e-12, vary=False)
            else:
                params.add('ts', value=50e-12, vary=False)

    if params["model"].value in [2,4,6,8,10,12]:
        params.add('Rex',value=((rand()*10)/(600**2)), min=0, max=3.3e-5, vary=True)

    return params


# ### MC params function
# def params_MC(nmodel, res_params_outputs, nb_MC):
#     rc = []
#     for k in range(nmodel):
#         bb = []
#         for p in range(len(res_params[0]["num"])):
#             ff = []
#             for r in range(nb_MC):
#                 tm = np.random.normal(res_params[k].iloc[p]["tm"], res_error[k].iloc[p]["E_tm"])
#                 s2s = np.random.normal(res_params[k].iloc[p]["S2s"], res_error[k].iloc[p]["E_S2s"])
#                 if "tf" in res_params[k].columns:
#                     tf = np.random.normal(res_params[k].iloc[p]["tf"], res_error[k].iloc[p]["E_tf"])
#                 else:
#                     tf = 0
#                 if "ts" in res_params[k].columns:
#                     ts = np.random.normal(res_params[k].iloc[p]["ts"], res_error[k].iloc[p]["E_ts"])
#                 else:
#                     ts = 50e-12
#                 if "S2f" in res_params[k].columns:
#                     s2f = np.random.normal(res_params[k].iloc[p]["S2f"], res_error[k].iloc[p]["E_S2f"])
#                 else:
#                     s2f = 0
#                 if "S3" in res_params[k].columns:
#                     s3 = np.random.normal(res_params[k].iloc[p]["S3"], res_error[k].iloc[p]["E_S3"])
#                 else:
#                     s3 = 0
#                 if "Rex" in res_params[k].columns:
#                     rex = np.random.normal(res_params[k].iloc[p]["Rex"], res_error[k].iloc[p]["E_Rex"])
#                 else:
#                     rex = 0
#                 ff.append([models[k], [tm, ts, tf], [s2s, s2f, s3], rex])
#             bb.append(ff)
#         rc.append(bb)
#     return rc
