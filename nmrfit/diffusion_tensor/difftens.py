# import numpy as np
# import pandas as pd
from .params import make_params_difftens, make_params_difftens_intmol
from .data_format import make_inputs_difftens, make_inputs_difftens2
from .fit import fit_difftens, fit_intmol
from .geometric import get_NH_bond_vectors


def run_fit_difftens(relaxation_input, pdb, isotopes, offset, start_residue, end_residue, excluded, fields, model, method, nb_iter):
    relax_input = []
    for k in relaxation_input:
        ri = k.loc[(k['num'] > start_residue) & (k['num'] < end_residue)]
        if excluded is not None:
            for i in excluded:
                ri = ri.drop(ri[ri["num"]==i].index)
        relax_input.append(ri)
        
    NH_bond_vectors_input = get_NH_bond_vectors(pdb, isotopes, offset, relaxation_input[0])  
    c = NH_bond_vectors_input.loc[(NH_bond_vectors_input['num'] > start_residue) & (NH_bond_vectors_input['num'] < end_residue)]

    if excluded is not None:
        for i in excluded:
            c = c.drop(c[c["num"]==i].index)
            
    ppar = make_params_difftens(model, c)
        
    if model in [1,2,3,4,5,6]:
        a, b = make_inputs_difftens(relax_input)
        m = fit_difftens(a, b, c, fields, ppar, method, nb_iter)
        
    return m

def run_fit_difftens_intmol(relaxation_input, pdb, isotopes, offset, params_difftens, start_residue, end_residue, excluded, fields, model, method, nb_iter):
    relax_input = []
    for k in relaxation_input:
        ri = k.loc[(k['num'] > start_residue) & (k['num'] < end_residue)]
        if excluded is not None:
            for i in excluded:
                ri = ri.drop(ri[ri["num"]==i].index)
        relax_input.append(ri)
    
    NH_bond_vectors_input = get_NH_bond_vectors(pdb, isotopes, offset, relaxation_input[0])
    c = NH_bond_vectors_input.loc[(NH_bond_vectors_input['num'] > start_residue) & (NH_bond_vectors_input['num'] < end_residue)]

    if excluded is not None:
        for i in excluded:
            c = c.drop(c[c["num"]==i].index)
            
    ppar = make_params_difftens_intmol(model, c, params_difftens)
        
    # if model in [1,2,3,4,5,6]:
    #     a, b = make_inputs_difftens(relax_input)
    #     m = fit_difftens(a, b, c, fields, ppar, method, nb_iter)

    # if model in [7,8,9,10,11,12]:
    a, b = make_inputs_difftens2(relax_input)
    m = fit_intmol(a, b, c, fields, ppar, method, nb_iter)
    
    return m