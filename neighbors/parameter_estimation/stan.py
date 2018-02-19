"""Module containing stan models"""
import pystan
import os, pickle
from glob import glob
#from hashlib import md5

__all__ = ["get_model"]

_modeldir = os.path.dirname(__file__)

def available():
    """Return a list of available models"""
    return list(map(os.path.basename, glob(_modeldir+"/*.stan")))


#TODO: hash model code inside file to auto-recompile
# when updated
def get_model(model_name, recompile=False, **kwargs):
    """Get cached stan model or compile it"""
    modelfn = _modeldir + "/" + model_name + ".stan"
    if not os.path.exists(modelfn):
        raise ValueError("{:s} does not exist!".format(modelfn))
    cachefn = modelfn.replace("stan", "pkl")
    cached = os.path.exists(cachefn)
    if recompile or not cached:
        sm = pystan.StanModel(modelfn)
        with open(cachefn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        sm = pickle.load(open(cachefn, 'rb'))
    return sm
