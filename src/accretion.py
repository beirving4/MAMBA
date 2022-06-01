# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb
from sys import prefix
import numpy as np
import numpy.typing as npt 

from scipy.stats import rv_continuous

from steady import SteadyModeAccretion
from triggered import MergerTriggeredAccretion 

from copy import copy
from constants import *
from population import *
from halo import *

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union

from dataclasses import dataclass, field 

# Reference __new__ vs __init__ video from mCoding on how to do thi propoerly 
# Also do this for Halo Class as well 

# Use a factory pattern for this .... 
@dataclass 
class AccretionModel: 
    # The leading _ makes the element private/protected to the class 
    # _registry = {} # Ex: burst -> BurstModel 

    triggeredModel : MergerTriggeredAccretion
    steadyModel : SteadyModeAccretion 

    '''

    gasDensity : float

    # Initialize Model Type 
    def __init_subclass__(cls, prefix: str, **kwargs) -> None:
        super().__init_subclass__(**kwargs)
        cls._registry[prefix] = cls
    
    def __new__(cls: type[self], key: Optional[int]=None) -> cls: 
        pass
    '''


# Make this a class method 
# class SimulatedAccretion: 
#     pass 

# Factory Function to Get the Accretion Model for Black Holes 
# def getAccretionModel() -> AccretionModel:
#     triggeredModel = ...
#     steadyModel = ... 
#     return AccretionModel(triggeredModel, steadyModel)


# Trinity Eddington Ratio Distribution 
class TrinityEddington(rv_continuous, AccretionModel, prefix="trinity"):
    def _pdf(self, eta: float, massBH: float, z: float):
        Mc, w = 7.827, 0.739
        f20, f2a = 0.308, 0.697 
        P0, eta0 = ...  # Normaliztion constants ... 
        c10, c1a = -0.071, 0.796 
        c20, c2a = 1.398, 0.126

        x = (np.log10(massBH) - Mc)/w

        a = 1.0/(1.0 + z)
        f1 = np.exp(x)/(1.0 + np.exp(x))
        f2 = f20 + f2a * (a - 1.0)
        f_duty = f1 * f2 

        c1 = c10 + c1a * (a - 1.0)
        c2 = c20 + c2a * (a - 1.0)
        pdf = (P0 * f_duty)/((eta/eta0)**c1 + (eta/eta0)**c2) + (1.0 + f_duty) 

        # Need to multiply the Dirac delta around eta == 0 on last term 
    
        return super()._pdf(x, *args)

    def random_state(self):
        return super().random_state
    
    def _cdf(self, x, *args):
        return super()._cdf(x, *args)


# if __name__ == "__main__": 

   