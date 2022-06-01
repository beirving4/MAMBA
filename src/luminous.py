# Put obersvational Fits in here 
import warnings
warnings.filterwarnings("ignore")

# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb, sys, os
import numpy as np
import numpy.typing as npt 
from copy import copy


# from cosmological import Cosmology

from dataclasses import dataclass, field  

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union, ClassVar 

@dataclass 
class Luminosities: 
    bolometricLuminosity: float 

    @property
    def xray(self) -> float: 
        raise NotImplementedError 

    @property 
    def uv(self) -> float: 
        raise NotImplementedError

    @property 
    def radio(self) -> float: 
        raise NotImplementedError 

    @property 
    def midIR(self) -> float:
        raise NotImplementedError

    @property 
    def nearIR(self) -> float:
        raise NotImplementedError 

    @property 
    def gamma(self) -> float: 
        raise NotImplementedError 


# Radio, IR, Mid-IR, Near-IR, 
# Optical (all the color magnitude bands)
# Soft X-ray, Hard X-ray, and Gamma Ray 

@dataclass
class QuasarLuminosities(Luminosities): 
    
    @property 
    def xray(self) -> float:
        # kBol = L/LX = Bolometric Correction Factor 
        kBol = 10.83 * (self.bolometricLuminosity/(1.0E10))**0.28 + 6.08 * (self.bolometricLuminosity/(1.0E10))**(-0.02)
        # LX = L/kBol 
        return self.bolometricLuminosity/kBol 

    @property 
    def uv(self, fuv: float=0.1) -> float: 
        return 0.5 * fuv * self.bolometricLuminosity 


# May need to put hard & Soft X-ray 
@dataclass
class GalaxyLuminosities(Luminosities):
    pass 