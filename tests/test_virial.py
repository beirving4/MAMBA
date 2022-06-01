import unittest

import numpy as np 
import astropy.units as u
from colossus.cosmology import cosmology 
from constants import h_Hubble, Omega_DE, Omega_m
from cosmological import getCosmology

class TestVirial(unittest.TestCase): 
    
    def test_initialize(self): 
        pass 