import unittest

import numpy as np 
import astropy.units as u
from colossus.cosmology import cosmology 
from src.constants import h_Hubble, Omega_DE, Omega_m
from cosmological import getCosmology

# separate src files from test files 

class TestCosmology(unittest.TestCase):

    # Test Factory Function 
    def test_getCosmology(self): 
        # Try with importing from parameter set 
        param_cosmo = getCosmology(parameters=(h_Hubble, Omega_m, Omega_DE))

        self.assertEqual(param_cosmo.h0, h_Hubble)
        self.assertEqual(param_cosmo.matter, Omega_m)
        self.assertEqual(param_cosmo.darkEng, Omega_DE)
        self.assertEqual(param_cosmo.curv, 0.0)
        self.assertEqual(param_cosmo.hubble, 100.0 * h_Hubble)

        cosmo_repr = lambda h0, Om, Ode: f"Cosmology(h0 = {h0:.3f}, Om0 = {Om:.3f}, Ode0 = {Ode:.3f})"

        self.assertEqual(repr(param_cosmo), cosmo_repr(h_Hubble, Omega_m, Omega_DE))

        # Test the Generation of premade cosmologies 
        for key in cosmology.cosmologies: 

            if (key == "powerlaw"):
                self.assertRaises(Exception, getCosmology, "powerlaw")
                continue 

            defined_cosmo = getCosmology(model=key)
            actual_cosmo = cosmology.setCosmology(key)

            self.assertEqual(defined_cosmo.h0, cosmology.cosmologies[key]["H0"]/100.0)
            self.assertEqual(defined_cosmo.matter, cosmology.cosmologies[key]["Om0"])
            self.assertEqual(defined_cosmo.baryon, cosmology.cosmologies[key]["Ob0"])
            self.assertEqual(defined_cosmo.darkEng, actual_cosmo.Ode0)

            # curve = 1.0 - cosmology.cosmologies[key]["Om0"] - cosmology.cosmologies[key]["Ode0"]
            self.assertEqual(defined_cosmo.curv, actual_cosmo.Ok0)
            self.assertEqual(defined_cosmo.hubble, cosmology.cosmologies[key]["H0"])

        # Test the Error Messsage is Raised for Defined Cosmology that Doesn't Exist 
        self.assertRaises(KeyError, getCosmology, model="Planck22")

        # Test for Loading Cosmology from Merger Tree Node !!!

    # Test Redshift/Time Conversions 
    def test_redshift_time(self): 

        z_to_t = lambda cosmo, z: (cosmo.age(z) * u.Gyr).to_value("Myr")
        z_to_lookback = lambda cosmo, z: (cosmo.lookbackTime(z) * u.Gyr).to_value("Myr")
        t_to_z = lambda cosmo, t: cosmo.age((t * u.Myr).to_value("Gyr"), inverse=True)
        lookback_to_z = lambda cosmo, t: cosmo.lookbackTime((t * u.Myr).to_value("Gyr"), inverse=True)

        for key in cosmology.cosmologies: 

            if (key == "powerlaw"):
                self.assertRaises(Exception, getCosmology, "powerlaw")
                continue 

            defined_cosmo = getCosmology(model=key)
            actual_cosmo = cosmology.setCosmology(key)

            for z in range(20): 
                actual = z_to_t(actual_cosmo, float(z))
                defined = defined_cosmo.redshift_to_time(float(z))
                self.assertEqual(defined, actual) 

                actual_lookback = z_to_lookback(actual_cosmo, float(z))
                defined_lookback = defined_cosmo.redshift_to_lookback_time(float(z))
                self.assertEqual(defined_lookback, actual_lookback) 

            for t in np.logspace(2, 3.5): 
                actual = t_to_z(actual_cosmo, t)
                defined = defined_cosmo.time_to_redshift(t)
                self.assertEqual(defined, actual) 

                actual_lookback = lookback_to_z(actual_cosmo, t)
                defined_lookback = defined_cosmo.lookback_time_to_redshift(t)
                self.assertEqual(defined_lookback, actual_lookback) 

# Update the test for this ... 
if __name__ == "__main__": 
    unittest.main() 