# May need to set virial parameters here to deal with change cosmology parameters 
from dataclasses import dataclass, field, InitVar
from typing import ClassVar
from constants import *

import numpy as np
import astropy.units as u
import astropy.cosmology.units as cu
import astropy.cosmology

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union 

import ytree

from cosmological import Cosmology 

# You want to mimimze the amount of time that this, Cosmologly, and UMRelations are created
# Put the creation in the population class and make sure that AstroObjects just call on these 
# objects to generate the relevant relations. 

# Don't use if node loaded - set virial to Optional[VirialProperties]
@dataclass(repr=False)
class VirialProperties: 
    # Set Cosmologty As Class Variable 
    cosmo : ClassVar[Cosmology] 

    # Instance Variables 
    M : InitVar[float]
    z : InitVar[float]

    # Gas composition for mu in virial temperature 
    # 0.59 = primordial, 0.61 =  ..., 1.22 = ...
    # set default to primodial 
    composition : float = field(default=0.59) 

    mass : float = field(init=False)
    radius : float = field(init=False)
    v_circ : float = field(init=False)
    dispersion : float = field(init=False)
    v_esc : float = field(init=False)
    temperature : float = field(init=False)


    # Generate Properties 
    def __post_init__(self, M : float, z : float) -> None: 
        # Virial Mass 
        self.mass = M / self.cosmo.h0 

        # Virial Radius 
        Om0, Ode0, Ok = self.cosmo.matter, self.cosmo.darkEng, self.cosmo.curve
        Omega_mz = (Om0 * (1 + z)**3.) / (Om0 * (1 + z)**3. + Ode0 + Ok * (1 + z)**2.)
        d = Omega_mz - 1 
        Delta_c = 18.0 * np.pi**2 + 82.0 * d - 39.0 * d**2 
        self.radius = 0.784 * ((self.mass/1.0E8) * (Omega_mz*18.0*np.pi**2) / (Om0 * Delta_c))**(1./3.) * (10./(1.+z))

        # Circular Velocity 
        self.v_circ = np.sqrt((G_Newton * (1.0E-6 * M_solar/kiloparsec)) * self.mass/self.radius) 

        # Dark Matter Velocity Dispersion 
        self.dispersion = self.v_circ / np.sqrt(2.0)

        # Halo Escape Velocity
        self.v_esc = 2.0 * self.v_circ 

        # Virial Temperature
        self.temperature = (self.composition * proton_mass)/(2.0 * kB) * (1000. * self.v_circ)**2

    
    def __repr__(self) -> str:
        return f"Virial Properties for {self.cosmo} and mu = {self.composition}"


# Parameters for UniverseMachine & Trinity Fits with proper central and satellite splits 
if __name__ == "__main__":
    testCosmo = Cosmology(Omega_m, Omega_DE, h0=h_Hubble)
    print(testCosmo)
    VirialProperties.cosmo = testCosmo
    print(VirialProperties)
    print(VirialProperties.cosmo)
    testVirial = VirialProperties(1.0E10, 2.0)
    # testVirial = VirialProperties(testCosmo, 1.0E10, 2.0)
    print(testVirial)
    print(f"{testVirial.mass = :.3e} Msolar/h")
    print(f"{testVirial.radius = :.3f} kpc/h")
    print(f"{testVirial.v_circ = :.3f} km/s")
    print(f"{testVirial.dispersion = :.3f} km/s")
    print(f"{testVirial.v_esc = :.3f} km/s")
    print(f"{testVirial.temperature = :.3e} K")
