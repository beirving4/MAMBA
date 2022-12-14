import numpy as np, pdb 

from typing import ClassVar
from astropy import constants as const
from dataclasses import dataclass, field, InitVar

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

        # Virial Radius [kpc]
        Om0, Ode0, Ok = self.cosmo.matter, self.cosmo.darkEng, self.cosmo.curve
        Omega_mz = (Om0 * (1 + z)**3.) / (Om0 * (1 + z)**3. + Ode0 + Ok * (1 + z)**2.)
        d = Omega_mz - 1 
        Delta_c = 18.0 * np.pi**2 + 82.0 * d - 39.0 * d**2 
        scaled_mass = (self.mass/1.0E8) * (Omega_mz*18.0*np.pi**2) / (Om0 * Delta_c)
        self.radius = 0.784 * (scaled_mass)**(1./3.) * (10./(1.+z))

        # Circular Velocity [km/s]
        v_circ = np.sqrt((const.GM_sun * self.mass) / (const.kpc * self.radius))
        self.v_circ = v_circ.to("km/s").value

        # Dark Matter Velocity Dispersion [km/s]
        self.dispersion = self.v_circ / np.sqrt(2.0)

        # Halo Escape Velocity
        self.v_esc = 2.0 * self.v_circ 

        # Virial Temperature
        factor = ((self.composition * const.m_p) / (2.0 * const.k_B)) 
        self.temperature = (factor * v_circ * v_circ).to("K").value


    
    def __repr__(self) -> str:
        return f"Virial Properties for {self.cosmo} and mu = {self.composition}"


# Parameters for UniverseMachine & Trinity Fits with proper central and satellite splits 
if __name__ == "__main__":
    testCosmo = Cosmology(0.674, 0.3158, 0.69)
    print(testCosmo)
    VirialProperties.cosmo = testCosmo
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
