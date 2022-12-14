from dataclasses import dataclass, field, InitVar

import numpy as np
import astropy.constants as const 
# Type Hinting Help for Debugging 
# from typing import Type, List, Optional, Tuple, Union, Dict, ClassVar

# Class to hold stellar properties for Halo 
# UMRelations and Gas Properties from Hydro Sim
# will be children classes 
@dataclass
class StellarProperties:
    starFormationRate : float
    mass : float
    bulgeMass : float

    # Instance Variables 
    Mhalo : InitVar[float]
    z : InitVar[float]   


    def dispersion(self, M_vir: float, r_vir: float, z: float) -> float: 

        oneFej = 1./0.58
        hernquistR = 1.81527 

        Re0 = 10.0**(-0.314) * self.mass**0.042 * (1.0 + (self.mass/10.0**10.537))**0.76
        gamma = max(0, (np.log10(self.mass) - 10.75)/0.85)
        f_z = (1.0 + z)**gamma
        Re = Re0 * f_z

        stellarProfile = lambda r: self.stellarMass * ((r/Re)/((r/Re) + (1./hernquistR)))**2
        dmProfile = lambda r: ((24. * (r/r_vir)**(4./9.))/(13. + 11. * (r/r_vir)**(4./9.)))**5 
        f_vir = (self.stellarMass/M_vir) * ((r_vir/Re)/((r_vir/Re) + (1.0/hernquistR)))**2
        f = f_vir * stellarProfile(Re)/dmProfile(Re) 

        # potential = G_Newton * (1.0E-6 * M_solar/kiloparsec)
        potential = ((const.GM_sun * self.mass) / (Re * const.kpc))

        # Return the Stellar Velocity Dispersion in km/s 
        return (0.389 * np.sqrt(potential * (oneFej + 0.86/f))).to("km/s").value 