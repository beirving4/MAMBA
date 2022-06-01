from dataclasses import dataclass, field


import astropy.units as u
import astropy.cosmology
from colossus.cosmology import cosmology 

# Type Hinting Help for Debugging 
from typing import Type, Optional, Tuple

import ytree 

arborType = Type[ytree.data_structures.arbor.Arbor]
definedCosmoType = Type[cosmology.Cosmology]

from constants import Omega_bary


@dataclass(repr=False) 
class Cosmology: 
    h0: float 
    matter: float
    darkEng: float 
    baryon: float = field(default=Omega_bary)
    curv: float = field(default=None)
    hubble: float = field(init=False)
    rad: float = field(default=0.0) 

    _cosmo : Optional[definedCosmoType] = field(default=None)

    def __post_init__(self) -> None:
        
        # Hubble Constant (km/ (Mpc * s))
        self.hubble = self.h0 * 100.0 

        # Set Cosmology Attribute 
        if self._cosmo is None:  
            astropy_cosmo = astropy.cosmology.FlatLambdaCDM(self.hubble, self.matter, Ob0=self.baryon)
            self._cosmo = cosmology.fromAstropy(astropy_cosmo, sigma8=0.8, ns=0.97, cosmo_name="customCosmo")

        # Cosmological Curvature Density  
        if self.curv is None: self.curv = self._cosmo.Ok0

    def __repr__(self) -> str:
        return f"Cosmology(h0 = {self.h0:.3f}, Om0 = {self.matter:.3f}, Ode0 = {self.darkEng:.3f})" 

    def redshift_to_time(self, redshift : float) -> float:
        return (self._cosmo.age(redshift) * u.Gyr).to_value("Myr") 

    def redshift_to_lookback_time(self, redshift : float) -> float:
        return self._cosmo.lookbackTime(redshift) * 1000.0

    def time_to_redshift(self, time : float) -> float:
        return self._cosmo.age((time * u.Myr).to_value("Gyr"), inverse=True)

    def lookback_time_to_redshift(self, lookback_time : float) -> float:
        return self._cosmo.lookbackTime((lookback_time * u.Myr).to_value("Gyr"), inverse=True) 

    # Load from Astropy
    @classmethod
    def fromDefinedCosmologies(cls, model : str): 

        defined_cosmo = cosmology.setCosmology(model)

        return cls(
            h0 = defined_cosmo.H0/100.0,
            matter = defined_cosmo.Om0, 
            darkEng = defined_cosmo.Ode0,  
            baryon = defined_cosmo.Ob0, 
            _cosmo = defined_cosmo
        )


    @classmethod
    def fromSimulation(cls, simulation : arborType):

        return cls(
            h0 = simulation.hubble_constant,
            matter = simulation.omega_matter, 
            darkEng = simulation.omega_lambda 
        )


# Factory to Return the Correct Cosmology Type
# in Python 3.10 will redo as ...
# input : Optional[str | arborType | Tuple[float]]
# use match case in here ... 
def getCosmology(model : Optional[str]=None, simulation : Optional[arborType]=None,
                 parameters : Optional[Tuple[float]]=None, printing: bool=False) -> Type[Cosmology]:

    # Load From Defined Cosmology 
    if (model is not None):
        if (model not in cosmology.cosmologies.keys()):
            raise KeyError("Input is not a defined cosmology")
        elif (model == "powerlaw"):
            raise Exception("Power Law Cosmologies Haven't Been Implemented Yet!")  
        else: 
            cosmo = Cosmology.fromDefinedCosmologies(model)
            if printing: print(f"Loading from the {model} cosmology...")
    
    # Load From Simulation
    if simulation is not None: 
        cosmo = Cosmology.fromSimulation(simulation)
        if printing: print(f"Loading the cosmology associated with the simulation ...")


    # Load from Parameters 
    if parameters is not None: 
        assert len(parameters) == 3, "Make sure you input (h0, Om0, Ode0)!"
        h0, Om0, Ode0 = parameters
        cosmo = Cosmology(h0, Om0, Ode0) 

    if printing: print(f"Generating {cosmo}\n")

    return cosmo 