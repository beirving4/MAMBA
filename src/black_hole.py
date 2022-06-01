# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb
import numpy as np
import numpy.typing as npt 
from copy import copy
from constants import *
from contribute import Contributions
# from population import *


# from halo import *
from objects import AstroObject
from accretion import AccretionModel  
from history import BlackHoleHistory 
from luminous import QuasarLuminosities 

from dataclasses import dataclass, field 

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union, ClassVar

# Black Hole Contributions 
@dataclass
class BlackHoleContributions(Contributions):

    def transferContributions(self, otherBH : BlackHole) -> None: 
        super().transferContributions(otherBH)


# Generate Class for Black Hole Phenomena 
@dataclass(order=True, repr=False) # will add slots=True after Python 3.10 upgrade 
class BlackHole(AstroObject): 
    # Basic Properties 
    redshift: float = field(default=0.0) 
    spin: float = field(default=np.random.uniform(-1.0, 1.0))
    radiative_efficiency: float = field(default=0.1)
    accretion_budget: float = field(default=0.0)
    massAccretionRatio: float = field(default=0.0)
    salpeter_time: float = field(init=False)
    eddingtonRate: float = field(init=False)

    # Associated Accretion Model 
    accretion: ClassVar[AccretionModel] = field(default=AccretionModel) 

    # For Monitoring the Mass Contributions During Evolution
    # May need to be a default_factory 
    gains: BlackHoleContributions = field(default_factory=BlackHoleContributions)

    # Seeding Statistics & History 
    history: BlackHoleHistory = field(default_factory=BlackHoleHistory)

    # Observational Luminosity 
    luminosity: QuasarLuminosities = field(default_factory=QuasarLuminosities)

    # Merger Probability
    merger_probability: float = field(default=0.1) 

    # Host Information (Halo ID to link to Halo Object)
    hostHalo : int = field(default=int)  
    haloHistory : List[int] = field(default_factory=list)

    schwarzschildRadius : float = field(default=0.0) 
    apparentRadius : float = field(default=0.0) 
    bondiRadius : float = field(default=0.0) 
    offsetDistance : float = field(default=0.0) 


    # Post-Initialization Routine 
    def __post_init__(self) -> None:
        # Populate Remaining Properties
        self.salpeter_time = t_eddington_Myr * (self.radiative_efficiency/(1.0 - self.radiative_efficiency))
        self.eddingtonRate = self.mass / (self.radiative_efficiency * t_eddington_Myr)

        # Set Seeding Statistics 
        self.history.seedMass = self.mass
        self.history.seedRedshift = self.redshift
        self.history.seedHaloID = self.hostHalo 
        
        # Initialize Tracking of Halo History 
        self.history.updateHistory(self.hostHalo)


    # Representer 
    def __repr__(self) -> str:
        return f"Black Hole(M = {self.mass:.3e} Msol, spin = {self.spin:.3f}, z = {self.redshift:.3f})"


    # Transfer Mass Conrtributions 
    def transferInfo(self, otherBH: BlackHole) -> None:
        # Transfer Mass Contribution Info 
        self.gains.transferContributions(otherBH)
        
        # Transfer Historical Information 
        self.history.transferHistory(otherBH)


    # Black Hole Growth CDFs 
    @property
    def cumulativeContributions(self) -> npt.NDArray[float]: 
        """
        Order: 
            1. Mass Contributions from Black Hole Merger 
            2. Mass Contributions from Merger-Triggered Accretion 
            3. Mass Contributions from Steady-Mode Accretion
            4. Total Accretion Mass Contributions 
            5. Total Merger Effects Contributions 

            --> Normalized by the Total Growth 
        """
        
        contributions = np.array([self.fromBHMerger, self.fromMergerTriggeredAccretion, 
                                  self.fromSteadyModeAccretion, self.fromAccretion,
                                  self.fromMergerEffects])/self.totalGrowth 


        return contributions 

    # Update Eddington Rate 
    def updateEddingtonRate(self) -> None:
        self.eddingtonRate = self.mass / (self.radiative_efficiency * t_eddington_Myr)
    
    # Eddington Mass Accretion Ratio from Models (Run at Eddington by default)
    # set Luminosites in here 
    def eddingtonRatio(self, dt: Optional[float]=None, ratio: Optional[float]=None) -> None: 
        if ((dt is None) and (ratio is None)): 
            self.massAccretionRatio = 1.0 
        elif (ratio):
            self.massAccretionRatio = ratio 
        else: 
            self.massAccretionRatio = min(1.0, np.log(1.0 + self.accretion_budget/self.mass) * self.salpeter_time/dt)
    
    # Black Hole Accretion Rate in (Msolar/Myr) 
    @property
    def accretionRate(self) -> float:
        return self.massAccretionRatio * self.eddingtonRate 

    # Bolometric Accretion Rate
    @property
    def bolometricLuminosity(self) -> float:
        return self.radiative_efficiency * self.accretionRate * (c_LIGHT**2 * M_solar)/(megaYear *L_solar) 

    # Generate Luminosities 
    def generateLuminosities(self) -> None:
        self.luminosity(self.bolometricLuminosity)


    # Mass Growth By Accretion (in solar masses)
    def accretionGrowth(self, dt: float, accretion_type: Optional[int]=None, ratio: Optional[float]=None) -> None: 
        """
        Accretion Types: 
    
            1. Merger Triggered 
            2. Steady Mode
        """

        # Only Run if AccretionBudget > 0.0 
        if (self.accretion_budget > 0.0): 

            # Set the Mass Accretion Rate 
            if ratio: 
                self.eddingtonRatio(ratio=ratio)
            else:
                self.eddingtonRatio(dt)

            # Update Mass 
            self.mass *= np.exp(self.massAccretionRatio * (dt/self.salpeter_time))

            # Deplete Budget 
            self.accretion_budget -= self.accretionRate * dt 

            # Update Contribution Lists 
            if (accretion_type == 1): self.gains.updateContributions(massTriggered=self.accretionRate * dt) 
            if (accretion_type == 2): self.gains.updateContributions(massSteady=self.accretionRate * dt) 

            # Update Eddington Rate
            self.updateEddingtonRate()  
    
    # Transition Values For Merger 
    def mergerGrowth(self, otherBH: BlackHole) -> BlackHole:
        # Add Masses 
        massBH = self.mass + otherBH.mass 
    
        # Add Spins 
        spinBH = np.random.uniform(self.spin, otherBH.spin)

        # Update Redshift 
        redshift = min(self.redshift, otherBH.redshift)

        # Find the Heavier & Lighter Partner 
        # (if the same, let the heavier partner be the one with the lighter seed mass)
        if np.isclose(self.mass, otherBH.mass): 
            lighterBH = max((self, otherBH), key=lambda bh: bh.seedMass)
            heavierBH = min((self, otherBH), key=lambda bh: bh.seedMass)
        else:
            lighterBH, heavierBH = min(self, otherBH), max(self, otherBH) 

        # Updated Black Hole with inherited Heavier BH values (maybe transfer over lighterBH contributions as well?) 
        newBH = self.__class__(massBH, redshift, spinBH) 
        newBH.transferInfo(heavierBH)
        newBH.gains.updateContributions(massMerge=lighterBH.mass)

        # Return Updated Black Hole 
        return newBH


    # Recoil Velocity for Black Hole Enounter 
    def recoilVelocity(self, otherBH: BlackHole) -> float: 

        # Spins 
        spin1, spin2 = abs(self.spin), abs(otherBH.spin)
        q, eta = self.massRatio(otherBH), self.symmetricMassRatio(otherBH)
        # alpha = np.sqrt(spin1 * spin1 + spin2 * spin2)
  
        # Fit Parameters 
        # Am, Bm
        # D, E, F = 3684.73, 0.0705, -0.6238
        # dphi, phi1 = np.radians(4678.90), np.radians(0.2447) 

        # Generate Vector Angles 
        theta = np.random.uniform(0.0, 2.0*np.pi)

        # Velocity Weights 
        # V1 = ((1 + E * alpha * np.cos(theta))/(1 + F * alpha * np.cos(theta))) * D * alpha * np.sin(theta)

        # Return the Recoil Velocity 
        # return abs(V1 * np.cos(dphi + phi1)) 

    # Compare Mass Ratio with Partner Black Hole
    def massRatio(self, otherBH: BlackHole) -> float: 
        return super().massRatio(otherBH)

    # Symmetric Mass Ratio
    def symmetricMassRatio(self, otherBH: BlackHole) -> float:
        q = self.massRatio(otherBH)
        return q/((1.0 + q)*(1.0 + q)) 

    # Track the Halo History of this Black Hole 
    def updateHaloHistory(self) -> None:
        pass 
        

    def __add__(self , otherBH: BlackHole) -> float:
        return self.mass + otherBH.mass 

    # Check Which Black Hole is Spinning Faster 
    def isFaster(self, otherBH: BlackHole) -> bool:
        return abs(self.spin) > abs(otherBH.spin) 

    # Check Which Black Hole is Spinning Faster 
    def sameSpeed(self, otherBH: BlackHole) -> bool:
        return abs(self.spin) == abs(otherBH.spin) 

    # Check Which Black Hole is Closer 
    def isCloser(self, otherBH: BlackHole) -> bool:
        return self.redshift < otherBH.redshift 

    # Check Which Black Hole is Closer 
    def sameDistance(self, otherBH: BlackHole) -> bool:
        return self.redshift == otherBH.redshift 


# 
# Factory Function to Set Class Parameters for BH 
# def setupBlackHoles() -> BlackHole: 
#     pass 


# Test Functionality 
# if __name__ == "__main__":
