# Store Annotations as String for Future Call Outs 
from __future__ import annotations

import warnings
warnings.filterwarnings("ignore")

# Import Libraries 
import pdb, sys, os
import numpy as np, ytree 
import numpy.typing as npt 
from copy import copy

# from constants import *
# from population import *

from black_hole import BlackHole 
from objects import AstroObject 
from virial import VirialProperties 
from stellar_props import StellarProperties
from cooling import CoolingModel 
from history import HaloHistory 
from galactic import GalacticProperties 
from seed import SeedingModel 
from contribute import Contributions 

from pympler import asizeof 


from dataclasses import dataclass, field 

from typing import Type, List, Optional, Tuple, Union, ClassVar

treeNode = Type[ytree.data_structures.tree_node.TreeNode]

# Halo Growth Contributions 
@dataclass
class HaloContributions(Contributions):

    def transferContributions(self, otherHalo : Halo) -> None: 
        super().transferContributions(otherHalo)

# Moving mass min, peakHalo, mass_min to Universe class 
# Have neighbors return the class object and store the link if possible
@dataclass(order=True, repr=False) # __init_subclass__ & new to load specific halo types 
class Halo(AstroObject): 
    # Basic Properties 
    spin: float = field(default=np.random.lognormal(0.035, sigma=0.5))
    time: float = field(init=False)
    scale_factor: float = field(init=False)

    # Class Variables (These may have to be default_factory)
    virial: ClassVar[VirialProperties] = field(default=VirialProperties)
    galaxy: ClassVar[GalacticProperties] = field(default=GalacticProperties)
    stellar: ClassVar[StellarProperties] = field(default=StellarProperties)
    cooling: ClassVar[CoolingModel] = field(default=CoolingModel) 
    seeding: ClassVar[SeedingModel] = field(default=SeedingModel)

    # Connection to Larger Simulation 
    peak_mass: float = field(default=None) 

    # Seeding Statistics & History 
    history: HaloHistory = field(default_factory=HaloHistory)

    # Initialize Halo w/ No Black Hole 
    BlackHole: Type[BlackHole] = field(default=None) 
    populationBH: dict = field(default_factory=dict) 

    def __post_init__(self) -> None:
        # Initialize Black Hole Population 
        self.populationBH = {"Host": self.BlackHole, "Bounded": [], "Wandering": []}

        # Set Virial Properties 
        self.virial(self.mass, self.redshift) 

        # Set Galactic Properties 
        self.galaxy(self.mass)

        # Initialize Seeding Model
        self.seeding.galaxy = self.galaxy 

        # Set Stellar Properties (Make this different for DM and Gaseuous Halo)
        # super().__post_init__(mass, z) vs. super().__post_init__(mass, z, peak_mass)

    # Mass Ratio 
    def massRatio(self, otherHalo: Halo) -> float:
        return super().massRatio(otherHalo)

    # Indicate if Major Merger Occured 
    @property 
    def isMajorMerger(self, otherHalo: Halo) -> bool: 
        return True if (self.massRatio(otherHalo) >= 0.1) else False

    # Generate Black Hole Seed
    def BlackHoleSeeding(self) -> None:
        """
        Seeding Options: 
            1. Heavy Seeds - Rare, High Mass Seeds According to (Natarajan, 2006)
            2. Medium Seed - Seeds from Runaway Collisions in NSCs
            3. Light Seeds - Look into IMF for low mass Pop III (Stacy, et. al. 2016 ??)
            ---> Also try M.Latif, 2021

        Update so that BH is initialized in the function or nothing happens at all
        """

        # Heavy Seeds 
        seedMass = self.seeding.makeHeavySeed(self.spin, self.virial.temperature)

        # Medium Seeds 
        if (seedMass == 0.0): 
            seedMass = self.seeding.makeMediumSeed() 

        # Light Seeds 
        if (seedMass == 0.0): 
            seedMass = self.seeding.makeLightSeed() 

        # Initialize Black Hole 
        self.BlackHole = BlackHole(seedMass, redshift = self.redshift)

        # Update Seeding History and Black Hole Population 
        if (seedMass != 0.0): 
            self.history.hasSeeded = True 
            self.populationBH["Host"] = self.BlackHole 

    @property
    def massCap(self) -> float: # sigmaBH = 0.285 
        betaBH = 8.47 + 0.553 * (1.0 - self.scale_factor) + 0.023 * self.redshift
        gammaBH = 1.082 - 0.221 * (1.0 - self.scale_factor) + 0.133 * self.redshift 
        logMassCap = betaBH + gammaBH * np.log10(self.bulgeMass/1.0E11)
        return 10.0**logMassCap 

    # Cooling Gas 
    def gasCooling(self) -> None:
        pass 

    # Gravitational Wave Recoil 
    def gravitationalWaveRecoil(self) -> None: 
        for otherBH in self.populationBH["Bounded"]:
            if otherBH is None:
                self.populationBH["Bounded"].remove(otherBH)
            else: 
                v_recoil = self.BlackHole.recoilVelocity(otherBH)
                
                # Eject Black Hole, and send it wandering
                if (v_recoil > self.escapeVelocity):
                    if otherBH in self.populationBH["Bounded"]: self.populationBH["Bounded"].remove(otherBH)
                    self.populationBH["Wandering"].append(otherBH)


if __name__ == "__main__": 
    pass 
    # M0vals = np.array([12.035, 12.081, 11.896, 12.021, 12.069, 12.054])
