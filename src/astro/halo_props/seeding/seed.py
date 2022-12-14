# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb, sys, os
import numpy as np, ytree
import numpy.typing as npt 
from copy import copy

# from virial import VirialProperties 
from ..stellar.galactic import GalacticProperties 
# from galactic import GalacticProperties 

from dataclasses import dataclass, field, InitVar 
from typing import Type, List, Optional, Tuple, Union, ClassVar

treeNode = Type[ytree.data_structures.tree_node.TreeNode]

# Down the line make a class of light, medium, and heavy seeding models 
# Then build a factory function to set a light, medium, and heavy attributes
# In the larger SeedingModel class

# Black Hole Seeding Model 
@dataclass
class SeedingModel: 
    # virialTemp: float
    # gasTemp: float 
    # maxTemp: float 
    # spin_ratio: float 

    # virial: ClassVar[VirialProperties] = field(default=VirialProperties)
    galaxy: ClassVar[GalacticProperties] = field(default=GalacticProperties)

    # def __post_init__(self) -> None: 
        # pass

    # Seeds from Pop III remnants 
    def makeLightSeed(self) -> float: 
        return 0.0 

    # Seeds from Runaway Collisions in Nuclear Star Clusters 
    def makeMediumSeed(self) -> float:
        return 0.0 

    # Seeds from Direct Collapse 
    def makeHeavySeed(self, spin: float, T_vir) -> float: 
        # 1. Check if Gas Collapses to Form Disc in Halo
        if (T_vir <= self.galaxy.gasTemperature): return 0.0 

        # 2. Check if Disc is Gravitationally Unstable
        spin_ratio = spin/self.galaxy.spinMax(T_vir)
        if (spin_ratio > 1.0): return 0.0 

        # 3. Check to See if Fragmentation Occurs 
        if (T_vir > self.galaxy.tempMax(spin, T_vir) > 1.0): return 0.0 

        # Seeding conditions have been met, form & return heavy DCBH seed
        return self.galaxy.discMass * (1.0 - np.sqrt(spin_ratio))