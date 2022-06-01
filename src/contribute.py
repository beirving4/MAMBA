# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb, sys, os
import numpy as np
import numpy.typing as npt 
from copy import copy

from dataclasses import dataclass, field

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union, ClassVar 

from objects import AstroObject
# from black_hole import BlackHole 
# from halo import Halo 

# Mass Contributions for AstroObjects 
@dataclass
class Contributions: 
    # Associated Mass Growth Contributions 
    fromMerger : float = field(default=0.0)
    fromMergerTriggeredAccretion : float=field(default=0.0)
    fromSteadyModeAccretion : float=field(default=0.0)
    fromAccretion : float=field(default=0.0)
    fromMergerEffects : float=field(default=0.0)
    totalGrowth : float=field(default=0.0)
    cumulatives : dict = field(default_factory=dict)

    @property
    def cumulativeContributions(self) -> npt.NDArray[float]: 
        """
        Order: 
            1. Mass Contributions from Merger with same AstroObject Type
            2. Mass Contributions from Merger-Triggered Accretion 
            3. Mass Contributions from Steady-Mode Accretion
            4. Total Accretion Mass Contributions 
            5. Total Merger Effects Contributions 

            --> Normalized by the Total Growth 
        """
        
        contributions = np.array([self.fromMerger, self.fromMergerTriggeredAccretion, 
                                  self.fromSteadyModeAccretion, self.fromAccretion,
                                  self.fromMergerEffects])/self.totalGrowth 


        return contributions 

    def updateContributions(self, massMerge: Optional[float]=None, massTriggered: Optional[float]=None,
                            massSteady: Optional[float]=None) -> None: 
        
        # Individual Contributions 
        self.fromMerger += massMerge or 0.0
        self.fromMergerTriggeredAccretion += massTriggered or 0.0 
        self.fromSteadyModeAccretion += massSteady or 0.0
        
        # Update Accretion Total 
        massAccretion = (massTriggered or 0.0) + (massSteady or 0.0)
        self.fromAccretion += massAccretion
    
        # Update Merger Total 
        mergerAccretion = (massMerge or 0.0) + (massTriggered or 0.0)
        self.fromMergerEffects += mergerAccretion 

        # Update Total Growth
        self.totalGrowth = (self.fromMerger + self.fromMergerTriggeredAccretion + self.fromSteadyModeAccretion)

        # Update Contributions (ith value looks at the contributions between i-1 -> i snapshot transition)
        if self.redshift not in self.cumulatives: 
            self.cumulatives.setdefault(self.redshift, self.cumulativeContributions) 
        else: 
            self.cumulatives[self.redshift] = self.cumulativeContributions 

    # Transfer Mass Contributions 
    def transferContributions(self, otherAstro: AstroObject) -> None: 
        self.fromMerger = copy(otherAstro.fromMerger)
        self.fromMergerTriggeredAccretion = copy(otherAstro.fromMergerTriggeredAccretion)
        self.fromSteadyModeAccretion = copy(otherAstro.fromSteadyModeAccretion)
        self.fromAccretion = copy(otherAstro.fromAccretion)
        self.fromMergerEffects = copy(otherAstro.fromMergerEffects) 
        self.cumulatives = copy(otherAstro.cumulatives) 

# Figure out how to specify subclass 
# @dataclass
# class BlackHoleContributions(Contributions):

#     def transferContributions(self, otherBH : BlackHole) -> None: 
#         super().transferContributions(otherBH)

# @dataclass
# class HaloContributions(Contributions):

#     def transferContributions(self, otherHalo : Halo) -> None: 
#         super().transferContributions(otherHalo)