# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb, sys, os
import numpy as np
import numpy.typing as npt 
from copy import copy

from typing import List

from dataclasses import dataclass, field

from objects import AstroObject 
from black_hole import BlackHole
from halo import Halo 

# Class Object to Hold Historical Information 
@dataclass 
class ObjectHistory: 
    seedMass: float = field(default=None)
    seedRedshift: float = field(default=None)
    seedHaloID: int = field(default=None)
    seedSpin: float = field(default=None)

    haloHistory : List[int] = field(default_factory=list)

    def updateHistory(self, haloID: int) -> None:
        self.haloHistory.append(haloID)

    def transferHistory(self, otherObject: AstroObject) -> None:
        self.seedMass = copy(otherObject.seedMass)
        self.seedRedshift = copy(otherObject.seedRedshift)

@dataclass
class BlackHoleHistory(ObjectHistory): 

    def transferHistory(self, otherBH: BlackHole) -> None:
        super().transferHistory(otherBH)
        self.updateHistory(otherBH.hostHalo)

@dataclass
class HaloHistory(ObjectHistory):
    hasSeeded: bool=field(default=False)

    def transferHistory(self, otherHalo: Halo) -> None:
        super().transferHistory(otherHalo)
        self.updateHistory(otherHalo)
        self.hasSeeded = copy(otherHalo.history.hasSeeded)
        # self.seedSpin = otherHalo.history.seedSpin