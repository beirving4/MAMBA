# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb, sys, os
import numpy as np
import numpy.typing as npt 
from copy import copy

from typing import List

from dataclasses import dataclass, field

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
