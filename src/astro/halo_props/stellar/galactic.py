# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb, sys, os
import numpy as np, ytree
import numpy.typing as npt 
from copy import copy

from dataclasses import dataclass, field, InitVar 
from typing import Type, List, Optional, Tuple, Union, ClassVar

treeNode = Type[ytree.data_structures.tree_node.TreeNode]

# Class to Hold Galactic Properties 
# May put stellar and/or cooling in here 
@dataclass 
class GalacticProperties: 
    # Instance Variables 
    haloMass: InitVar[float]
    discMass: float = field(init=False)

    gasTemperature: float = field(default=5000.0)
    discMassFraction: float = field(default=0.05)

    # Critical Values for Stability Conditions 
    criticalTorque: float = field(default=0.6)
    criticalToomre: float = field(default=3.0)

    def __post_init__(self, haloMass: float) -> None:
        self.discMass = haloMass * self.discMassFraction 

    @classmethod
    def fromSimulation(cls, node: treeNode):
        return cls(float(node["mass"]))

    @property
    def spinMax(self, T_vir: float) -> float: 
        m_d, j_d = self.discMassFraction, self.discAngularMomentumFraction
        Q_c, T_gas = self.criticalToomre, self.gasTemperature
        return (m_d * Q_c)/(8.0 * (m_d/j_d))*np.sqrt(T_vir/T_gas)

    @property
    def tempMax(self, spin: float, T_vir: float) -> float:
        # Necessary Parameters
        m_d, alpha_c = self.discMassFraction, self.criticalTorque
        spin_ratio = spin/self.spinMax(T_vir)
        if (spin_ratio > 1.0):
            return self.gasTemperature * (4.0 * alpha_c/m_d)**(2.0/3)
        else:
            return self.gasTemperature * ((4*alpha_c/m_d) / (2 - np.sqrt(spin_ratio)))**(2./3)