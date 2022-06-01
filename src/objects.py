# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb

from cosmological import Cosmology, getCosmology
# from contribute import Contributions 

from dataclasses import dataclass, field

# Type Hinting Help for Debugging 
from typing import Optional, ClassVar 


# Astrophysical Object Metaclass 
@dataclass
class AstroObject: 
    # Basic Properties for All Objects 
    mass : float
    redshift : float
    spin : float
    time : float = field(default=float, repr=False) 
    scale_factor : float = field(default=float, repr=False)

    # Mass Contributions 
    # gains: Contributions = field(default=None)

    # Associated Cosmological Model 
    cosmo : ClassVar[Cosmology] = field(default=Cosmology)

    # ID Parameter (Get ID)
    _id : Optional[int] = field(default=None)

    # Comparison Attribute  
    sort_index : float = field(default=float, repr=False)

    # Don't create New Object if Mass <= 0.0 
    def __new__(cls, mass : float, *args, **kwargs):
        if mass <= 0.0: return None
        return super().__new__(cls)

    # Set Remaining Key Parameters 
    def __post_init__(self) -> None: 
        # Setting Time Parameters  
        self.time = self.cosmo.redshift_to_time(self.redshift)
        self.scale_factor = 1.0/(1.0 + self.redshift) 

        # Compare by Mass -> Age (__setattr__ more robust)
        object.__setattr__(self, "sort_index", self.mass)
        object.__setattr__(self, "sort_index", self.redshift) 

    # Mass Ratio 
    def massRatio(self, otherAstro: AstroObject) -> float:
        return min(self.mass, otherAstro.mass)/max(self.mass, otherAstro.mass)

    # def getID(self) -> int:
    #     object.__getattribute__(self, "_id", )

    # def mergerGrowth(self, otherAstro: AstroObject) -> AstroObject: 
    #     NotImplemented 
    
# Factory Function to Set Class Parameters for Astro Objects 
# def setupAstroObjects() -> BlackHole: 
#     pass 






if __name__ == "__main__": 
    AstroObject.cosmo = getCosmology("Planck18") 
    goodAstro = AstroObject(1.0E6, 0.0, 0.9)
    print(f"{goodAstro = }\n")
    badAstro = AstroObject(0.0, 0.0, 0.9)
    print(f"{badAstro = }\n")