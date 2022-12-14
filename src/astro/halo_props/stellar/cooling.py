# Put Universe MAchine Relations 
from dataclasses import dataclass, field, InitVar
# from typing import ClassVar
 

# import numpy as np
# import numpy.typing as npt
# import astropy.units as u
# import astropy.cosmology.units as cu
# import astropy.cosmology

# Type Hinting Help for Debugging 
# from typing import Type, List, Optional, Tuple, Union, Dict
#
# import ytree

# from cosmological import Cosmology 


# Cooling Rate Considerations that will be passed to halo class 
@dataclass
class CoolingModel: 
    rate: float 
    
    
    # Get Cooling Function 
    def getRate(self, T: float) -> float:
        pass 

    # Read in cooling function from data 
    @classmethod
    def fromData(cls, filename : str):
        pass 