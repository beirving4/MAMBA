# Store Annotations as String for Future Call Outs 
from __future__ import annotations

# Import Libraries 
import pdb
from sys import prefix
import numpy as np
import numpy.typing as npt 

from scipy.stats import rv_continuous

# from copy import copy
# from constants import *
# from population import *
# from halo import *

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union

from dataclasses import dataclass, field 


# Put merger triggered models in their own file 
class MergerTriggeredAccretion: 
    pass 

class BurstModel(MergerTriggeredAccretion):
    pass 