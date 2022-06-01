# Put Universe MAchine Relations 
from dataclasses import dataclass, field, InitVar
from constants import *

import numpy as np

# Type Hinting Help for Debugging 
from typing import Type, List, Optional, Tuple, Union, Dict, ClassVar

from stellar_props import StellarProperties 
from cosmological import Cosmology 

# UniverseMachine Stellar Mass Fit 
def stellarMassUM(Mhalo : float, z : float, Mpeak : float, **fit_parameters) -> float:
    # Load Stellar Mass Fit Parameters 
    M0, Ma, Mlna, Mz = fit_parameters["M"]
    eps0, epsA, epsLnA, epsZ = fit_parameters["eps"]
    alpha0, alphaA, alphaLnA, alphaZ = fit_parameters["alpha"]
    beta0, betaA, betaZ = fit_parameters["beta"]
    gamma0, gammaA, gammaZ = fit_parameters["gamma"]
    delta0 = fit_parameters["delta"]

    # Convert to Scale Factor
    a = 1.0/(1.0 + z)

    # Fit Relations 
    log10M1 = M0 + Ma * (a - 1.) - Mlna * np.log(a) + Mz * z 
    epsilon = eps0 + epsA * (a - 1) - epsLnA * np.log(a) + epsZ * z 
    alpha = alpha0 + alphaA * (a - 1.) - alphaLnA * np.log(a) + alphaZ * z
    beta = beta0 + betaA * (a - 1.) + betaZ * z
    delta = delta0 
    logGamma = gamma0 + gammaA * (a - 1.) + gammaZ * z

    x = np.log10(Mpeak) - log10M1
    gamma = 10.0**logGamma

    logM = log10M1 - np.log10(10.0**(-alpha * x) + 10.0**(-beta * x))
    logStellarMass = logM + epsilon + gamma * np.exp(-0.5 * (x/delta)**2)

    return 10.0**logStellarMass

# Universe Machine Star Formation Rate 
def sfrUM(z : float, vMpeak : float, **fit_parameters) -> float:
    # Star Formation Rate Fit Parameters 
    V0, Va, Vz1, Vz2 = fit_parameters["V"]
    eps0_sfr, epsA_sfr, epsLnA_sfr, epsZ_sfr = fit_parameters["eps"]
    alpha0_sfr, alphaA_sfr, alphaLnA_sfr, alphaZ_sfr = fit_parameters["alpha"]
    beta0_sfr, betaA_sfr, betaZ_sfr = fit_parameters["beta"]
    gamma0_sfr, gammaA_sfr, gammaZ_sfr = fit_parameters["gamma"]
    delta0_sfr = fit_parameters["delta"]

    # Convert to Scale Factor
    a = 1.0/(1.0 + z)

    # Star Formation Rate Fits 
    log10V = V0 + Va * (a - 1.) - Vz1 * np.log(a) + Vz2 * z
    log10eps = eps0_sfr + epsA_sfr * (a - 1) - epsLnA_sfr * np.log(a) + epsZ_sfr * z 
    alpha_sfr = alpha0_sfr + alphaA_sfr * (a - 1.) - alphaLnA_sfr * np.log(a) + alphaZ_sfr * z
    beta_sfr = beta0_sfr + betaA_sfr* (a - 1.) + betaZ_sfr * z
    delta_sfr = delta0_sfr 
    logGamma_sfr = gamma0_sfr + gammaA_sfr * (a - 1.) + gammaZ_sfr * z

    # Star Formation Rate 
    v = vMpeak / 10.0**(log10V)
    epsilon_sfr, gamma_sfr = 10.0**log10eps, 10.0**logGamma_sfr

    sfr = epsilon_sfr * (1./(v**alpha_sfr + v**beta_sfr)) 
    sfr += epsilon_sfr * (gamma_sfr * np.exp(- 0.5 * (np.log10(v)/delta_sfr)**2))

    return sfr 


@dataclass
class UniverseMachineRelations(StellarProperties): # See if you can use tuple instead
    # Instance Variables 
    # Mhalo : InitVar[float]
    # z : InitVar[float]   
    Mpeak : InitVar[float]

    # Fit Parameters (Look into Dictionary as Class Variables)
    stellarMassFitParameters : ClassVar[Dict[str, float]] = field(default=dict)
    sfrFitParameters : ClassVar[Dict[str, float]] = field(default={
        "V" : (2.151, 1.658, 1.68, -0.233),
        "eps" : (0.109, 3.441, 5.079, -0.781), 
        "alpha" : (-5.598, 20.731, 13.455, -1.321),
        "beta" : (-1.911, -0.395, 0.747), 
        "gamma" : (-1.699, -4.206, -0.809), 
        "delta" : 0.055
    })

    # Stellar Relations 
    vMpeak : float = field(default=0.0)

    # Load Values at Creation, to reduce computation during runtime 
    def __post_init__(self, Mhalo : float, z : float, Mpeak : float) -> None:

        # Load Stellar Mass Fit Parameters 
        M0, Ma, Mlna, Mz = self.stellarMassFitParameters["M"]
        eps0, epsA, epsLnA, epsZ = self.stellarMassFitParameters["eps"]
        alpha0, alphaA, alphaLnA, alphaZ = self.stellarMassFitParameters["alpha"]
        beta0, betaA, betaZ = self.stellarMassFitParameters["beta"]
        gamma0, gammaA, gammaZ = self.stellarMassFitParameters["gamma"]
        delta0 = self.stellarMassFitParameters["delta"]

        # Convert to Scale Factor
        a = 1.0/(1.0 + z)

        # Fit Relations 
        log10M1 = M0 + Ma * (a - 1.) - Mlna * np.log(a) + Mz * z 
        epsilon = eps0 + epsA * (a - 1) - epsLnA * np.log(a) + epsZ * z 
        alpha = alpha0 + alphaA * (a - 1.) - alphaLnA * np.log(a) + alphaZ * z
        beta = beta0 + betaA * (a - 1.) + betaZ * z
        delta = delta0 
        logGamma = gamma0 + gammaA * (a - 1.) + gammaZ * z

        x = np.log10(Mpeak) - log10M1
        gamma = 10.0**logGamma

        logM = log10M1 - np.log10(10.0**(-alpha * x) + 10.0**(-beta * x))
        logStellarMass = logM + epsilon + gamma * np.exp(-0.5 * (x/delta)**2)

        self.stellarMass = 10.0**(logStellarMass)

        # Bulge Mass 
        f2 = 0.5 * (z + 2.)/(z + 1.)
        self.bulgeMass = f2 * self.stellarMass/(1. + np.exp(-1.13 * (np.log10(self.stellarMass) - 10.2)))

        # Velocity of Peak Mass 
        M200 = 1.64E12/((a/0.378)**(-0.142) + (a/0.378)**(-1.79))
        self.vMpeak = 200.0 * (Mhalo/M200)**(1./3.)

        # Star Formation Rate Fit Parameters 
        V0, Va, Vz1, Vz2 = self.sfrFitParameters["V"]
        eps0_sfr, epsA_sfr, epsLnA_sfr, epsZ_sfr =self.sfrFitParameters["eps"]
        alpha0_sfr, alphaA_sfr, alphaLnA_sfr, alphaZ_sfr = self.sfrFitParameters["alpha"]
        beta0_sfr, betaA_sfr, betaZ_sfr = self.sfrFitParameters["beta"]
        gamma0_sfr, gammaA_sfr, gammaZ_sfr = self.sfrFitParameters["gamma"]
        delta0_sfr = self.sfrFitParameters["delta"]

        # Star Formation Rate Fits 
        log10V = V0 + Va * (a - 1.) - Vz1 * np.log(a) + Vz2 * z
        log10eps = eps0_sfr + epsA_sfr * (a - 1) - epsLnA_sfr * np.log(a) + epsZ_sfr * z 
        alpha_sfr = alpha0_sfr + alphaA_sfr * (a - 1.) - alphaLnA_sfr * np.log(a) + alphaZ_sfr * z
        beta_sfr = beta0_sfr + betaA_sfr* (a - 1.) + betaZ_sfr * z
        delta_sfr = delta0_sfr 
        logGamma_sfr = gamma0_sfr + gammaA_sfr * (a - 1.) + gammaZ_sfr * z

        # Star Formation Rate 
        v = self.vMpeak / 10.0**(log10V)
        epsilon_sfr, gamma_sfr = 10.0**log10eps, 10.0**logGamma_sfr

        self.starFormationRate = epsilon_sfr * (1./(v**alpha_sfr + v**beta_sfr) + gamma_sfr * np.exp(- 0.5 * (np.log10(v)/delta_sfr)**2))


@dataclass
class CentralsSatellites(UniverseMachineRelations):
    # Stellar Mass Fit Parameters (Class Variables)
    stellarMassFitParameters : ClassVar[Dict[str, float]] = field(default={
        "M" : (12.035, 4.556, 4.417, -0.731), 
        "eps" : (-1.435, 1.831, 1.368, -0.217),
        "alpha" : (1.963, -2.316, -1.732, 0.178),
        "beta" : (0.482, -0.841, -0.471), 
        "gamma" : (-1.034, -3.1, -1.055),
        "delta" : 0.411 
    })
    

@dataclass
class CentralsOnly(UniverseMachineRelations):
    # Fit Parameters (Class Variables)
    stellarMassFitParameters : ClassVar[Dict[str, float]] = field(default={
        "M" : (12.081, 4.696, 4.458, -0.740),
        "eps" : (-1.435, 1.813, 1.353, -0.214),
        "alpha" : (1.957, -2.65, -1.953, 0.204),
        "beta" : (0.474, -0.903, -0.492),
        "gamma" : (-1.065, -3.243, -1.107),
        "delta" : 0.386
    })
    

@dataclass
class SatellitesOnly(UniverseMachineRelations):
    # Fit Parameters (Class Variables)
    stellarMassFitParameters : ClassVar[Dict[str, float]] = field(default={
        "M" : (11.896, 3.284, 3.413, -0.58),
        "eps" : (-1.449, -1.256, -1.031, 0.108),
        "alpha" : (1.949, -4.096, -3.226, 0.401),
        "beta" : (0.477, 0.046, -0.214),
        "gamma" : (-0.755, 0.461, 0.025),
        "delta" : 0.357 
    })
    


if __name__ == "__main__":
    testAllClass = CentralsSatellites
    print(f"{testAllClass = }")
    print(f"{testAllClass.stellarMassFitParameters = }\n")
    # print(f"{testAllClass.sfrFitParameters = }\n")

    testCentralClass = CentralsOnly
    print(f"{testCentralClass = }")
    print(f"{testCentralClass.stellarMassFitParameters = }\n")
    # print(f"{testCentralClass.sfrFitParameters = }\n")

    testSatelliteClass = SatellitesOnly
    print(f"{testSatelliteClass = }")
    print(f"{testSatelliteClass.stellarMassFitParameters = }\n")
    # print(f"{testSatelliteClass.sfrFitParameters = }\n")

    assert (testAllClass.sfrFitParameters == testCentralClass.sfrFitParameters == testSatelliteClass.sfrFitParameters)

    # testAll = CentralsSatellites(1.0E10, 2.0, 1.0E12)
    # print(f"{testAll.stellarMassFitParameters = }")
    # print(f"{testAll.sfrFitParameters = }")
    # testCentrals = CentralsOnly()
    # print(f"{testCentrals = }")
    # testSatellites = SatellitessOnly() 
    # print(f"{testSatellites = }")