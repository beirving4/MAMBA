'''
       Log of Important Values, Constants, and Parameters Needed for Computation
    ------------------------------------------------------------------------------
'''

# Import Necessary Libraries
from math import pi, exp, sqrt

# Fundamental Physics Constats in SI
c_LIGHT = 2.99792458e8                  # Speed of Light [m/s]
h_PLANCK = 6.6260755e-34                # Planck Constant [J-s]
G_Newton = 6.67259e-11                  # Gravitation Constant [m^3/(kg*s^2)]
kB = 1.380658e-23                       # Boltzmann Constant [J/K]
R_gas = 8.31451                         # Molar Gas Constant [J/K]
N_Avg = 6.0221367e23                    # Avogadro's Number [1/mol]
e_Volt = 1.602176634e-19                # Electron Volt [J]/[C]
eps_Vacuum = 8.854187812813E-12         # Permittivity of Vacuum [F/m]
K_coulomb = 1/(4*pi*eps_Vacuum)         # Coulomb's Constant [(N*m^2)/C^2]
mu_Vacuum = (4*pi)*1.0E-7               # Permeability of Vacuum [H/m]

# Important Mass Scales
amu = 1.6605402E-27                     # Atomic Mass Unit [kg]
electron_mass = .0005486 * amu          # Electron Rest Mass [kg]
proton_mass = 1.0072765 * amu           # Proton Rest Mass [kg]
neutron_mass = 1.0086649 * amu          # Nuetron Mass [kg]
muon_mass = 0.1134288 * amu             # Muon Mass [kg]
hydrogen_mass = 1.00784 * amu           # Hydrogen Mass 
M_solar = 1.98892e30                    # Solar Mass [kg]
M_Earth = 5.972e24                      # Earth Mass [kg]
M_Lunar = 7.34767309E22                 # Lunar Mass [kg]
Chandrasekhar = 2.765E30                # Chandrasekhar Limit [kg]
Chandra = Chandrasekhar/M_solar         # Chandrasekhar Limit [Msol]

# Earth Surface/Atmospheric Values
P_atm = 101325                          # Atmospheric Pressure [Pa/m^2]
g_Earth = 9.80665                       # Earth's Surface Gravity [m/s^2]

# Important Length Scales
ast_u = 1.495978707e11                  # Astronomical Unit [m]
light_yr = 9.460532621E15               # Light Year [m]
parsec = 3.08567808E16                  # Parsec [m]
kiloparsec = 1.0E3*parsec               # Kiloparsec [m]
megaparsec = 1.0E6*parsec               # Megaparsec [m]
R_solar = 6.957E8                       # Solar Radius [m]

# Cosmological Parameters (Reported by Planck 2018)
h_Hubble = 0.674                        # Little Hubble
Hubble = h_Hubble*(1.0E5/megaparsec)    # Hubble Constant [1/s]
Omega_bary = 0.0493                     # Cosmological Baryon Density
Omega_DM = 0.264                        # Cosmological Dark Matter Density
Omega_DE = 0.69                         # Cosmological Dark Energy Density
Omega_k =  0.001                        # Cosmological Curvature Density
Omega_m = 0.3158                        # Cosmological Matter Density
sigma_8 = 0.811                         # Cosmological Matter Fluctuation Amplitude
w_EOS = -1.03                           # Cosmological Equation of State Parameter
n_spectral = 0.965                      # Scalar Spectral Index
f_bary = Omega_bary/Omega_m             # Cosmological Baryonic Fraction 

# Important Time Scales
julian_day = 86400                      # Julian Day [s]
julian_year = julian_day*365.25         # Julian Year [s]
megaYear = 1.0E6 * julian_year          # Megayear [Myr] Timescale 
gigaYear = 1.0E9 * julian_year          # Gigayear [Myr] Timescale

# Important Quantum Values
b_wein = 2.897756E-3                                    # Wein Displacement Constant [m-K]
Bohr_magneton = 9.2740154e-28                           # Bohr Magneton [(J-G)/s]
thomp_cross = 6.6524587158E-29                          # Thomson Cross-Section [m^2]
hbar = h_PLANCK/(2*pi)                                  # Reduced Planck Constant [J-s]
stfB = (pi**2*kB**4)/(60*hbar**3*c_LIGHT*2)             # Stefan-Boltzmann Constant [erg/(s*cm^2*K^4)]
Phi0 = h_PLANCK/(2*exp(1))                              # Magnetic Flux Quantum [(m^(3/2)kg^(1/2))/s]
a0 = (hbar**2)/(electron_mass*e_Volt**2)                # Bohr Radius [m]
Rinf = e_Volt**2/(2*a0)                                 # Rydberg Constant [1/cm]
F_const = e_Volt*N_Avg                                  # Faraday's Constant [(m^(3/2)kg^(1/2))/(s-mol)]
fine_struct = (K_coulomb*e_Volt**2)/(hbar*c_LIGHT)      # Fine Structure Constant

# Planck Scale Values
l_planck = 1.616E-39                        # Planck Length [m]
m_planck = sqrt(hbar*c_LIGHT/G_Newton)      # Planck Mass [kg]
t_planck = 5.391E-44                        # Planck Time [s]
q_planck = 5.623E-9                         # Planck Charge [Q]
T_planck = 1.417E32                         # Planck Temperature [K]

# Important Radiative Transfer Values
thomp_absorp = thomp_cross/proton_mass                  # Thomson Absorption Coefficient for Pure Hydrogen
t_eddington = (thomp_absorp*c_LIGHT)/(4*pi*G_Newton)    # Eddington Accretion Time Scale [s]
t_eddington_Myr = t_eddington/megaYear                  # Eddington Accretion Time Scale [Myr]
L_solar = 3.828E26                                      # Solar Luminosity [W]

# SI -> Geometrized Units (G = c = 1) - All in meters!!
time_gu = c_LIGHT                                       # Time -> Length Conversion (~4.925E-6 s for 1 M_Solar) -> Test
velocity_gu = 1/c_LIGHT                                 # Velocity/Angular Velocity -> Length Conversion
acceleration_gu = 1/c_LIGHT**2                          # Acceleration -> Length Conversion
mass_gu = G_Newton/c_LIGHT**2                           # Mass/Mass Density -> Length Conversion (50% of Sun Swarzchild Radius ~1.48km)
energy_gu = G_Newton/c_LIGHT**4                         # Energy/Energy Density -> Length Conversion
momentum_gu = G_Newton/c_LIGHT**3                       # Momentum/Angular Momentum -> Length Conversion
pressure_gu = G_Newton/c_LIGHT**4                       # Pressure -> Length Conversion
force_gu = G_Newton/c_LIGHT**4                          # Force -> Length Conversion

# Power -> Length Conversion
power_gu = G_Newton/c_LIGHT**5                          # Electric Charge -> Length Conversion
charge_gu = sqrt(G_Newton*K_coulomb)/c_LIGHT**2         # Electric Field -> Length Conversion
eField_gu = sqrt(G_Newton/K_coulomb)/c_LIGHT**2         # Electric Potential -> Length Conversion
Vele_gu = sqrt(G_Newton/K_coulomb)/c_LIGHT**2           # Magnetic Field & Potential -> Length Conversion
bField_gu = sqrt(G_Newton/K_coulomb)/c_LIGHT            # Magnetic Potential -> Length Conversion
Vele_gu = sqrt(G_Newton/K_coulomb)/c_LIGHT              # Energy Potential -> Length Conversion
potential_gu = sqrt(G_Newton/K_coulomb)/c_LIGHT         # Sun Swarzchild Radius (~2.95E3m for 1 M_Solar) -> Test
Rs_solar = 2*mass_gu                                    # Sun Swarzchild Radius (~2.95km for 1 M_Solar) -> Test
Rs_solar_km = Rs_solar/1E3

# Gravitational Units (G = c = M = 1) - Input Mass := M_solar
length_gru = mass_gu*M_solar           # Length Conversion [m]
area_gru = length_gru**2               # Area Conversion [m^2]
volume_gru = length_gru**3             # Volume Conversion [m^3]
spin_gru = 1/(M_solar*momentum_gu)     # Angular Momentum Conversion [rad/s]
# acceleration_gru = M_solar/mass_gu   # Acceleration Conversion

# In Case You Need to Check a Value
if __name__ == '__main__':
    # print(t_eddington_Myr)
    print(Hubble)
    # print('time_gu = {:.3e} [s]\n'.format(time_gu*length_gu*M_solar))
    # print('f = 1/2dt = {:.3e} \n'.format(1.0/(4*julian_day)))
    # print('H0 = {:.4e}\n'.format(G_Newton * (M_solar / 1.0E6) / kiloparsec))
    # print((1.0E5/megaparsec))
    # print('G^(5/3)/c^(4) = {:.3e} \n'.format(G_Newton**(5./3)/c_LIGHT**(4)))
    # print('t_eddington = {:.3e} [Myrs]\n'.format(1.0E-6*t_eddington/julian_year))
    
# Unit Conversion Factors
# MeV_2_ergs = 1.60217733E-6                               # MeV to ergs Unit Conversion
# MeVfm3_2_gcm3 = 1.7827E12                                # MeV/fm3 to g/cm3 Conversion
# MeVfm3_2_dynecm2 = 1.60217733E33      		             # MeV/fm3 to dyne/cm2 Conversion 
# fm_2_cm = 1E-13                                          # femtometers (1E-15 m) to Centimeter (1E-2 m) Conversion
# cm_2_km = 1E-5                                           # cm to km 
# cm2_2_km2 = 1E-10                                        # cm^2 to km^2 
# cm3_2_km3 = 1E-15                                        # cm^3 to km^3
# # GigaeV-femtometer/(dyne-cm)
# GeVfm_per_dynecm = (1E19*e_Volt)/(fm_2_cm)**3
# MeVfm3_2_km2 = (G/c**2)*(MeVfm3_2_gcm3/cm2_2_km2)        # MeV/fm3 to km^-2 for Geometrized Units [g/cm^3]
# # MeVfm3_2_km2_2 = (G/c**4)*(MeVfm3_2_dynecm2/cm2_2_km2) # MeV/fm3 to km^-2 for Geometrized Units [dyne/cm^2]
# MeVfm3_2_gkm3 = MeVfm3_2_gcm3/cm3_2_km3                  # MeV/fm3 to g/km3 Conversion  
# km_2_Msol = (1/Msolar)*(c**2/G)*1E5                      # km to m/Msolar Conversion 

# Nuclear Units (c = 1, hbar = 197 MeV*fm, Energy := MeV, Length := femtometers)
# hbar_nu = (hbar*c)/(MeV_2_ergs*fm_2_cm)   # Reduced Planck Constant Conversion
# mass_nu = MeV_2_ergs/c**2                 # Mass Conversion 
# length_nu = fm_2_cm                       # Length Conversion 
# density_nu = fm_2_cm**(-3)                # Density Conversion 
# temperature_nu = MeV_2_ergs/kB            # Temperature Conversion
# mn_nu = (mn*c**2)/MeV_2_ergs              # Neutron Mass Conversion 
# me_nu = (mn*c**2)/MeV_2_ergs              # Electron Mass Conversion 
# mp_nu = (mn*c**2)/MeV_2_ergs              # Proton Mass Conversion 
# mm_nu = (mn*c**2)/MeV_2_ergs              # Muon Mass Conversion 
# amu_nu = (mn*c**2)/MeV_2_ergs             # Atomic Mass Unit Conversion 

# Nuclear Units to Gravitational Units Conversion 
# mn_nu2gu = mn*(density_nu/density_gu)     # Neutron Mass Conversion
# amu_nu2gu = amu*(density_nu/density_gu)   # Atomic Mass Unit Conversion

# In Case You Need to Check a Value
# print('MeVfm3_2_km2 [g/cc] = ',MeVfm3_2_km2)#_1,' MeVfm3_2_km2 [dyne/cm^2] = ',MeVfm3_2_km2_2)
# print(' G/c^2 [g.units] ', (G/c**2)/cm2_2_km2, '\n') 
# print('MeVfm3_2_km2/G/c^2 [g.units] [g/cc] = ',MeVfm3_2_km2)
# print("{:.3e}".format(10E-12/Gcc_2_km2)) 
# print(Chandra)
# print(Gcc_2_km2*1E9)
