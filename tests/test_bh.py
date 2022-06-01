# Unit Tests for Class Functionality 
# Import Libraries 
import pdb
import numpy as np
# import numpy.typing as npt 
# from copy import deepcopy 
# from constants import *
# from population import *
# from halo import *
from black_hole import BlackHole
from constants import * 

import astropy.units as u
import astropy.cosmology.units as cu
from astropy.cosmology import Planck18 as cosmo

import sys 
from pympler import asizeof 

# Red -> Green -> Refactor !
# UnitTest Class from unittest library 


def blackHoleTests():

    print("\n-------------------------")
    print("Tests for BlackHole Class")
    print("-------------------------\n")

    n_test = 7

    # Generate a Black Hole 
    M87 = BlackHole(6.5E9, 0.00428, spin=0.9)
    SgA = BlackHole(4.0E6, 0.00, spin=0.1)
    failedBH1 = BlackHole(0.0)
    failedBH2 = BlackHole(-1.0E6)
    q = 4.0E6/6.5E9 
    eta = q/(1 + q)**2

    # Print 
    print(f"\n\tM87* -> {M87}")
    print(f"\t{sys.getsizeof(M87)=}")
    print(f"\t{asizeof.asizeof(M87)=}")

    print(f"\tSgA* -> {SgA}") 

    assert failedBH1 is None, "Check the Mass Condition for AstroObject for Zero Mass."
    assert failedBH2 is None, "Check the Mass Condition for AstroObject for Negative Mass."

    print(f"\n\tTest (1/{n_test}): Class Initialization & Presentation --- PASS")

    # Check Operator Overloads 
    assert M87 > SgA, "Check (>) overload!"
    assert M87 >= SgA, "Check (>=) overload!"
    assert SgA < M87, "Check (<) overload!"
    assert SgA <= M87, "Check (<=) overload!"
    assert M87 != SgA, "Check (!=) overload!"
    assert M87 == M87, "Check (==) overload!"
    assert SgA == SgA, "Check (==) overload!"
    assert (SgA + SgA) == (SgA.mass + SgA.mass), "Check (+) operator overload!" 
    assert (M87 + M87) == (M87.mass + M87.mass), "Check (+) operator overload!" 
    assert M87.isFaster(SgA), "Check isFaster() function!"
    assert SgA.isCloser(M87), "Check isCloser() function!"
    assert SgA.sameSpeed(SgA), "Check sameSpeed() function!"
    assert M87.sameDistance(M87), "Check sameDistance() function!"
    assert M87.massRatio(SgA) == q, "Check Mass Ratio Function"
    assert SgA.massRatio(M87) == q, "Check Mass Ratio Function"
    assert M87.symmetricMassRatio(SgA) == eta, "Check Mass Ratio Function"
    assert SgA.symmetricMassRatio(M87) == eta, "Check Mass Ratio Function"

    BH1 = BlackHole(1.0E6, 1.0, 0.4)
    BH2 = BlackHole(1.5E6, 1.0, 0.4)
    BH3 = BlackHole(1.0E6, 2.0, 0.4)

    assert BH2 > BH1, "Check order overloads for mass comparison"
    assert BH3 > BH1, "Check order overloads for age comparison"

    print(f"\n\tTest (2/{n_test}): Operator Overloads --- PASS")

    t_grow = lambda t_sp, n: t_sp * n * np.log(10.0)

    # Check Accretion Growth (at Eddington Limit)
    z_test = 6.0

    for n in range(1, 8): 
        light_seed = BlackHole(1.0E3, z_test)
        light_seed.accretion_budget = 1.0E10 
        dt = t_grow(light_seed.salpeter_time, n)
        light_seed.accretionGrowth(dt, ratio=1.0)
        assert abs(np.log10(light_seed.mass/light_seed.seedMass) - n) <= 1.0E-15, "Check Accretion Growth" 


    print(f"\n\tTest (3/{n_test}): accretionGrowth() --- PASS")

    # Check Merger 
    seedBH = BlackHole(1.0E5, z_test) 
    small_partner = BlackHole(1.0E3, z_test)
    same_partner = BlackHole(1.0E5, z_test)
    large_partner = BlackHole(1.0E7, z_test)

    small_merge = seedBH.mergerGrowth(small_partner)
    same_merge = seedBH.mergerGrowth(same_partner)
    large_merge = seedBH.mergerGrowth(large_partner)

    inBetween = lambda x, a, b: min(a, b) <= x <= max(a, b)

    # Check to See if Mass and Spins Handled Correctly 
    assert (small_merge.mass == (seedBH + small_partner)), "Check Merger for smaller mass addition"
    assert (small_merge.seedMass == seedBH.seedMass), "Check Seed Mass Transfer for smaller mass addition"
    assert (inBetween(small_merge.spin, seedBH.spin, small_partner.spin)), "Check Spin Addition for smaller mass additon"
    assert (small_merge.fromBHMerger == small_partner.mass), "Check BH Merger Contribution Transfer for smaller mass addition"
    assert np.all(small_merge.cumulatives[z_test] == np.array([1., 0., 0., 0., 1.])), "Check Cumulative Totals Transfer"

    assert (same_merge.mass == (same_partner + seedBH) == (2.0 * seedBH.mass)), "Check Merger for same mass addition"
    assert (same_merge.seedMass == seedBH.seedMass), "Check Seed Mass Transfer for same mass addition"
    assert (inBetween(same_merge.spin, seedBH.spin, same_partner.spin)), "Check Spin Addition for same mass additon"
    assert (same_merge.fromBHMerger == same_partner.mass == seedBH.mass), "Check BH Merger Contribution Transfer for same mass addition"
    assert np.all(same_merge.cumulatives[z_test] == np.array([1., 0., 0., 0., 1.])), "Check Cumulative Totals Transfer"

    assert (large_merge.mass == (seedBH + large_partner)), "Check Merger for larger mass addition"
    assert (large_merge.seedMass == large_partner.seedMass), "Check Seed Mass Transfer for larger mass addition"
    assert (inBetween(large_merge.spin, seedBH.spin, large_partner.spin)), "Check Spin Addition for larger mass additon"
    assert (large_merge.fromBHMerger == seedBH.mass), "Check BH Merger Contribution Transfer for larger mass addition"
    assert np.all(large_merge.cumulatives[z_test] == np.array([1., 0., 0., 0., 1.])), "Check Cumulative Totals Transfer"


    print(f"\n\tTest (4/{n_test}): mergerGrowth() --- PASS")

    # Check Accretion then Merger  (This is to see if the cumulative values were handled correctly...)
    seedBH2 = BlackHole(1.0E4, z_test) 

    n_total = 2.0 
    n_trigger = n_total/2.0
    n_steady = n_total/2.0 

    seedBH2.accretion_budget = 1.0E10 

    trigMedd = seedBH2.eddingtonRate 
    dt_trigger = t_grow(seedBH2.salpeter_time, n_trigger)
    trigMass = dt_trigger * trigMedd
    seedBH2.accretionGrowth(dt_trigger, ratio=1.0, accretion_type=1)

    assert np.all(seedBH2.cumulatives[z_test] == np.array([0., trigMass, 0., trigMass, trigMass])/trigMass), "Check Cumulative Setup"
    assert (seedBH2.fromMergerTriggeredAccretion == trigMass), "Check Merger Triggered Mass Accretion Update"

    steadyMedd = seedBH2.eddingtonRate 
    dt_steady = t_grow(seedBH2.salpeter_time, n_steady)
    steadyMass = dt_steady * steadyMedd
    seedBH2.accretionGrowth(dt_steady, ratio=1.0, accretion_type=2)

    total_accretion = trigMass + steadyMass
    accretionCumulatives = np.array([0.0, trigMass, steadyMass, total_accretion, trigMass])/total_accretion

    assert (seedBH2.fromSteadyModeAccretion == steadyMass), "Check Steady Mass Accretion Update" 
    assert (seedBH2.fromAccretion == total_accretion), "Check fromAccretion calculation"
    assert abs(np.log10(seedBH2.mass/seedBH2.seedMass) - n_total) <= 1.0E-15, "Check Accretion Growth" 
    assert np.allclose(seedBH2.cumulatives[z_test], accretionCumulatives), "Check Cumulative Calculations"

    same_partner_2 = BlackHole(1.0E6, 6.0)

    assert np.allclose(same_partner_2.mass, seedBH2.mass)

    small_merge_2 = seedBH2.mergerGrowth(small_partner)
    same_merge_2 = seedBH2.mergerGrowth(same_partner_2)
    large_merge_2 = seedBH2.mergerGrowth(large_partner)  

    small_totals = total_accretion + small_partner.mass 
    same_totals = total_accretion + seedBH2.mass
    assert np.allclose(same_totals, (total_accretion + same_partner_2.mass))
    large_totals = seedBH2.mass 

    small_cumulatives = np.array([small_partner.mass , trigMass, steadyMass, total_accretion, trigMass + small_partner.mass])/small_totals
    same_cumulatives = np.array([same_partner_2.mass , trigMass, steadyMass, total_accretion, trigMass + same_partner_2.mass])/same_totals 
    large_cumulatives = np.array([seedBH2.mass , 0., 0., 0., seedBH2.mass])/large_totals

    assert (small_merge_2.mass == (seedBH2 + small_partner)), "Check Merger for smaller mass addition"
    assert (small_merge_2.seedMass == seedBH2.seedMass), "Check Seed Mass Transfer for smaller mass addition"
    assert (inBetween(small_merge_2.spin, seedBH2.spin, small_partner.spin)), "Check Spin Addition for smaller mass additon"
    assert (small_merge_2.fromBHMerger == small_partner.mass), "Check BH Merger Contribution Transfer for smaller mass addition"
    assert (small_merge_2.fromMergerTriggeredAccretion == seedBH2.fromMergerTriggeredAccretion), "Check Mass Transfer for smaller mass addition"
    assert (small_merge_2.fromSteadyModeAccretion == seedBH2.fromSteadyModeAccretion), "Check Mass Transfer for smaller mass addition"
    assert (small_merge_2.fromAccretion == seedBH2.fromAccretion), "Check Mass Transfer for smaller mass addition"
    assert (small_merge_2.fromMergerEffects == (seedBH2.fromMergerTriggeredAccretion + small_partner.mass)), "Check fromMergerEffects calculation for smaller mass addition"
    assert np.allclose(small_merge_2.cumulatives[z_test], small_cumulatives), "Check Cumulative Generation and Transfer for smaller mass addition"

    assert np.isclose(same_merge_2.mass, (seedBH2 + same_partner_2)), "Check Merger for same mass addition"
    assert np.isclose(same_merge_2.mass, 2.0 * seedBH2.mass), "Check Merger for same mass addition"
    assert np.isclose(same_merge_2.seedMass, seedBH2.seedMass), "Check Seed Mass Transfer for same mass addition"
    assert np.isclose(same_merge_2.seedRedshift, seedBH2.seedRedshift), "Check Seed Redshift Transfer for same mass addition"
    assert (inBetween(same_merge_2.spin, seedBH2.spin, same_partner_2.spin)), "Check Spin Addition for same mass additon"
    assert np.isclose(same_merge_2.fromBHMerger, seedBH2.mass), "Check BH Merger Contribution Transfer for same mass addition"
    assert (same_merge_2.fromMergerTriggeredAccretion == seedBH2.fromMergerTriggeredAccretion), "Check Mass Transfer for same mass addition"
    assert (same_merge_2.fromSteadyModeAccretion == seedBH2.fromSteadyModeAccretion), "Check Mass Transfer for same mass addition"
    assert (same_merge_2.fromAccretion == seedBH2.fromAccretion), "Check Mass Transfer for same mass addition"
    assert np.isclose(same_merge_2.fromMergerEffects, (seedBH2.fromMergerTriggeredAccretion + seedBH2.mass)), "Check fromMergerEffects calculation for same mass addition"
    assert np.allclose(same_merge_2.cumulatives[z_test], same_cumulatives), "Check Cumulative Generation and Transfer for same mass addition"

    assert (large_merge_2.mass == (seedBH2 + large_partner)), "Check Merger for larger mass addition"
    assert (large_merge_2.seedMass == large_partner.seedMass), "Check Seed Mass Transfer for larger mass addition"
    assert (large_merge_2.seedRedshift == large_partner.seedRedshift), "Check Seed Redshift Transfer for larger mass addition"
    assert (inBetween(large_merge_2.spin, seedBH2.spin, large_partner.spin)), "Check Spin Addition for larger mass additon"
    assert (large_merge_2.fromBHMerger == seedBH2.mass), "Check BH Merger Contribution Transfer for larger mass addition"
    assert (large_merge_2.fromMergerTriggeredAccretion == large_partner.fromMergerTriggeredAccretion), "Check Mass Transfer for larger mass addition"
    assert (large_merge_2.fromSteadyModeAccretion == large_partner.fromSteadyModeAccretion), "Check Mass Transfer for large mass addition"
    assert (large_merge_2.fromAccretion == large_partner.fromAccretion), "Check Mass Transfer for large mass addition"
    assert (large_merge_2.fromMergerEffects == (large_partner.fromMergerTriggeredAccretion + seedBH2.mass)), "Check fromMergerEffects calculation for smaller mass addition"
    assert np.allclose(large_merge_2.cumulatives[z_test], large_cumulatives), "Check Cumulative Generation and Transfer for larger mass addition"


    print(f"\n\tTest (5/{n_test}): accretionGrowth -> mergerGrowth() -> Information Transfers --- PASS")

    # Check For Time Evolution 
    redshifts = np.arange(10.0, 5.0, -1.)
    seedBH3 = BlackHole(1.0E4, redshifts[0])
    seedBH3.accretion_budget = 1.0E10 

    # Run the Case for Same BH Time Evolution 
    for i, z in enumerate(redshifts): 
        if (i > 0):
            dt = 1000.0 * abs(float((cosmo.age(z) - cosmo.age(redshifts[i-1]))/u.Gyr)) 
            seedBH3.redshift = z
            seedBH3.accretionGrowth(dt, ratio=1.0, accretion_type=(1 + i % 2)) 
            if (i == 4): 
                merge_result = seedBH3.mergerGrowth(small_partner) 
                assert (merge_result.redshift == z), "Check Transfer of Redshift Info"
                assert (merge_result.seedMass == seedBH3.seedMass), "Check Seed Info Transfer"
                assert (merge_result.seedRedshift == seedBH3.seedRedshift), "Check Seed Info Transfer"
                assert (np.count_nonzero(merge_result.cumulatives[z]) == 5.0)
            
            if (i == 1): assert (np.count_nonzero(seedBH3.cumulatives[z]) == 2.0)
            if (i in [2,3]): assert (np.count_nonzero(seedBH3.cumulatives[z]) == 4.0)
            if (i < 4): assert (np.array(list(seedBH3.cumulatives.keys())) == redshifts[1:(i+1)]).all()


    assert (np.array(list(merge_result.cumulatives.keys())) == redshifts[1:]).all()

    print(f"\n\tTest (6/{n_test}): Time-Evolution Test --- PASS")


    # Verify Bolometric & X-ray Luminosity Conversions Are Correct 
    kx = lambda x: 10.83 * x**0.28 + 6.08 * x**(-0.02) 
    mass_val = lambda frac: frac * ((1.0E10 * L_solar * t_eddington)/(M_solar * c_LIGHT**2)) 
    bolLumin = lambda mass: (mass * M_solar * c_LIGHT**2/t_eddington)/L_solar  
    uvLumin = lambda fuv, mass: 0.5 * fuv * bolLumin(mass)

    # Check Parameters Based on Model Relations for Varying Values 
    for frac_val in np.logspace(-5, 5, 100):
        # Initialize BH 
        MBH = BlackHole(mass_val(frac_val), massAccretionRatio=1.0) 

        # Check Bolometric, UV & X-Ray Luminosity Calculations
        assert MBH.massAccretionRatio == 1.0, "Check f_Edd update rule"
        assert np.isclose(MBH.bolometricLuminosity, bolLumin(MBH.mass)), "Check Bolometric Lumonisty Calculation"
        assert np.isclose(MBH.uvLuminosity, uvLumin(0.1, MBH.mass)), "Check UV Lumonisty Calculation"
        assert np.isclose(MBH.bolometricLuminosity/MBH.xrayLuminosity, kx(frac_val)), "Check X-Ray Luminosity Calculation"
    
        # 
    
    # Gravitational Wave Recoil Velocity Tests

    # Check Drawing Parameters for Different Accretion Distributions 

    print(f"\n\tTest (7/{n_test}): Model Relations Test --- PASS")

    print("\n")


    if __name__ == "__main__":
        blackHoleTests() 