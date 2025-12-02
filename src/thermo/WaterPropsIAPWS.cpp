/**
 * @file WaterPropsIAPWS.cpp
 * Definitions for a class for calculating the equation of state of water
 * from the IAPWS 1995 Formulation based on the steam tables thermodynamic
 * basis (See class @link Cantera::WaterPropsIAPWS WaterPropsIAPWS@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{
// Critical Point values of water in mks units

//! Critical Temperature value (kelvin)
const CanteraDouble T_c = 647.096;
//! Critical Pressure (Pascals)
static const CanteraDouble P_c = 22.064E6;
//! Value of the Density at the critical point (kg m-3)
const CanteraDouble Rho_c = 322.;

static const CanteraDouble R_water = 461.51805; // J/kg/K (Eq. 6.3)

void WaterPropsIAPWS::calcDim(CanteraDouble temperature, CanteraDouble rho)
{
    tau = T_c / temperature;
    delta = rho / Rho_c;

    // Determine the internal state
    if (temperature > T_c) {
        iState = WATER_SUPERCRIT;
    } else {
        if (delta < 1.0) {
            iState = WATER_GAS;
        } else {
            iState = WATER_LIQUID;
        }
    }
}

CanteraDouble WaterPropsIAPWS::pressure() const
{
    CanteraDouble retn = m_phi.pressureM_rhoRT(tau, delta);
    CanteraDouble rho = delta * Rho_c;
    CanteraDouble temperature = T_c / tau;
    return retn * rho * R_water * temperature;
}

CanteraDouble WaterPropsIAPWS::density(CanteraDouble temperature, CanteraDouble pressure,
                                int phase, CanteraDouble rhoguess)
{
    if (fabs(pressure - P_c) / P_c < 1.e-8 &&
        fabs(temperature - T_c) / T_c < 1.e-8) {
        // Catch critical point, as no solution is found otherwise
        setState_TD(temperature, Rho_c);
        return Rho_c;
    }
    CanteraDouble deltaGuess = 0.0;
    if (rhoguess == -1.0) {
        if (phase != -1) {
            if (temperature > T_c) {
                rhoguess = pressure / (R_water * temperature);
            } else {
                if (phase == WATER_GAS || phase == WATER_SUPERCRIT) {
                    rhoguess = pressure / (R_water * temperature);
                } else if (phase == WATER_LIQUID) {
                    // Provide a guess about the liquid density that is
                    // relatively high -> convergence from above seems robust.
                    rhoguess = 1000.;
                } else if (phase == WATER_UNSTABLELIQUID || phase == WATER_UNSTABLEGAS) {
                    throw CanteraError("WaterPropsIAPWS::density",
                                       "Unstable Branch finder is untested");
                } else {
                    throw CanteraError("WaterPropsIAPWS::density",
                                       "unknown state: {}", phase);
                }
            }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to
            // the routine
            rhoguess = pressure / (R_water * temperature);
        }
    }
    CanteraDouble p_red = pressure / (R_water * temperature * Rho_c);
    deltaGuess = rhoguess / Rho_c;
    setState_TD(temperature, rhoguess);
    CanteraDouble delta_retn = m_phi.dfind(p_red, tau, deltaGuess);
    if (delta_retn <= 0) {
        // No solution found for first initial guess; perturb initial guess once
        // to avoid spurious failures (band-aid fix)
        delta_retn = m_phi.dfind(p_red, tau, 0.9 * deltaGuess);
    }
    CanteraDouble density_retn;
    if (delta_retn > 0.0) {
        delta = delta_retn;

        // Dimensionalize the density before returning
        density_retn = delta_retn * Rho_c;

        // Set the internal state -> this may be a duplication. However, let's
        // just be sure.
        setState_TD(temperature, density_retn);
    } else {
        density_retn = -1.0;
    }
    return density_retn;
}

CanteraDouble WaterPropsIAPWS::density() const
{
    return delta * Rho_c;
}

CanteraDouble WaterPropsIAPWS::temperature() const
{
    return T_c / tau;
}

CanteraDouble WaterPropsIAPWS::psat_est(CanteraDouble temperature) const
{
    // Formula and constants from: "NBS/NRC Steam Tables: Thermodynamic and
    // Transport Properties and Computer Programs for Vapor and Liquid States of
    // Water in SI Units". L. Haar, J. S. Gallagher, G. S. Kell. Hemisphere
    // Publishing. 1984.
    static const CanteraDouble A[8] = {
        -7.8889166E0,
        2.5514255E0,
        -6.716169E0,
        33.2239495E0,
        -105.38479E0,
        174.35319E0,
        -148.39348E0,
        48.631602E0
    };
    CanteraDouble ps;
    if (temperature < 314.) {
        CanteraDouble pl = 6.3573118E0 - 8858.843E0 / temperature
                        + 607.56335E0 * pow(temperature, -0.6);
        ps = 0.1 * exp(pl);
    } else {
        CanteraDouble v = temperature / 647.25;
        CanteraDouble w = fabs(1.0-v);
        CanteraDouble b = 0.0;
        for (int i = 0; i < 8; i++) {
            CanteraDouble z = i + 1;
            b += A[i] * pow(w, ((z+1.0)/2.0));
        }
        CanteraDouble q = b / v;
        ps = 22.093*exp(q);
    }

    // Original correlation was in cgs. Convert to mks
    ps *= 1.0E6;
    return ps;
}

CanteraDouble WaterPropsIAPWS::isothermalCompressibility() const
{
    CanteraDouble dpdrho_val = dpdrho();
    CanteraDouble dens = delta * Rho_c;
    return 1.0 / (dens * dpdrho_val);
}

CanteraDouble WaterPropsIAPWS::dpdrho() const
{
    CanteraDouble retn = m_phi.dimdpdrho(tau, delta);
    CanteraDouble temperature = T_c/tau;
    return retn * R_water * temperature;
}

CanteraDouble WaterPropsIAPWS::coeffPresExp() const
{
    return m_phi.dimdpdT(tau, delta);
}

CanteraDouble WaterPropsIAPWS::coeffThermExp() const
{
    CanteraDouble kappa = isothermalCompressibility();
    CanteraDouble beta = coeffPresExp();
    CanteraDouble dens = delta * Rho_c;
    return kappa * dens * R_water * beta;
}

void WaterPropsIAPWS::corr(CanteraDouble temperature, CanteraDouble pressure,
    CanteraDouble& densLiq, CanteraDouble& densGas, CanteraDouble& delGRT)
{
    densLiq = density(temperature, pressure, WATER_LIQUID, densLiq);
    if (densLiq <= 0.0) {
        throw CanteraError("WaterPropsIAPWS::corr",
            "Error occurred trying to find liquid density at (T,P) = {}  {}",
            temperature, pressure);
    }
    setState_TD(temperature, densLiq);
    CanteraDouble gibbsLiqRT = m_phi.gibbs_RT();

    densGas = density(temperature, pressure, WATER_GAS, densGas);
    if (densGas <= 0.0) {
        throw CanteraError("WaterPropsIAPWS::corr",
            "Error occurred trying to find gas density at (T,P) = {}  {}",
            temperature, pressure);
    }
    setState_TD(temperature, densGas);
    CanteraDouble gibbsGasRT = m_phi.gibbs_RT();

    delGRT = gibbsLiqRT - gibbsGasRT;
}

CanteraDouble WaterPropsIAPWS::psat(CanteraDouble temperature, int waterState)
{
    CanteraDouble densLiq = -1.0, densGas = -1.0, delGRT = 0.0;
    if (temperature >= T_c) {
        densGas = density(temperature, P_c, WATER_SUPERCRIT);
        setState_TD(temperature, densGas);
        return P_c;
    }
    CanteraDouble p = psat_est(temperature);
    for (int i = 0; i < 30; i++) {
        corr(temperature, p, densLiq, densGas, delGRT);
        CanteraDouble delV = 1.0/densLiq - 1.0/densGas;
        p -= delGRT * R_water * temperature / delV;
        if (delGRT < 1.0E-8) {
            break;
        }
    }
    // Put the fluid in the desired end condition
    if (waterState == WATER_LIQUID) {
        setState_TD(temperature, densLiq);
    } else if (waterState == WATER_GAS) {
        setState_TD(temperature, densGas);
    } else {
        throw CanteraError("WaterPropsIAPWS::psat",
                           "unknown water state input: {}", waterState);
    }
    return p;
}

int WaterPropsIAPWS::phaseState(bool checkState) const
{
    if (checkState) {
        if (tau <= 1.0) {
            iState = WATER_SUPERCRIT;
        } else {
            CanteraDouble T = T_c / tau;
            CanteraDouble rho = delta * Rho_c;
            CanteraDouble rhoMidAtm = 0.5 * (OneAtm / (R_water * 373.15) + 1.0E3);
            CanteraDouble rhoMid = Rho_c + (T - T_c) * (Rho_c - rhoMidAtm) / (T_c - 373.15);
            int iStateGuess = WATER_LIQUID;
            if (rho < rhoMid) {
                iStateGuess = WATER_GAS;
            }
            CanteraDouble kappa = isothermalCompressibility();
            if (kappa >= 0.0) {
                iState = iStateGuess;
            } else {
                // When we are here we are between the spinodal curves
                CanteraDouble rhoDel = rho * 1.000001;
                CanteraDouble deltaSave = delta;
                CanteraDouble deltaDel = rhoDel / Rho_c;
                delta = deltaDel;
                m_phi.tdpolycalc(tau, deltaDel);

                CanteraDouble kappaDel = isothermalCompressibility();
                CanteraDouble d2rhodp2 = (rhoDel * kappaDel - rho * kappa) / (rhoDel - rho);
                if (d2rhodp2 > 0.0) {
                    iState = WATER_UNSTABLELIQUID;
                } else {
                    iState = WATER_UNSTABLEGAS;
                }
                delta = deltaSave;
                m_phi.tdpolycalc(tau, delta);
            }
        }
    }
    return iState;
}

void WaterPropsIAPWS::setState_TD(CanteraDouble temperature, CanteraDouble rho)
{
    calcDim(temperature, rho);
    m_phi.tdpolycalc(tau, delta);
}

CanteraDouble WaterPropsIAPWS::gibbs_mass() const
{
    return m_phi.gibbs_RT() * R_water * T_c / tau;
}

CanteraDouble WaterPropsIAPWS::enthalpy_mass() const
{
    return m_phi.enthalpy_RT() * R_water * T_c / tau;
}

CanteraDouble WaterPropsIAPWS::intEnergy_mass() const
{
    return m_phi.intEnergy_RT() * R_water * T_c / tau;
}

CanteraDouble WaterPropsIAPWS::entropy_mass() const
{
    return m_phi.entropy_R() * R_water;
}

CanteraDouble WaterPropsIAPWS::cv_mass() const
{
    return m_phi.cv_R() * R_water;
}

CanteraDouble WaterPropsIAPWS::cp_mass() const
{
    return m_phi.cp_R() * R_water;
}

}
