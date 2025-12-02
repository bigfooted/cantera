/**
 *  @file WaterProps.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/WaterProps.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
WaterProps::WaterProps():
    m_waterIAPWS(new WaterPropsIAPWS()),
    m_own_sub(true)
{
}

WaterProps::WaterProps(PDSS_Water* wptr)
{
    if (wptr) {
        // object in slave mode; it doesn't own its own water evaluator.
        m_waterIAPWS = wptr->getWater();
    } else {
        m_waterIAPWS = new WaterPropsIAPWS();
        m_own_sub = true;
    }
}

WaterProps::WaterProps(WaterPropsIAPWS* waterIAPWS)
{
    if (waterIAPWS) {
        m_waterIAPWS = waterIAPWS;
    } else {
        m_waterIAPWS = new WaterPropsIAPWS();
        m_own_sub = true;
    }
}

WaterProps::~WaterProps()
{
    if (m_own_sub) {
        delete m_waterIAPWS;
    }
}

CanteraDouble WaterProps::density_T(CanteraDouble T, CanteraDouble P, int ifunc)
{
    static const CanteraDouble Tc = T - 273.15;
    static const CanteraDouble U1 = 288.9414;
    static const CanteraDouble U2 = 508929.2;
    static const CanteraDouble U3 = 68.12963;
    static const CanteraDouble U4 = -3.9863;

    CanteraDouble tmp1 = Tc + U1;
    CanteraDouble tmp4 = Tc + U4;
    CanteraDouble t4t4 = tmp4 * tmp4;
    CanteraDouble tmp3 = Tc + U3;
    CanteraDouble rho = 1000. * (1.0 - tmp1*t4t4/(U2 * tmp3));

    // Impose an ideal gas lower bound on rho. We need this to ensure positivity
    // of rho, even though it is grossly unrepresentative.
    CanteraDouble rhomin = P / (GasConstant * T);
    if (rho < rhomin) {
        rho = rhomin;
        if (ifunc == 1) {
            return - rhomin / T;
        } else if (ifunc == 3) {
            return rhomin / P;
        } else if (ifunc == 2) {
            return 2.0 * rhomin / (T * T);
        }
    }

    if (ifunc == 1) {
        CanteraDouble drhodT = 1000./U2 * (
                                - tmp4 * tmp4 / (tmp3)
                                - tmp1 * 2 * tmp4 / (tmp3)
                                + tmp1 * t4t4 / (tmp3*tmp3)
                            );
        return drhodT;
    } else if (ifunc == 3) {
        return 0.0;
    } else if (ifunc == 2) {
        CanteraDouble t3t3 = tmp3 * tmp3;
        CanteraDouble d2rhodT2 = 1000./U2 *
                               ((-4.0*tmp4-2.0*tmp1)/tmp3 +
                                (2.0*t4t4 + 4.0*tmp1*tmp4)/t3t3
                                - 2.0*tmp1 * t4t4/(t3t3*tmp3));
        return d2rhodT2;
    }
    return rho;
}

CanteraDouble WaterProps::relEpsilon(CanteraDouble T, CanteraDouble P_pascal, int ifunc)
{
    static const CanteraDouble U1 = 3.4279E2;
    static const CanteraDouble U2 = -5.0866E-3;
    static const CanteraDouble U3 = 9.4690E-7;
    static const CanteraDouble U4 = -2.0525;
    static const CanteraDouble U5 = 3.1159E3;
    static const CanteraDouble U6 = -1.8289E2;
    static const CanteraDouble U7 = -8.0325E3;
    static const CanteraDouble U8 = 4.2142E6;
    static const CanteraDouble U9 = 2.1417;
    CanteraDouble T2 = T * T;

    CanteraDouble eps1000 = U1 * exp(U2 * T + U3 * T2);
    CanteraDouble C = U4 + U5/(U6 + T);
    CanteraDouble B = U7 + U8/T + U9 * T;
    CanteraDouble Pbar = P_pascal * 1.0E-5;
    CanteraDouble tmpBpar = B + Pbar;
    CanteraDouble tmpB1000 = B + 1000.0;
    CanteraDouble ltmp = log(tmpBpar/tmpB1000);
    CanteraDouble epsRel = eps1000 + C * ltmp;

    if (ifunc == 1 || ifunc == 2) {
        CanteraDouble tmpC = U6 + T;
        CanteraDouble dCdT = - U5/(tmpC * tmpC);
        CanteraDouble dBdT = - U8/(T * T) + U9;
        CanteraDouble deps1000dT = eps1000 * (U2 + 2.0 * U3 * T);
        CanteraDouble dltmpdT = (dBdT/tmpBpar - dBdT/tmpB1000);
        if (ifunc == 1) {
            return deps1000dT + dCdT * ltmp + C * dltmpdT;
        }
        CanteraDouble T3 = T2 * T;
        CanteraDouble d2CdT2 = - 2.0 * dCdT / tmpC;
        CanteraDouble d2BdT2 = 2.0 * U8 / (T3);
        CanteraDouble d2ltmpdT2 = (d2BdT2*(1.0/tmpBpar - 1.0/tmpB1000) +
                                dBdT*dBdT*(1.0/(tmpB1000*tmpB1000) - 1.0/(tmpBpar*tmpBpar)));
        CanteraDouble d2eps1000dT2 = (deps1000dT * (U2 + 2.0 * U3 * T) + eps1000 * (2.0 * U3));

        if (ifunc == 2) {
            CanteraDouble d2epsReldT2 = (d2eps1000dT2 + d2CdT2 * ltmp + 2.0 * dCdT * dltmpdT
                                      + C * d2ltmpdT2);
            return d2epsReldT2;
        }
    }
    if (ifunc == 3) {
        CanteraDouble dltmpdP = 1.0E-5 / tmpBpar;
        return C * dltmpdP;
    }
    return epsRel;
}

CanteraDouble WaterProps::ADebye(CanteraDouble T, CanteraDouble P_input, int ifunc)
{
    CanteraDouble psat = satPressure(T);
    CanteraDouble P;
    if (psat > P_input) {
        P = psat;
    } else {
        P = P_input;
    }
    CanteraDouble epsRelWater = relEpsilon(T, P, 0);
    CanteraDouble epsilon = epsilon_0 * epsRelWater;
    CanteraDouble dw = density_IAPWS(T, P);
    CanteraDouble tmp = sqrt(2.0 * Avogadro * dw / 1000.);
    CanteraDouble tmp2 = ElectronCharge * ElectronCharge * Avogadro /
                      (epsilon * GasConstant * T);
    CanteraDouble tmp3 = tmp2 * sqrt(tmp2);
    CanteraDouble A_Debye = tmp * tmp3 / (8.0 * Pi);

    // dAdT = - 3/2 Ad/T + 1/2 Ad/dw d(dw)/dT - 3/2 Ad/eps d(eps)/dT
    // dAdT = - 3/2 Ad/T - 1/2 Ad/Vw d(Vw)/dT - 3/2 Ad/eps d(eps)/dT
    if (ifunc == 1 || ifunc == 2) {
        CanteraDouble dAdT = - 1.5 * A_Debye / T;

        CanteraDouble depsRelWaterdT = relEpsilon(T, P, 1);
        dAdT -= A_Debye * (1.5 * depsRelWaterdT / epsRelWater);

        // calculate d(lnV)/dT _constantP, that is, the cte
        CanteraDouble cte = coeffThermalExp_IAPWS(T, P);
        CanteraDouble contrib2 = - A_Debye * (0.5 * cte);
        dAdT += contrib2;

        if (ifunc == 1) {
            return dAdT;
        }

        if (ifunc == 2) {
            // Get the second derivative of the dielectric constant wrt T
            // -> we will take each of the terms in dAdT and differentiate
            //    it again.
            CanteraDouble d2AdT2 = 1.5 / T * (A_Debye/T - dAdT);
            CanteraDouble d2epsRelWaterdT2 = relEpsilon(T, P, 2);
            d2AdT2 += 1.5 * (- dAdT * depsRelWaterdT / epsRelWater
                             - A_Debye / epsRelWater *
                             (d2epsRelWaterdT2 - depsRelWaterdT * depsRelWaterdT / epsRelWater));
            CanteraDouble deltaT = -0.1;
            CanteraDouble Tdel = T + deltaT;
            CanteraDouble cte_del = coeffThermalExp_IAPWS(Tdel, P);
            CanteraDouble dctedT = (cte_del - cte) / Tdel;
            CanteraDouble contrib3 = 0.5 * (-(dAdT * cte) -(A_Debye * dctedT));
            d2AdT2 += contrib3;
            return d2AdT2;
        }
    }

    // A_Debye = (1/(8 Pi)) sqrt(2 Na dw / 1000)
    //                          (e e/(epsilon R T))^3/2
    //
    // dAdP =  + 1/2 Ad/dw d(dw)/dP - 3/2 Ad/eps d(eps)/dP
    // dAdP =  - 1/2 Ad/Vw d(Vw)/dP - 3/2 Ad/eps d(eps)/dP
    // dAdP =  + 1/2 Ad * kappa  - 3/2 Ad/eps d(eps)/dP
    //
    // where kappa = - 1/Vw d(Vw)/dP_T (isothermal compressibility)
    if (ifunc == 3) {
        CanteraDouble dAdP = 0.0;
        CanteraDouble depsRelWaterdP = relEpsilon(T, P, 3);
        dAdP -= A_Debye * (1.5 * depsRelWaterdP / epsRelWater);
        CanteraDouble kappa = isothermalCompressibility_IAPWS(T,P);
        dAdP += A_Debye * (0.5 * kappa);
        return dAdP;
    }
    return A_Debye;
}

CanteraDouble WaterProps::satPressure(CanteraDouble T)
{
    return m_waterIAPWS->psat(T);
}

CanteraDouble WaterProps::density_IAPWS(CanteraDouble temp, CanteraDouble press)
{
    return m_waterIAPWS->density(temp, press, WATER_LIQUID);
}

CanteraDouble WaterProps::density_IAPWS() const
{
    return m_waterIAPWS->density();
}

CanteraDouble WaterProps::coeffThermalExp_IAPWS(CanteraDouble temp, CanteraDouble press)
{
    CanteraDouble dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
        throw CanteraError("WaterProps::coeffThermalExp_IAPWS",
            "Unable to solve for density at T = {} and P = {}", temp, press);
    }
    return m_waterIAPWS->coeffThermExp();
}

CanteraDouble WaterProps::isothermalCompressibility_IAPWS(CanteraDouble temp, CanteraDouble press)
{
    CanteraDouble dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
        throw CanteraError("WaterProps::isothermalCompressibility_IAPWS",
            "Unable to solve for density at T = {} and P = {}", temp, press);
    }
    return m_waterIAPWS->isothermalCompressibility();
}

}
