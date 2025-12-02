/**
 * @file CarbonDioxide.cpp representation of substance Carbon Dioxide.
 *
 * Values and functions are from "Thermodynamic Properties in SI" by W.C.
 * Reynolds AUTHOR: me@rebeccahhunt.com: GCEP, Stanford University
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "CarbonDioxide.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{

/*
 * Carbon Dioxide constants
 */
static const CanteraDouble Tmn = 216.54; // [K] minimum temperature for which calculations are valid
static const CanteraDouble Tmx = 1500.0; // [K] maximum temperature for which calculations are valid
static const CanteraDouble Tc=304.21; // [K] critical temperature
static const CanteraDouble Roc=464.00; // [kg/m^3] critical density
static const CanteraDouble To=216.54; // [K] reference Temperature
static const CanteraDouble R=188.918; // [] gas constant for CO2 J/kg/K
static const CanteraDouble Gamma=5.0E-6; // [??]
static const CanteraDouble u0=3.2174105E5; // [] internal energy at To
static const CanteraDouble s0=2.1396056E3; // [] entropy at To
static const CanteraDouble Tp=250; // [K] ??
static const CanteraDouble Pc=7.38350E6; // [Pa] critical pressure
static const CanteraDouble M=44.01; // [kg/kmol] molar density

// array Acarbdi is used by the function named Pp
static const CanteraDouble Acarbdi[]= {
    2.2488558E-1,
    -1.3717965E2,
    -1.4430214E4,
    -2.9630491E6,
    -2.0606039E8,
    4.5554393E-5,
    7.7042840E-2,
    4.0602371E1,
    4.0029509E-7,
    -3.9436077E-4,
    1.2115286E-10,
    1.0783386E-7,
    4.3962336E-11,
    -3.6505545E4,
    1.9490511E7,
    -2.9186718E9,
    2.4358627E-2,
    -3.7546530E1,
    1.1898141E4
};

// array F is used by the function named Psat
static const CanteraDouble F[]= {
    -6.5412610,
    -2.7914636E-1,
    -3.4716202,
    -3.4989637,
    -1.9770948E1,
    1.3922839E2,
    -2.7670389E2,
    -7.0510251E3
};

// array D is used by the function ldens
static const CanteraDouble D[]= {
    4.6400009E2,
    6.7938129E2,
    1.4776836E3,
    -3.1267676E3,
    3.6397656E3,
    -1.3437098E3
};

// array G is used by the function sp
static const CanteraDouble G[]= {
    8.726361E3,
    1.840040E2,
    1.914025,
    -1.667825E-3,
    7.305950E-7,
    -1.255290E-10,
};

CanteraDouble CarbonDioxide::C(int j,CanteraDouble Tinverse, CanteraDouble T2inverse, CanteraDouble T3inverse, CanteraDouble T4inverse)
{
    switch (j) {
    case 0:
        return Acarbdi[0]*T +
               Acarbdi[1] +
               Acarbdi[2] * Tinverse +
               Acarbdi[3] * T2inverse +
               Acarbdi[4] * T3inverse;
    case 1:
        return Acarbdi[5] *T +
               Acarbdi[6] +
               Acarbdi[7] * Tinverse;
    case 2:
        return Acarbdi[8]*T + Acarbdi[9];
    case 3:
        return Acarbdi[10]*T + Acarbdi[11];
    case 4:
        return Acarbdi[12];
    case 5:
        return Acarbdi[13] *T2inverse +
               Acarbdi[14] *T3inverse +
               Acarbdi[15] *T4inverse;
    case 6:
        return Acarbdi[16] *T2inverse +
               Acarbdi[17] *T3inverse +
               Acarbdi[18] *T4inverse;
    default:
        return 0.0;
    }
}

CanteraDouble CarbonDioxide::Cprime(int j, CanteraDouble T2inverse, CanteraDouble T3inverse, CanteraDouble T4inverse)
{
    switch (j) {
    case 0:
        return Acarbdi[0] +
               - Acarbdi[2] * T2inverse +
               -2 * Acarbdi[3] * T3inverse +
               -3 * Acarbdi[4] * T4inverse;
    case 1:
        return Acarbdi[5] -
               Acarbdi[7] * T2inverse;
    case 2:
        return Acarbdi[8];
    case 3:
        return Acarbdi[10];
    case 4:
        return 0;
    case 5:
        return
            -2 *Acarbdi[13] *T3inverse +
            -3 *Acarbdi[14] *T4inverse +
            -4 *Acarbdi[15]* pow(T,-5);
    case 6:
        return
            -2 *Acarbdi[16] *T3inverse +
            -3 *Acarbdi[17] *T4inverse +
            -4 *Acarbdi[18] *pow(T,-5);
    default:
        return 0.0;
    }
}

CanteraDouble CarbonDioxide::I(int j, CanteraDouble ergho, CanteraDouble Gamma)
{
    switch (j) {
    case 0:
        return Rho;
    case 1:
        return pow(Rho, 2)/2;
    case 2:
        return pow(Rho, 3)/ 3;
    case 3:
        return pow(Rho, 4)/ 4;
    case 4:
        return pow(Rho, 5)/ 5;
    case 5:
        return (1 - ergho) / CanteraDouble(2 * Gamma);
    case 6:
        return (1 - ergho * CanteraDouble(Gamma * pow(Rho,2) + CanteraDouble(1)))/ CanteraDouble(2 * Gamma * Gamma);
    default:
        return 0.0;
    }
}

CanteraDouble CarbonDioxide::H(int i, CanteraDouble egrho)
{
    if (i < 5) {
        return pow(Rho,i+2);
    } else if (i == 5) {
        return pow(Rho,3)*egrho;
    } else if (i == 6) {
        return pow(Rho,5)*egrho;
    } else {
        return 0;
    }
}

CanteraDouble CarbonDioxide::up()
{
    CanteraDouble Tinverse = 1.0/T;
    CanteraDouble T2inverse = pow(T, -2);
    CanteraDouble T3inverse = pow(T, -3);
    CanteraDouble T4inverse = pow(T, -4);
    CanteraDouble egrho = exp(-Gamma*Rho*Rho);

    CanteraDouble sum = 0.0;
    // Equation C-6 integrated
    sum += G[0]*log(T/To);
    int i;
    for (i=1; i<=5; i++) {
        sum += G[i]*(pow(T,i) - pow(To,i))/CanteraDouble(i);
    }
    for (i=0; i<=6; i++) {
        sum += I(i,egrho, Gamma) *
               (C(i, Tinverse, T2inverse, T3inverse, T4inverse) - T*Cprime(i,T2inverse, T3inverse, T4inverse));
    }
    sum += u0;
    return sum + m_energy_offset;
}

CanteraDouble CarbonDioxide::sp()
{
    CanteraDouble T2inverse = pow(T, -2);
    CanteraDouble T3inverse = pow(T, -3);
    CanteraDouble T4inverse = pow(T, -4);
    CanteraDouble egrho = exp(-Gamma*Rho*Rho);

    CanteraDouble sum = 0.0;
    for (int i=2; i<=5; i++) {
        sum += G[i]*(pow(T,i-1) - pow(To,i-1))/CanteraDouble(i-1);
    }
    sum += G[1]*log(T/To);
    sum -= G[0]*(1.0/T - 1.0/To);
    for (int i=0; i<=6; i++) {
        sum -= Cprime(i,T2inverse, T3inverse, T4inverse)*I(i,egrho,Gamma);
    }
    sum += s0 - R*log(Rho);
    return sum + m_entropy_offset;
}

CanteraDouble CarbonDioxide::Pp()
{
    CanteraDouble Tinverse = pow(T,-1);
    CanteraDouble T2inverse = pow(T, -2);
    CanteraDouble T3inverse = pow(T, -3);
    CanteraDouble T4inverse = pow(T, -4);
    CanteraDouble egrho = exp(-Gamma*Rho*Rho);
    CanteraDouble P = Rho*R*T;

    // when i=0 we are on second sum of equation (where rho^2)
    for (int i=0; i<=6; i++) {
        P += C(i,Tinverse, T2inverse, T3inverse, T4inverse)*H(i,egrho);
    }
    return P;
}

CanteraDouble CarbonDioxide::Psat()
{
    CanteraDouble log, sum=0,P;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("CarbonDixoide::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=8; i++) {
        sum += F[i-1] * pow((T/Tp -1),CanteraDouble(i-1));
    }

    log = ((Tc/T)-1)*sum;
    P=exp(log)*Pc;
    return P;
}

CanteraDouble CarbonDioxide::ldens()
{
    CanteraDouble xx=1-(T/Tc), sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("CarbonDixoide::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=6; i++) {
        sum+=D[i-1]*pow(xx,CanteraDouble(i-1)/3.0);
    }
    return sum;
}

// The following functions allow users to get the properties of CarbonDioxide
// that are not dependent on the state

CanteraDouble CarbonDioxide::Tcrit()
{
    return Tc;
}
CanteraDouble CarbonDioxide::Pcrit()
{
    return Pc;
}
CanteraDouble CarbonDioxide::Vcrit()
{
    return 1.0/Roc;
}
CanteraDouble CarbonDioxide::Tmin()
{
    return Tmn;
}
CanteraDouble CarbonDioxide::Tmax()
{
    return Tmx;
}
CanteraDouble CarbonDioxide::MolWt()
{
    return M;
}

}
