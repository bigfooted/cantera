//! @file Sub.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_SUB_H
#define TPX_SUB_H

#include "cantera/base/ctexceptions.h"

namespace tpx
{

namespace PropertyPair
{
enum type {
    TV = 12, HP = 34, SP = 54, PV = 42, TP = 14, UV = 62, ST = 51,
    SV = 52, UP = 64, VH = 23, TH = 13, SH = 53, PX = 47, TX = 17,
    VT = -12, PH = -34, PS = -54, VP = -42, PT = -14, VU = -62, TS = -51,
    VS = -52, PU = -64, HV = -23, HT = -13, HS = -53, XP = -47, XT = -17
};
}

const int Pgiven = 0, Tgiven = 1;

namespace propertyFlag
{
enum type { H, S, U, V, P, T };
}

const CanteraDouble Undef = 999.1234;

/**
 * Base class from which all pure substances are derived
 */
class Substance
{
public:
    Substance() = default;

    virtual ~Substance() = default;

    void setStdState(CanteraDouble h0 = 0.0, CanteraDouble s0 = 0.0,
                     CanteraDouble t0 = 298.15, CanteraDouble p0 = 1.01325e5);

    //! @name Information about a substance
    //! @{

    //! Molecular weight [kg/kmol]
    virtual CanteraDouble MolWt()=0;

    //! Critical temperature [K]
    virtual CanteraDouble Tcrit()=0;

    //! Critical pressure [Pa]
    virtual CanteraDouble Pcrit()=0;

    //! Critical specific volume [m^3/kg]
    virtual CanteraDouble Vcrit()=0;

    //! Minimum temperature for which the equation of state is valid
    virtual CanteraDouble Tmin()=0;

    //! Maximum temperature for which the equation of state is valid
    virtual CanteraDouble Tmax()=0;

    //! Name of the substance
    const char* name() {
        return m_name.c_str();
    }

    //! Chemical formula for the substance
     const char* formula() {
        return m_formula.c_str();
    }

    //! @}

    //! @name Properties
    //! @{

    //! Pressure [Pa]. If two phases are present, return the saturation
    //! pressure; otherwise return the pressure computed directly from the
    //! underlying eos.
    CanteraDouble P();

    //! Temperature [K]
    CanteraDouble Temp() {
        return T;
    }

    //! Specific volume [m^3/kg]
    CanteraDouble v() {
        return prop(propertyFlag::V);
    }

    //! Internal energy [J/kg]
    CanteraDouble u() {
        return prop(propertyFlag::U);
    }

    //! Enthalpy [J/kg]
    CanteraDouble h() {
        return prop(propertyFlag::H);
    }

    //! Entropy [J/kg/K]
    CanteraDouble s() {
        return prop(propertyFlag::S);
    }

    //! Helmholtz function [J/kg]
    CanteraDouble f() {
        return u() - T*s();
    }

    //! Gibbs function [J/kg]
    CanteraDouble g() {
        return h() - T*s();
    }

    //! Specific heat at constant volume [J/kg/K]
    virtual CanteraDouble cv();

    //! Specific heat at constant pressure [J/kg/K]
    virtual CanteraDouble cp();

    virtual CanteraDouble thermalExpansionCoeff();

    virtual CanteraDouble isothermalCompressibility();

    //! @}
    //! @name Saturation Properties
    //! @{

    CanteraDouble Ps();

    //! The derivative of the saturation pressure with respect to temperature.
    virtual CanteraDouble dPsdT();

    //! Saturation temperature at pressure *p*.
    CanteraDouble Tsat(CanteraDouble p);

    //! Vapor mass fraction. If T >= Tcrit, 0 is returned for v < Vcrit, and 1
    //! is returned if v > Vcrit.
    CanteraDouble x();

    //! Returns 1 if the current state is a liquid/vapor mixture, 0 otherwise.
    //! By default, saturated vapor and saturated liquid are included; setting
    //! the flag *strict* to true will exclude the boundaries.
    int TwoPhase(bool strict=false);
    //! @}

    virtual CanteraDouble Pp()=0;

    //! Enthalpy of a single-phase state
    CanteraDouble hp() {
        return up() + Pp()/Rho;
    }

    //! Gibbs function of a single-phase state
    CanteraDouble gp() {
        return hp() - T*sp();
    }

    CanteraDouble prop(propertyFlag::type ijob);

    //! set T and P
    void set_TPp(CanteraDouble t0, CanteraDouble p0);

    //! Function to set or change the state for a property pair *XY* where
    //! *x0* is the value of first property and *y0* is the value of the
    //! second property.
    void Set(PropertyPair::type XY, CanteraDouble x0, CanteraDouble y0);

protected:
    CanteraDouble T = Undef;
    CanteraDouble Rho = Undef;
    CanteraDouble Tslast = Undef;
    CanteraDouble Rhf = Undef;
    CanteraDouble Rhv = Undef;
    CanteraDouble Pst = Undef;
    CanteraDouble m_energy_offset = 0.0;
    CanteraDouble m_entropy_offset = 0.0;
    std::string m_name;
    std::string m_formula;

    virtual CanteraDouble ldens()=0;

    //! Saturation pressure, Pa
    virtual CanteraDouble Psat()=0;

    //! Internal energy of a single-phase state
    virtual CanteraDouble up()=0;

    //! Entropy of a single-phase state
    virtual CanteraDouble sp()=0;

    virtual int ideal() {
        return 0;
    }

    CanteraDouble vp() {
        return 1.0/Rho;
    }

    //! Uses the lever rule to set state in the dome. Returns 1 if in dome,
    //! 0 if not, in which case state not set.
    int Lever(int itp, CanteraDouble sat, CanteraDouble val, propertyFlag::type ifunc);

    //! Update saturated liquid and vapor densities and saturation pressure
    void update_sat();

private:
    void set_Rho(CanteraDouble r0);
    void set_T(CanteraDouble t0);
    void set_v(CanteraDouble v0);
    void BracketSlope(CanteraDouble p);
    CanteraDouble vprop(propertyFlag::type ijob);
    void set_xy(propertyFlag::type if1, propertyFlag::type if2,
                CanteraDouble X, CanteraDouble Y,
                CanteraDouble atx, CanteraDouble aty, CanteraDouble rtx, CanteraDouble rty);

    int kbr = 0;
    CanteraDouble Vmin, Vmax;
    CanteraDouble Pmin, Pmax;
    CanteraDouble dvbf, dv;
    CanteraDouble v_here, P_here;
};

}

#endif
