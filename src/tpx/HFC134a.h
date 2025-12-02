//! @file HFC134a.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_HFC134_H
#define TPX_HFC134_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! Equation of state for HFC-134a.
//!
//! Implements the equation of state given in Tillner-Roth and Baehr
//! @cite tillner-roth1994.
class HFC134a : public Substance
{
public:
    HFC134a() {
        m_name = "HFC-134a";
        m_formula = "C2F4H2";
    }

    CanteraDouble MolWt() override;
    CanteraDouble Tcrit() override;
    CanteraDouble Pcrit() override;
    CanteraDouble Vcrit() override;
    CanteraDouble Tmin() override;
    CanteraDouble Tmax() override;

    CanteraDouble Pp() override;
    CanteraDouble fp();
    CanteraDouble up() override;
    CanteraDouble sp() override {
        return ((up() - m_energy_offset) - fp())/T + m_entropy_offset;
    }
    CanteraDouble Psat() override;

protected:
    CanteraDouble ldens() override;
};
}
#endif // ! HFC134_H
