//! @file Oxygen.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_OXYGEN_H
#define TPX_OXYGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! Pure species representation of oxygen. Values and functions are
//! from Reynolds @cite reynolds1979.
class oxygen : public Substance
{
public:
    oxygen() {
        m_name="oxygen";
        m_formula="O2";
    }

    CanteraDouble MolWt() override;
    CanteraDouble Tcrit() override;
    CanteraDouble Pcrit() override;
    CanteraDouble Vcrit() override;
    CanteraDouble Tmin() override;
    CanteraDouble Tmax() override;

    CanteraDouble Pp() override;
    CanteraDouble up() override;
    CanteraDouble sp() override;

    //! Saturation pressure. Equation S4 from Reynolds TPSI.
    CanteraDouble Psat() override;

protected:
    //! Liquid density. Equation D2 from Reynolds TPSI.
    CanteraDouble ldens() override;

private:
    //! Equation P4 from Reynolds TPSI.
    CanteraDouble C(int i, CanteraDouble rt, CanteraDouble rt2);
    CanteraDouble Cprime(int i, CanteraDouble rt, CanteraDouble rt2, CanteraDouble rt3);
    CanteraDouble I(int i, CanteraDouble egrho);
    CanteraDouble H(int i, CanteraDouble egrho);
    CanteraDouble W(int i, CanteraDouble egrho);
};

}
#endif // ! OXYGEN_H
