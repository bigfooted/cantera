//! @file Nitrogen.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_NITROGEN_H
#define TPX_NITROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of nitrogen. Values and functions are
//! from Reynolds @cite reynolds1979.
class nitrogen : public Substance
{
public:
    nitrogen() {
        m_name = "nitrogen";
        m_formula = "N2";
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

#endif // ! TPX_NITROGEN_H
