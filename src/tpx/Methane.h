//! @file Methane.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_METHANE_H
#define TPX_METHANE_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of methane. Values and functions are
//! from Reynolds @cite reynolds1979.
class methane : public Substance
{
public:
    methane() {
        m_name = "methane";
        m_formula = "CH4";
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

    //! Saturation pressure. Equation S3 from Reynolds TPSI.
    CanteraDouble Psat() override;

protected:
    //! Liquid density. Equation D3 from Reynolds TPSI.
    CanteraDouble ldens() override;

private:
    CanteraDouble C(int i, CanteraDouble rt, CanteraDouble rt2);
    CanteraDouble Cprime(int i, CanteraDouble rt, CanteraDouble rt2, CanteraDouble rt3);
    CanteraDouble I(int i, CanteraDouble egrho);
    CanteraDouble H(int i, CanteraDouble egrho);
    CanteraDouble W(int i, CanteraDouble egrho);
};
}
#endif // ! METHANE_H
