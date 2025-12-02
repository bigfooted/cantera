//! @file Hydrogen.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_HYDROGEN_H
#define TPX_HYDROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of hydrogen. Values and functions are
//! from Reynolds @cite reynolds1979.
class hydrogen : public Substance
{
public:
    hydrogen() {
        m_name = "hydrogen";
        m_formula = "H2";
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

    //! Saturation pressure. Equation s3 in Reynolds TPSI.
    CanteraDouble Psat() override;

protected:
    //! Liquid density. Equation D4 in Reynolds TPSI.
    CanteraDouble ldens() override;

private:
    CanteraDouble C(int i, CanteraDouble rt, CanteraDouble rt2);
    CanteraDouble Cprime(int i, CanteraDouble rt, CanteraDouble rt2, CanteraDouble rt3);
    CanteraDouble I(int i, CanteraDouble egrho);
    CanteraDouble H(int i, CanteraDouble egrho);
    CanteraDouble W(int i, CanteraDouble egrho);
    CanteraDouble icv(int i, CanteraDouble x, CanteraDouble xlg);
};

}

#endif // ! HYDROGEN_H
