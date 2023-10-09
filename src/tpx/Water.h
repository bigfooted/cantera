//! @file Water.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_WATER_H
#define TPX_WATER_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! Pure species representation of water. Values and functions are from
//! from Reynolds @cite reynolds1979.
class water : public Substance
{
public:
    water() {
        m_name = "water";
        m_formula = "H2O";
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
    CanteraDouble Psat() override;
    CanteraDouble dPsatdT();

protected:
    CanteraDouble ldens() override;

private:
    CanteraDouble C(int i);
    CanteraDouble Cprime(int i);
    CanteraDouble I(int i);
    CanteraDouble H(int i);
};

}
#endif // ! WATER_H
