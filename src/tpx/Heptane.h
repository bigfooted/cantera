//! @file Heptane.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_HEPTANE_H
#define TPX_HEPTANE_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! Pure species representation of heptane. Values and functions are
//! from Reynolds @cite reynolds1979.
class Heptane : public Substance
{
public:
    Heptane() {
        m_name = "heptane";
        m_formula = "C7H16";
    }

    CanteraDouble MolWt() override;
    CanteraDouble Tcrit() override;
    CanteraDouble Pcrit() override;
    CanteraDouble Vcrit() override;
    CanteraDouble Tmin() override;
    CanteraDouble Tmax() override;

    //! Pressure. Equation P-2 in Reynolds.
    CanteraDouble Pp() override;

    /**
     * Internal energy.
     * See Reynolds eqn (15) section 2
     *  u = (the integral from T to To of co(T)dT) +
     *         sum from i to N ([C(i) - T*Cprime(i)] + uo
     */
    CanteraDouble up() override;

    //! Entropy. See Reynolds eqn (16) section 2
    CanteraDouble sp() override;

    //! Pressure at Saturation. Equation S-2 in Reynolds.
    CanteraDouble Psat() override;

private:
    //! liquid density. Equation D2 in Reynolds.
    CanteraDouble ldens() override;

    /**
     * C returns a multiplier in each term of the sum
     * in P-2, used in conjunction with C in the function Pp
     * - j is used to represent which of the values in the summation to calculate
     * - j=0 is the second additive in the formula in reynolds
     * - j=1 is the third...
     */
    CanteraDouble C(int jm, CanteraDouble, CanteraDouble, CanteraDouble, CanteraDouble);

    //! derivative of C(i)
    CanteraDouble Cprime(int i, CanteraDouble, CanteraDouble, CanteraDouble);

    /**
     * I = integral from o-rho { 1/(rho^2) * H(i, rho) d rho }
     * ( see section 2 of Reynolds TPSI )
     */
    CanteraDouble I(int i, CanteraDouble, CanteraDouble);

    /**
     * H returns a multiplier in each term of the sum in P-2.
     * this is used in conjunction with C in the function Pp
     * this represents the product rho^n
     * - i=0 is the second additive in the formula in reynolds
     * - i=1 is the third ...
     */
    CanteraDouble H(int i, CanteraDouble egrho);
};

}

#endif // ! TPX_HEPTANE_H
