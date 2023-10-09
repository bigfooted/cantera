/**
 *  @file MetalPhase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_METALPHASE_H
#define CT_METALPHASE_H

#include "ThermoPhase.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 *
 * Class MetalPhase represents electrons in a metal.
 */
class MetalPhase : public ThermoPhase
{
public:
    MetalPhase() {}

    // Overloaded methods of class ThermoPhase

    string type() const override {
        return "electron-cloud";
    }

    bool isCompressible() const override {
        return false;
    }

    CanteraDouble enthalpy_mole() const override {
        return 0.0;
    }
    CanteraDouble intEnergy_mole() const override {
        return - pressure() * molarVolume();
    }
    CanteraDouble entropy_mole() const override {
        return 0.0;
    }
    CanteraDouble gibbs_mole() const override {
        return 0.0;
    }
    CanteraDouble cp_mole() const override {
        return 0.0;
    }
    CanteraDouble cv_mole() const override {
        return 0.0;
    }

    void setPressure(CanteraDouble pres) override {
        m_press = pres;
    }
    CanteraDouble pressure() const override {
        return m_press;
    }

    void getChemPotentials(CanteraDouble* mu) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu[n] = 0.0;
        }
    }

    void getEnthalpy_RT(CanteraDouble* hrt) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            hrt[n] = 0.0;
        }
    }

    void getEntropy_R(CanteraDouble* sr) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            sr[n] = 0.0;
        }
    }

    void getStandardChemPotentials(CanteraDouble* mu0) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            mu0[n] = 0.0;
        }
    }

    void getActivityConcentrations(CanteraDouble* c) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            c[n] = 1.0;
        }
    }
    void getPartialMolarEnthalpies(CanteraDouble *h) const override {
        for (size_t n = 0; n < nSpecies(); n++) {
            h[n] = 0.0;
        }
    }

    Units standardConcentrationUnits() const override {
        return Units(1.0);
    }

    CanteraDouble standardConcentration(size_t k=0) const override {
        return 1.0;
    }

    CanteraDouble logStandardConc(size_t k=0) const override {
        return 0.0;
    }

    void initThermo() override {
        if (m_input.hasKey("density")) {
            assignDensity(m_input.convert("density", "kg/m^3"));
        }
    }

    void getParameters(AnyMap& phaseNode) const override {
        ThermoPhase::getParameters(phaseNode);
        phaseNode["density"].setQuantity(density(), "kg/m^3");
    }

private:
    CanteraDouble m_press;
};
}

#endif
