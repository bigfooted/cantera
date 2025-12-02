//! @file IonFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IONFLOW_H
#define CT_IONFLOW_H

#include "cantera/oneD/Flow1D.h"

namespace Cantera
{
/**
 * This class models the ion transportation in a flame. There are three
 * stages of the simulation.
 *
 * The first stage turns off the diffusion of ions due to the fast
 * diffusion rate of electron without internal electric forces (ambi-
 * polar diffusion effect).
 *
 * The second stage evaluates drift flux from electric field calculated from
 * Poisson's equation, which is solved together with other equations. Poisson's
 * equation is coupled because the total charge densities depends on the species'
 * concentration. See Pedersen and Brown @cite pedersen1993 for details.
 *
 * @ingroup flowGroup
 */
class IonFlow : public Flow1D
{
public:
    //! Create a new IonFlow domain.
    //! @param phase  Solution object used to evaluate all thermodynamic, kinetic, and
    //!     transport properties
    //! @param id  name of flow domain
    //! @param points  initial number of grid points
    IonFlow(shared_ptr<Solution> phase, const string& id="", size_t points = 1);

    string domainType() const override;

    void resize(size_t components, size_t points) override;
    bool componentActive(size_t n) const override;

    void solveElectricField() override;
    void fixElectricField() override;
    bool doElectricField() const override {
        return m_do_electric_field;
    }

    /**
     * Sometimes it is desired to carry out the simulation using a specified
     * electron transport profile, rather than assuming it as a constant (0.4).
     * See Bisetti and El Morsli @cite bisetti2012.
     * If in the future the class GasTransport is improved, this method may
     * be discarded. This method specifies this profile.
     */
    void setElectronTransport(vector<CanteraDouble>& tfix,
                              vector<CanteraDouble>& diff_e,
                              vector<CanteraDouble>& mobi_e);

protected:

    /**
     * Evaluate the electric field equation residual by Gauss's law.
     *
     * The function calculates the electric field equation as:
     * @f[
     *    \frac{dE}{dz} = \frac{e}{\varepsilon_0} \sum (q_k \cdot n_k)
     * @f]
     *
     * and
     *
     * @f[
     *    E = -\frac{dV}{dz}
     * @f]
     *
     * The electric field equation is based on Gauss's law,
     * accounting for charge density and permittivity of free space
     * (@f$ \varepsilon_0 @f$).
     * The zero electric field is first evaluated and if the solution state is 2,
     * then the alternative form the electric field equation is evaluated.
     *
     * For argument explanation, see evalContinuity() base class.
     */
    void evalElectricField(CanteraDouble* x, CanteraDouble* rsd, int* diag,
                           CanteraDouble rdt, size_t jmin, size_t jmax) override;

    /**
     * Evaluate the species equations' residual. This function overloads the
     * original species function.
     *
     * A Neumann boundary for the charged species at the
     * left boundary is added, and the default boundary condition from the overloaded
     * method is left the same for the right boundary.
     *
     * For argument explanation, see evalContinuity() base class.
     */
    void evalSpecies(CanteraDouble* x, CanteraDouble* rsd, int* diag,
                     CanteraDouble rdt, size_t jmin, size_t jmax) override;
    void updateTransport(CanteraDouble* x, size_t j0, size_t j1) override;
    void updateDiffFluxes(const CanteraDouble* x, size_t j0, size_t j1) override;
    //! Solving phase one: the fluxes of charged species are turned off and the electric
    //! field is not solved.
    void frozenIonMethod(const CanteraDouble* x, size_t j0, size_t j1);
    //! Solving phase two: the electric field equation is added coupled
    //! by the electrical drift
    void electricFieldMethod(const CanteraDouble* x, size_t j0, size_t j1);
    //! flag for solving electric field or not
    bool m_do_electric_field = false;

    //! flag for importing transport of electron
    bool m_import_electron_transport = false;

    //! electrical properties
    vector<CanteraDouble> m_speciesCharge;

    //! index of species with charges
    vector<size_t> m_kCharge;

    //! index of neutral species
    vector<size_t> m_kNeutral;

    //! Coefficients of polynomial fit for electron mobility as a function of
    //! temperature.
    //! @see setElectronTransport
    vector<CanteraDouble> m_mobi_e_fix;

    //! Coefficients of polynomial fit for electron diffusivity as a function of
    //! temperature.
    //! @see setElectronTransport
    vector<CanteraDouble> m_diff_e_fix;

    //! mobility
    vector<CanteraDouble> m_mobility;

    //! index of electron
    size_t m_kElectron = npos;

    //! electric field [V/m]
    CanteraDouble E(const CanteraDouble* x, size_t j) const {
        return x[index(c_offset_E, j)];
    }

    //! Axial gradient of the electric field [V/m²]
    CanteraDouble dEdz(const CanteraDouble* x, size_t j) const {
        return (E(x,j)-E(x,j-1))/(z(j)-z(j-1));
    }

    //! number density [molecules/m³]
    CanteraDouble ND(const CanteraDouble* x, size_t k, size_t j) const {
        return Avogadro * m_rho[j] * Y(x,k,j) / m_wt[k];
    }

    //! total charge density
    CanteraDouble rho_e(CanteraDouble* x, size_t j) const {
        CanteraDouble chargeDensity = 0.0;
        for (size_t k : m_kCharge) {
            chargeDensity += m_speciesCharge[k] * ElectronCharge * ND(x,k,j);
        }
        return chargeDensity;
    }
};

}

#endif
