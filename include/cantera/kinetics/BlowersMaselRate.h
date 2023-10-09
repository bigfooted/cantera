//! @file BlowersMaselRate.h   Header for Blowers-Masel reaction rates

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BLOWERSMASELRATE_H
#define CT_BLOWERSMASELRATE_H

#include "Arrhenius.h"

namespace Cantera
{

//! Data container holding shared data specific to BlowersMaselRate
/**
 * The data container `BlowersMaselData` holds precalculated data common to
 * all `BlowersMaselRate` objects.
 */
struct BlowersMaselData : public ReactionData
{
    BlowersMaselData() = default;

    void update(CanteraDouble T) override;
    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;

    void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        partialMolarEnthalpies.resize(nSpecies, 0.);
        ready = true;
    }

    bool ready = false; //!< boolean indicating whether vectors are accessible
    CanteraDouble density = NAN; //!< used to determine if updates are needed
    vector<CanteraDouble> partialMolarEnthalpies; //!< partial molar enthalpies

protected:
    int m_state_mf_number = -1; //!< integer that is incremented when composition changes
};


//! Blowers Masel reaction rate type depends on the enthalpy of reaction
/**
 * The Blowers Masel approximation @cite blowers2004 adjusts the activation energy
 * based on enthalpy change of a reaction:
 *
 *   @f{eqnarray*}{
 *        E_a &=& 0\; &\text{if }\Delta H < -4E_0 \\
 *        E_a &=& \Delta H\; &\text{if }\Delta H > 4E_0 \\
 *        E_a &=& \frac{(w + \Delta H / 2)(V_P - 2w +
 *               \Delta H)^2}{(V_P^2 - 4w^2 + (\Delta H)^2)}\; &\text{otherwise}
 *   @f}
 * where
 *   @f[
 *        V_P = \frac{2w (w + E_0)}{w - E_0},
 *   @f]
 * @f$ w @f$ is the average bond dissociation energy of the bond breaking
 * and that being formed in the reaction. Since the expression is
 * very insensitive to @f$ w @f$ for @f$ w >= 2 E_0 @f$, @f$ w @f$
 * can be approximated to an arbitrary high value like 1000 kJ/mol.
 *
 * After the activation energy is determined by Blowers-Masel approximation,
 * it can be plugged into Arrhenius function to calculate the rate constant.
 *   @f[
 *        k_f =  A T^b \exp (-E_a/RT)
 *   @f]
 *
 * @ingroup arrheniusGroup
 */
class BlowersMaselRate : public ArrheniusBase
{
public:
    //! Default constructor.
    BlowersMaselRate();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea0  Intrinsic activation energy in energy units [J/kmol]
     *  @param w  Average bond dissociation energy of the bond being formed and
     *      broken in the reaction, in energy units [J/kmol]
     */
    BlowersMaselRate(CanteraDouble A, CanteraDouble b, CanteraDouble Ea0, CanteraDouble w);

    explicit BlowersMaselRate(const AnyMap& node,
                              const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<BlowersMaselRate, BlowersMaselData>>();
    }

    const string type() const override {
        return "Blowers-Masel";
    }

    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate reaction rate
    CanteraDouble evalRate(CanteraDouble logT, CanteraDouble recipT) const {
        CanteraDouble Ea_R = effectiveActivationEnergy_R(m_deltaH_R);
        return m_A * std::exp(m_b * logT - Ea_R * recipT);
    }

    //! Update information specific to reaction
    void updateFromStruct(const BlowersMaselData& shared_data) {
        if (shared_data.ready) {
            m_deltaH_R = 0.;
            for (const auto& [k, multiplier] : m_stoich_coeffs) {
                m_deltaH_R += shared_data.partialMolarEnthalpies[k] * multiplier;
            }
            m_deltaH_R /= GasConstant;
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    CanteraDouble evalFromStruct(const BlowersMaselData& shared_data) const {
        CanteraDouble Ea_R = effectiveActivationEnergy_R(m_deltaH_R);
        return m_A * std::exp(m_b * shared_data.logT - Ea_R * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  This method does not consider potential changes due to a changed reaction
     *  enthalpy. A corresponding warning is raised.
     *  @param shared_data  data shared by all reactions of a given type
     */
    CanteraDouble ddTScaledFromStruct(const BlowersMaselData& shared_data) const;

protected:
    //! Return the effective activation energy (a function of the delta H of reaction)
    //! divided by the gas constant (that is, the activation temperature) [K]
    CanteraDouble effectiveActivationEnergy_R(CanteraDouble deltaH_R) const {
        if (deltaH_R < -4 * m_Ea_R) {
            return 0.;
        }
        if (deltaH_R > 4 * m_Ea_R) {
            return deltaH_R;
        }
        // m_E4_R is the bond dissociation energy "w" (in temperature units)
        CanteraDouble vp = 2 * m_E4_R * ((m_E4_R + m_Ea_R) / (m_E4_R - m_Ea_R)); // in Kelvin
        CanteraDouble vp_2w_dH = (vp - 2 * m_E4_R + deltaH_R); // (Vp - 2 w + dH)
        return (m_E4_R + deltaH_R / 2) * (vp_2w_dH * vp_2w_dH) /
            (vp * vp - 4 * m_E4_R * m_E4_R + deltaH_R * deltaH_R); // in Kelvin
    }

public:
    CanteraDouble activationEnergy() const override {
        return effectiveActivationEnergy_R(m_deltaH_R) * GasConstant;
    }

    //! Return the bond dissociation energy *w* [J/kmol]
    CanteraDouble bondEnergy() const {
        return m_E4_R * GasConstant;
    }

    //! Return current enthalpy change of reaction [J/kmol]
    CanteraDouble deltaH() const {
        return m_deltaH_R * GasConstant;
    }

    //! Set current enthalpy change of reaction [J/kmol]
    /*!
     *  @note  used for testing purposes only; this quantity is not an
     *      independent variable and will be overwritten during an update of the state.
     *
     *  @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    void setDeltaH(CanteraDouble deltaH) {
        m_deltaH_R = deltaH / GasConstant;
    }

protected:
    //! Pairs of species indices and multipliers to calculate enthalpy change
    vector<pair<size_t, CanteraDouble>> m_stoich_coeffs;

    CanteraDouble m_deltaH_R = 0.0; //!< enthalpy change of reaction (in temperature units)
};

}

#endif
