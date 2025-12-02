/**
 * @file BulkKinetics.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BULKKINETICS_H
#define CT_BULKKINETICS_H

#include "Kinetics.h"
#include "ThirdBodyCalc.h"

namespace Cantera
{

//! Specialization of Kinetics for chemistry in a single bulk phase
//! @ingroup kineticsmgr
class BulkKinetics : public Kinetics
{
public:
    //! @name Constructors and General Information
    //! @{
    BulkKinetics();

    string kineticsType() const override {
        return "bulk";
    }

    //! @}

    //! @name Reaction Mechanism Setup Routines
    //! @{
    bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;
    void addThirdBody(shared_ptr<Reaction> r);
    void modifyReaction(size_t i, shared_ptr<Reaction> rNew) override;
    void resizeSpecies() override;
    void resizeReactions() override;
    void setMultiplier(size_t i, CanteraDouble f) override;
    void invalidateCache() override;
    //! @}

    //! @name Reaction rate constants, rates of progress, and thermodynamic properties
    //! @{
    void getFwdRateConstants(CanteraDouble* kfwd) override;
    void getEquilibriumConstants(CanteraDouble* kc) override;
    void getRevRateConstants(CanteraDouble* krev, bool doIrreversible=false) override;

    void getDeltaGibbs(CanteraDouble* deltaG) override;
    void getDeltaEnthalpy(CanteraDouble* deltaH) override;
    void getDeltaEntropy(CanteraDouble* deltaS) override;

    void getDeltaSSGibbs(CanteraDouble* deltaG) override;
    void getDeltaSSEnthalpy(CanteraDouble* deltaH) override;
    void getDeltaSSEntropy(CanteraDouble* deltaS) override;
    //! @}

    //! @name Derivatives of rate constants and rates of progress
    //! @{
    void getDerivativeSettings(AnyMap& settings) const override;
    void setDerivativeSettings(const AnyMap& settings) override;
    void getFwdRateConstants_ddT(CanteraDouble* dkfwd) override;
    void getFwdRatesOfProgress_ddT(CanteraDouble* drop) override;
    void getRevRatesOfProgress_ddT(CanteraDouble* drop) override;
    void getNetRatesOfProgress_ddT(CanteraDouble* drop) override;
    void getFwdRateConstants_ddP(CanteraDouble* dkfwd) override;
    void getFwdRatesOfProgress_ddP(CanteraDouble* drop) override;
    void getRevRatesOfProgress_ddP(CanteraDouble* drop) override;
    void getNetRatesOfProgress_ddP(CanteraDouble* drop) override;
    void getFwdRateConstants_ddC(CanteraDouble* dkfwd) override;
    void getFwdRatesOfProgress_ddC(CanteraDouble* drop) override;
    void getRevRatesOfProgress_ddC(CanteraDouble* drop) override;
    void getNetRatesOfProgress_ddC(CanteraDouble* drop) override;
    Eigen::SparseMatrix<CanteraDouble> fwdRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<CanteraDouble> revRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<CanteraDouble> netRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<CanteraDouble> fwdRatesOfProgress_ddCi() override;
    Eigen::SparseMatrix<CanteraDouble> revRatesOfProgress_ddCi() override;
    Eigen::SparseMatrix<CanteraDouble> netRatesOfProgress_ddCi() override;
    //! @}

    //! @name Rate calculation intermediate methods
    //! @{

    void updateROP() override;

    void getThirdBodyConcentrations(CanteraDouble* concm) override;
    const vector<CanteraDouble>& thirdBodyConcentrations() const override {
        return m_concm;
    }

    //! @}

protected:
    //! @name Internal service methods
    //!
    //! @note These methods are for internal use, and seek to avoid code duplication
    //! while evaluating terms used for rate constants, rates of progress, and
    //! their derivatives.
    //! @{

    //! Multiply rate with third-body collider concentrations
    void processThirdBodies(CanteraDouble* rop);

    //! Multiply rate with inverse equilibrium constant
    void applyEquilibriumConstants(CanteraDouble* rop);

    //! Multiply rate with scaled temperature derivatives of the inverse
    //! equilibrium constant
    /*!
     *  This (scaled) derivative is handled by a finite difference.
     */
    void applyEquilibriumConstants_ddT(CanteraDouble* drkcn);

    //! Process temperature derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddT(const vector<CanteraDouble>& in, CanteraDouble* drop);

    //! Process pressure derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddP(const vector<CanteraDouble>& in, CanteraDouble* drop);

    //! Process concentration (molar density) derivative
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    //! @param mass_action  boolean indicating whether law of mass action applies
    void process_ddC(StoichManagerN& stoich, const vector<CanteraDouble>& in,
                     CanteraDouble* drop, bool mass_action=true);

    //! Process derivatives
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param ddX true: w.r.t mole fractions false: w.r.t species concentrations
    //! @return a sparse matrix of derivative contributions for each reaction of
    //! dimensions nTotalReactions by nTotalSpecies
    Eigen::SparseMatrix<CanteraDouble> calculateCompositionDerivatives(
        StoichManagerN& stoich, const vector<CanteraDouble>& in, bool ddX=true);

    //! Helper function ensuring that all rate derivatives can be calculated
    //! @param name  method name used for error output
    //! @throw CanteraError if ideal gas assumption does not hold
    void assertDerivativesValid(const string& name);

    //! @}

    //! Difference between the global reactants order and the global products
    //! order. Of type "CanteraDouble" to account for the fact that we can have real-
    //! valued stoichiometries.
    vector<CanteraDouble> m_dn;

    ThirdBodyCalc m_multi_concm; //!< used with MultiRate evaluator

    //! Third body concentrations
    vector<CanteraDouble> m_concm;

    //! Activity concentrations, as calculated by ThermoPhase::getActivityConcentrations
    vector<CanteraDouble> m_act_conc;

    //! Physical concentrations, as calculated by ThermoPhase::getConcentrations
    vector<CanteraDouble> m_phys_conc;

    //! Derivative settings
    bool m_jac_skip_third_bodies;
    bool m_jac_skip_falloff;
    CanteraDouble m_jac_rtol_delta;

    bool m_ROP_ok = false;

    //! Buffers for partial rop results with length nReactions()
    vector<CanteraDouble> m_rbuf0;
    vector<CanteraDouble> m_rbuf1;
    vector<CanteraDouble> m_rbuf2;
    vector<CanteraDouble> m_kf0; //!< Forward rate constants without perturbation
    vector<CanteraDouble> m_sbuf0;
    vector<CanteraDouble> m_state;
    vector<CanteraDouble> m_grt; //!< Standard chemical potentials for each species
};

}

#endif
