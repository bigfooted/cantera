/**
 *  @file ThirdBodyCalc.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_THIRDBODYCALC_H
#define CT_THIRDBODYCALC_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Calculate and apply third-body effects on reaction rates, including non-
//! unity third-body efficiencies.
//! @ingroup rateEvaluators
class ThirdBodyCalc
{
public:
    //! Install reaction that uses third-body effects in ThirdBodyCalc manager
    void install(size_t rxnNumber, const map<size_t, CanteraDouble>& efficiencies,
                 CanteraDouble default_efficiency, bool mass_action) {
        m_reaction_index.push_back(rxnNumber);
        m_default.push_back(default_efficiency);

        if (mass_action) {
            m_mass_action_index.push_back(m_reaction_index.size() - 1);
        } else {
            m_no_mass_action_index.push_back(m_reaction_index.size() - 1);
        }

        m_species.emplace_back();
        m_eff.emplace_back();
        for (const auto& [k, efficiency] : efficiencies) {
            AssertTrace(k != npos);
            m_species.back().push_back(k);
            m_eff.back().push_back(efficiency - default_efficiency);
            m_efficiencyList.emplace_back(
                static_cast<int>(rxnNumber),
                static_cast<int>(k), efficiency - default_efficiency);
        }
    }

    //! Resize the sparse coefficient matrix
    void resizeCoeffs(size_t nSpc, size_t nRxn) {
        // Sparse Efficiency coefficient matrix
        Eigen::SparseMatrix<CanteraDouble> efficiencies;
        efficiencies.setZero();
        efficiencies.resize(nRxn, nSpc);
        efficiencies.reserve(m_efficiencyList.size());
        efficiencies.setFromTriplets(
            m_efficiencyList.begin(), m_efficiencyList.end());

        // derivative matrix multipliers
        vector<Eigen::Triplet<CanteraDouble>> triplets;
        triplets.reserve(m_reaction_index.size() * nSpc);
        for (size_t i = 0; i < m_default.size(); i++) {
            if (m_default[i] != 0) {
                for (size_t j = 0; j < nSpc; j++) {
                    triplets.emplace_back(
                        static_cast<int>(m_reaction_index[i]),
                        static_cast<int>(j), m_default[i]);
                }
            }
        }
        Eigen::SparseMatrix<CanteraDouble> defaults(nRxn, nSpc);
        defaults.reserve(triplets.size());
        defaults.setFromTriplets(triplets.begin(), triplets.end());
        m_multipliers = efficiencies + defaults;
    }

    //! Update third-body concentrations in full vector
    void update(const vector<CanteraDouble>& conc, CanteraDouble ctot, CanteraDouble* concm) const {
        for (size_t i = 0; i < m_reaction_index.size(); i++) {
            CanteraDouble sum = 0.0;
            for (size_t j = 0; j < m_species[i].size(); j++) {
                sum += m_eff[i][j] * conc[m_species[i][j]];
            }
            concm[m_reaction_index[i]] = m_default[i] * ctot + sum;
        }
    }

    //! Multiply output with effective third-body concentration
    void multiply(CanteraDouble* output, const CanteraDouble* concm) {
        for (size_t i = 0; i < m_mass_action_index.size(); i++) {
            size_t ix = m_reaction_index[m_mass_action_index[i]];
            output[ix] *= concm[ix];
        }
    }

    //! Calculate derivatives with respect to species concentrations.
    /*!
     *  @param product   Product of law of mass action and rate terms.
     */
    Eigen::SparseMatrix<CanteraDouble> derivatives(const CanteraDouble* product) {
        Eigen::Map<const Cantera::VectorXd> mapped(product, m_multipliers.rows());
        return mapped.asDiagonal() * m_multipliers;
    }

    //! Scale entries involving third-body collider in law of mass action by factor
    void scale(const CanteraDouble* in, CanteraDouble* out, CanteraDouble factor) const {
        for (size_t i = 0; i < m_mass_action_index.size(); i++) {
            size_t ix = m_reaction_index[m_mass_action_index[i]];
            out[ix] = factor * in[ix];
        }
    }

    //! Scale entries involving third-body collider in rate expression
    //! by third-body concentration and factor
    void scaleM(const CanteraDouble* in, CanteraDouble* out,
                const CanteraDouble* concm, CanteraDouble factor) const
    {
        for (size_t i = 0; i < m_no_mass_action_index.size(); i++) {
            size_t ix = m_reaction_index[m_no_mass_action_index[i]];
            out[ix] = factor * concm[ix] * in[ix];
        }
    }

    //! Return boolean indicating whether ThirdBodyCalc is empty
    bool empty() const {
        return m_reaction_index.empty();
    }

protected:
    //! Indices of reactions that use third-bodies within vector of concentrations
    vector<size_t> m_reaction_index;

    //! Indices within m_reaction_index of reactions that consider third-body effects
    //! in the law of mass action
    vector<size_t> m_mass_action_index;

    //! Indices within m_reaction_index of reactions that consider third-body effects
    //! in the rate expression
    vector<size_t> m_no_mass_action_index;

    //! m_species[i][j] is the index of the j-th species in reaction i.
    vector<vector<size_t>> m_species;

    //! m_eff[i][j] is the efficiency of the j-th species in reaction i.
    vector<vector<CanteraDouble>> m_eff;

    //! The default efficiency for each reaction
    vector<CanteraDouble> m_default;

    //! Sparse efficiency matrix (compensated for defaults)
    //! Each triplet corresponds to (reaction index, species index, efficiency)
    vector<Eigen::Triplet<CanteraDouble>> m_efficiencyList;

    //! Sparse derivative multiplier matrix
    Eigen::SparseMatrix<CanteraDouble> m_multipliers;
};

}

#endif
