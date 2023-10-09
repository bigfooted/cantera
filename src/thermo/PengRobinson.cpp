//! @file PengRobinson.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>

namespace bmt = boost::math::tools;

namespace Cantera
{

const CanteraDouble PengRobinson::omega_a = 4.5723552892138218E-01;
const CanteraDouble PengRobinson::omega_b = 7.77960739038885E-02;
const CanteraDouble PengRobinson::omega_vc = 3.07401308698703833E-01;

PengRobinson::PengRobinson(const string& infile, const string& id_)
{
    initThermoFile(infile, id_);
}

void PengRobinson::setSpeciesCoeffs(const string& species, CanteraDouble a, CanteraDouble b, CanteraDouble w)
{
    size_t k = speciesIndex(species);
    if (k == npos) {
        throw CanteraError("PengRobinson::setSpeciesCoeffs",
            "Unknown species '{}'.", species);
    }

    // Calculate value of kappa (independent of temperature)
    // w is an acentric factor of species
    if (w <= 0.491) {
        m_kappa[k] = 0.37464 + 1.54226*w - 0.26992*w*w;
    } else {
        m_kappa[k] = 0.374642 + 1.487503*w - 0.164423*w*w + 0.016666*w*w*w;
    }
    m_acentric[k] = w; // store the original acentric factor to enable serialization

    // Calculate alpha (temperature dependent interaction parameter)
    CanteraDouble critTemp = speciesCritTemperature(a, b);
    CanteraDouble sqt_T_r = sqrt(temperature() / critTemp);
    CanteraDouble sqt_alpha = 1 + m_kappa[k] * (1 - sqt_T_r);
    m_alpha[k] = sqt_alpha*sqt_alpha;

    m_a_coeffs(k,k) = a;
    CanteraDouble aAlpha_k = a*m_alpha[k];
    m_aAlpha_binary(k,k) = aAlpha_k;

    // standard mixing rule for cross-species interaction term
    for (size_t j = 0; j < m_kk; j++) {
        if (k == j) {
            continue;
        }
        CanteraDouble a0kj = sqrt(m_a_coeffs(j,j) * a);
        CanteraDouble aAlpha_j = a*m_alpha[j];
        CanteraDouble a_Alpha = sqrt(aAlpha_j*aAlpha_k);
        if (m_a_coeffs(j, k) == 0) {
            m_a_coeffs(j, k) = a0kj;
            m_aAlpha_binary(j, k) = a_Alpha;
            m_a_coeffs(k, j) = a0kj;
            m_aAlpha_binary(k, j) = a_Alpha;
        }
    }
    m_b_coeffs[k] = b;
}

void PengRobinson::setBinaryCoeffs(const string& species_i,
        const string& species_j, CanteraDouble a0)
{
    size_t ki = speciesIndex(species_i);
    if (ki == npos) {
        throw CanteraError("PengRobinson::setBinaryCoeffs",
            "Unknown species '{}'.", species_i);
    }
    size_t kj = speciesIndex(species_j);
    if (kj == npos) {
        throw CanteraError("PengRobinson::setBinaryCoeffs",
            "Unknown species '{}'.", species_j);
    }

    m_a_coeffs(ki, kj) = m_a_coeffs(kj, ki) = a0;
    m_binaryParameters[species_i][species_j] = a0;
    m_binaryParameters[species_j][species_i] = a0;
    // Calculate alpha_ij
    CanteraDouble alpha_ij = m_alpha[ki] * m_alpha[kj];
    m_aAlpha_binary(ki, kj) = m_aAlpha_binary(kj, ki) = a0*alpha_ij;
}

// ------------Molar Thermodynamic Properties -------------------------

CanteraDouble PengRobinson::cp_mole() const
{
    _updateReferenceStateThermo();
    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble vpb = mv + (1 + Sqrt2) * m_b;
    CanteraDouble vmb = mv + (1 - Sqrt2) * m_b;
    calculatePressureDerivatives();
    CanteraDouble cpref = GasConstant * mean_X(m_cp0_R);
    CanteraDouble dHdT_V = cpref + mv * m_dpdT - GasConstant
                    + 1.0 / (2.0 * Sqrt2 * m_b) * log(vpb / vmb) * T * d2aAlpha_dT2();
    return dHdT_V - (mv + T * m_dpdT / m_dpdV) * m_dpdT;
}

CanteraDouble PengRobinson::cv_mole() const
{
    _updateReferenceStateThermo();
    CanteraDouble T = temperature();
    calculatePressureDerivatives();
    return (cp_mole() + T * m_dpdT * m_dpdT / m_dpdV);
}

CanteraDouble PengRobinson::pressure() const
{
    _updateReferenceStateThermo();
    // Get a copy of the private variables stored in the State object
    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble denom = mv * mv + 2 * mv * m_b - m_b * m_b;
    return GasConstant * T / (mv - m_b) - m_aAlpha_mix / denom;
}

CanteraDouble PengRobinson::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void PengRobinson::getActivityCoefficients(CanteraDouble* ac) const
{
    CanteraDouble mv = molarVolume();
    CanteraDouble vpb2 = mv + (1 + Sqrt2) * m_b;
    CanteraDouble vmb2 = mv + (1 - Sqrt2) * m_b;
    CanteraDouble vmb = mv - m_b;
    CanteraDouble pres = pressure();

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
        }
    }
    CanteraDouble num = 0;
    CanteraDouble denom = 2 * Sqrt2 * m_b * m_b;
    CanteraDouble denom2 = m_b * (mv * mv + 2 * mv * m_b - m_b * m_b);
    CanteraDouble RT_ = RT();
    for (size_t k = 0; k < m_kk; k++) {
        num = 2 * m_b * m_pp[k] - m_aAlpha_mix * m_b_coeffs[k];
        ac[k] = (-RT_ * log(pres * mv/ RT_) + RT_ * log(mv / vmb)
                 + RT_ * m_b_coeffs[k] / vmb
                 - (num /denom) * log(vpb2/vmb2)
                 - m_aAlpha_mix * m_b_coeffs[k] * mv/denom2
                );
    }
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(ac[k]/ RT_);
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void PengRobinson::getChemPotentials(CanteraDouble* mu) const
{
    getGibbs_ref(mu);
    CanteraDouble RT_ = RT();
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT_ * (log(xx));
    }

    CanteraDouble mv = molarVolume();
    CanteraDouble vmb = mv - m_b;
    CanteraDouble vpb2 = mv + (1 + Sqrt2) * m_b;
    CanteraDouble vmb2 = mv + (1 - Sqrt2) * m_b;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
        }
    }
    CanteraDouble pres = pressure();
    CanteraDouble refP = refPressure();
    CanteraDouble denom = 2 * Sqrt2 * m_b * m_b;
    CanteraDouble denom2 = m_b * (mv * mv + 2 * mv * m_b - m_b * m_b);

    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble num = 2 * m_b * m_pp[k] - m_aAlpha_mix * m_b_coeffs[k];

        mu[k] += RT_ * log(pres/refP) - RT_ * log(pres * mv / RT_)
                 + RT_ * log(mv / vmb) + RT_ * m_b_coeffs[k] / vmb
                 - (num /denom) * log(vpb2/vmb2)
                 - m_aAlpha_mix * m_b_coeffs[k] * mv/denom2;
    }
}

void PengRobinson::getPartialMolarEnthalpies(CanteraDouble* hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar, hbar+m_kk, hbar, RT());
    vector<CanteraDouble> tmp;
    tmp.resize(m_kk,0.0);

    // We calculate m_dpdni
    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble vmb = mv - m_b;
    CanteraDouble vpb2 = mv + (1 + Sqrt2) * m_b;
    CanteraDouble vmb2 = mv + (1 - Sqrt2) * m_b;
    CanteraDouble daAlphadT = daAlpha_dT();

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        tmp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            CanteraDouble grad_aAlpha = m_dalphadT[i]/m_alpha[i] + m_dalphadT[k]/m_alpha[k];
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
            tmp[k] +=moleFractions_[i] * m_aAlpha_binary(k, i) * grad_aAlpha;
        }
    }

    CanteraDouble denom = mv * mv + 2 * mv * m_b - m_b * m_b;
    CanteraDouble denom2 = denom * denom;
    CanteraDouble RT_ = RT();
    for (size_t k = 0; k < m_kk; k++) {
        m_dpdni[k] = RT_ / vmb + RT_ * m_b_coeffs[k] / (vmb * vmb)
            - 2.0 * m_pp[k] / denom + 2 * vmb * m_aAlpha_mix * m_b_coeffs[k] / denom2;
    }

    CanteraDouble fac = T * daAlphadT - m_aAlpha_mix;
    calculatePressureDerivatives();
    CanteraDouble fac2 = mv + T * m_dpdT / m_dpdV;
    CanteraDouble fac3 = 2 * Sqrt2 * m_b * m_b;
    CanteraDouble fac4 = 0;
    for (size_t k = 0; k < m_kk; k++) {
        fac4 = T*tmp[k] -2 * m_pp[k];
        CanteraDouble hE_v = mv * m_dpdni[k] - RT_
                     - m_b_coeffs[k] / fac3 * log(vpb2 / vmb2) * fac
                     + (mv * m_b_coeffs[k]) / (m_b * denom) * fac
                     + 1 / (2 * Sqrt2 * m_b) * log(vpb2 / vmb2) * fac4;
        hbar[k] = hbar[k] + hE_v;
        hbar[k] -= fac2 * m_dpdni[k];
    }
}

void PengRobinson::getPartialMolarEntropies(CanteraDouble* sbar) const
{
    // Using the identity : (hk - T*sk) = gk
    CanteraDouble T = temperature();
    getPartialMolarEnthalpies(sbar);
    getChemPotentials(m_tmpV.data());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] = (sbar[k] - m_tmpV[k])/T;
    }
}

void PengRobinson::getPartialMolarIntEnergies(CanteraDouble* ubar) const
{
    // u_i = h_i - p*v_i
    CanteraDouble p = pressure();
    getPartialMolarEnthalpies(ubar);
    getPartialMolarVolumes(m_tmpV.data());
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = ubar[k] - p*m_tmpV[k];
    }
}

void PengRobinson::getPartialMolarCp(CanteraDouble* cpbar) const
{
    throw NotImplementedError("PengRobinson::getPartialMolarCp");
}

void PengRobinson::getPartialMolarVolumes(CanteraDouble* vbar) const
{
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            m_pp[k] += moleFractions_[i] * m_aAlpha_binary(k, i);
        }
    }

    CanteraDouble mv = molarVolume();
    CanteraDouble vmb = mv - m_b;
    CanteraDouble vpb = mv + m_b;
    CanteraDouble fac = mv * mv + 2 * mv * m_b - m_b * m_b;
    CanteraDouble fac2 = fac * fac;
    CanteraDouble RT_ = RT();

    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble num = RT_ + RT_ * m_b/ vmb + RT_ * m_b_coeffs[k] / vmb
                     + RT_ * m_b * m_b_coeffs[k] /(vmb * vmb) - 2 * mv * m_pp[k] / fac
                     + 2 * mv * vmb * m_aAlpha_mix * m_b_coeffs[k] / fac2;
        CanteraDouble denom = pressure() + RT_ * m_b / (vmb * vmb) + m_aAlpha_mix/fac
                       - 2 * mv* vpb * m_aAlpha_mix / fac2;
        vbar[k] = num / denom;
    }
}

CanteraDouble PengRobinson::speciesCritTemperature(CanteraDouble a, CanteraDouble b) const
{
    if (b <= 0.0) {
        return 1000000.;
    } else if (a <= 0.0) {
        return 0.0;
    } else {
        return a * omega_b / (b * omega_a * GasConstant);
    }
}

bool PengRobinson::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        m_a_coeffs.resize(m_kk, m_kk, 0.0);
        m_b_coeffs.push_back(0.0);
        m_aAlpha_binary.resize(m_kk, m_kk, 0.0);
        m_kappa.push_back(0.0);
        m_acentric.push_back(0.0);
        m_alpha.push_back(0.0);
        m_dalphadT.push_back(0.0);
        m_d2alphadT2.push_back(0.0);
        m_pp.push_back(0.0);
        m_partialMolarVolumes.push_back(0.0);
        m_dpdni.push_back(0.0);
        m_coeffSource.push_back(CoeffSource::EoS);
    }
    return added;
}

void PengRobinson::initThermo()
{
    // Contents of 'critical-properties.yaml', loaded later if needed
    AnyMap critPropsDb;
    std::unordered_map<string, AnyMap*> dbSpecies;

    for (auto& [name, species] : m_species) {
        auto& data = species->input;
        size_t k = speciesIndex(name);
        if (m_a_coeffs(k, k) != 0.0) {
            continue;
        }
        bool foundCoeffs = false;
        if (data.hasKey("equation-of-state") &&
            data["equation-of-state"].hasMapWhere("model", "Peng-Robinson"))
        {
            // Read a and b coefficients and acentric factor w_ac from species input
            // information, specified in a YAML input file.
            auto eos = data["equation-of-state"].getMapWhere(
                "model", "Peng-Robinson");
            if (eos.hasKey("a") && eos.hasKey("b") && eos.hasKey("acentric-factor")) {
                CanteraDouble a0 = eos.convert("a", "Pa*m^6/kmol^2");
                CanteraDouble b = eos.convert("b", "m^3/kmol");
                // unitless acentric factor:
                CanteraDouble w = eos["acentric-factor"].asDouble();
                setSpeciesCoeffs(name, a0, b, w);
                foundCoeffs = true;
            }

            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& [name2, coeff] : binary_a) {
                    CanteraDouble a0 = units.convert(coeff, "Pa*m^6/kmol^2");
                    setBinaryCoeffs(name, name2, a0);
                }
            }
            if (foundCoeffs) {
                m_coeffSource[k] = CoeffSource::EoS;
                continue;
            }
        }

        // Coefficients have not been populated from model-specific input
        CanteraDouble Tc = NAN, Pc = NAN, omega_ac = NAN;
        if (data.hasKey("critical-parameters")) {
            // Use critical state information stored in the species entry to
            // calculate a, b, and the acentric factor.
            auto& critProps = data["critical-parameters"].as<AnyMap>();
            Tc = critProps.convert("critical-temperature", "K");
            Pc = critProps.convert("critical-pressure", "Pa");
            omega_ac = critProps["acentric-factor"].asDouble();
            m_coeffSource[k] = CoeffSource::CritProps;
        } else {
            // Search 'crit-properties.yaml' to find Tc and Pc. Load data if needed.
            if (critPropsDb.empty()) {
                critPropsDb = AnyMap::fromYamlFile("critical-properties.yaml");
                dbSpecies = critPropsDb["species"].asMap("name");
            }

            // All names in critical-properties.yaml are upper case
            auto ucName = boost::algorithm::to_upper_copy(name);
            if (dbSpecies.count(ucName)) {
                auto& spec = *dbSpecies.at(ucName);
                auto& critProps = spec["critical-parameters"].as<AnyMap>();
                Tc = critProps.convert("critical-temperature", "K");
                Pc = critProps.convert("critical-pressure", "Pa");
                omega_ac = critProps["acentric-factor"].asDouble();
                m_coeffSource[k] = CoeffSource::Database;
            }
        }

        // Check if critical properties were found in either location
        if (!isnan(Tc)) {
            CanteraDouble a = omega_a * std::pow(GasConstant * Tc, 2) / Pc;
            CanteraDouble b = omega_b * GasConstant * Tc / Pc;
            setSpeciesCoeffs(name, a, b, omega_ac);
        } else {
            throw InputFileError("PengRobinson::initThermo", data,
            "No Peng-Robinson model parameters or critical properties found for "
            "species '{}'", name);
        }
    }
}

void PengRobinson::getSpeciesParameters(const string& name, AnyMap& speciesNode) const
{
    MixtureFugacityTP::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name);
    checkSpeciesIndex(k);

    // Pure species parameters
    if (m_coeffSource[k] == CoeffSource::EoS) {
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Peng-Robinson", true);
        eosNode["a"].setQuantity(m_a_coeffs(k, k), "Pa*m^6/kmol^2");
        eosNode["b"].setQuantity(m_b_coeffs[k], "m^3/kmol");
        eosNode["acentric-factor"] = m_acentric[k];
    } else if (m_coeffSource[k] == CoeffSource::CritProps) {
        auto& critProps = speciesNode["critical-parameters"];
        CanteraDouble Tc = speciesCritTemperature(m_a_coeffs(k, k), m_b_coeffs[k]);
        CanteraDouble Pc = omega_b * GasConstant * Tc / m_b_coeffs[k];
        critProps["critical-temperature"].setQuantity(Tc, "K");
        critProps["critical-pressure"].setQuantity(Pc, "Pa");
        critProps["acentric-factor"] = m_acentric[k];
    }
    // Nothing to do in the case where the parameters are from the database

    if (m_binaryParameters.count(name)) {
        // Include binary parameters regardless of where the pure species parameters
        // were found
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Peng-Robinson", true);
        AnyMap bin_a;
        for (const auto& [species, coeff] : m_binaryParameters.at(name)) {
            bin_a[species].setQuantity(coeff, "Pa*m^6/kmol^2");
        }
        eosNode["binary-a"] = std::move(bin_a);
    }
}

CanteraDouble PengRobinson::sresid() const
{
    CanteraDouble molarV = molarVolume();
    CanteraDouble hh = m_b / molarV;
    CanteraDouble zz = z();
    CanteraDouble alpha_1 = daAlpha_dT();
    CanteraDouble vpb = molarV + (1.0 + Sqrt2) * m_b;
    CanteraDouble vmb = molarV + (1.0 - Sqrt2) * m_b;
    CanteraDouble fac = alpha_1 / (2.0 * Sqrt2 * m_b);
    CanteraDouble sresid_mol_R = log(zz*(1.0 - hh)) + fac * log(vpb / vmb) / GasConstant;
    return GasConstant * sresid_mol_R;
}

CanteraDouble PengRobinson::hresid() const
{
    CanteraDouble molarV = molarVolume();
    CanteraDouble zz = z();
    CanteraDouble aAlpha_1 = daAlpha_dT();
    CanteraDouble T = temperature();
    CanteraDouble vpb = molarV + (1 + Sqrt2) * m_b;
    CanteraDouble vmb = molarV + (1 - Sqrt2) * m_b;
    CanteraDouble fac = 1 / (2.0 * Sqrt2 * m_b);
    return GasConstant * T * (zz - 1.0)
        + fac * log(vpb / vmb) * (T * aAlpha_1 - m_aAlpha_mix);
}

CanteraDouble PengRobinson::liquidVolEst(CanteraDouble T, CanteraDouble& presGuess) const
{
    CanteraDouble v = m_b * 1.1;
    CanteraDouble atmp, btmp, aAlphatmp;
    calculateAB(atmp, btmp, aAlphatmp);
    CanteraDouble pres = std::max(psatEst(T), presGuess);
    CanteraDouble Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
        int nsol = solveCubic(T, pres, atmp, btmp, aAlphatmp, Vroot);
        if (nsol == 1 || nsol == 2) {
            CanteraDouble pc = critPressure();
            if (pres > pc) {
                foundLiq = true;
            }
            pres *= 1.04;
        } else {
            foundLiq = true;
        }
    }

    if (foundLiq) {
        v = Vroot[0];
        presGuess = pres;
    } else {
        v = -1.0;
    }
    return v;
}

CanteraDouble PengRobinson::densityCalc(CanteraDouble T, CanteraDouble presPa, int phaseRequested,
                                 CanteraDouble rhoGuess)
{
    // It's necessary to set the temperature so that m_aAlpha_mix is set correctly.
    setTemperature(T);
    CanteraDouble tcrit = critTemperature();
    CanteraDouble mmw = meanMolecularWeight();
    if (rhoGuess == -1.0) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            CanteraDouble lqvol = liquidVolEst(T, presPa);
            rhoGuess = mmw / lqvol;
        }
    } else {
        // Assume the Gas phase initial guess, if nothing is specified to the routine
        rhoGuess = presPa * mmw / (GasConstant * T);
    }

    CanteraDouble volGuess = mmw / rhoGuess;
    m_NSolns = solveCubic(T, presPa, m_a, m_b, m_aAlpha_mix, m_Vroot);

    CanteraDouble molarVolLast = m_Vroot[0];
    if (m_NSolns >= 2) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            molarVolLast = m_Vroot[0];
        } else if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = m_Vroot[2];
        } else {
            if (volGuess > m_Vroot[1]) {
                molarVolLast = m_Vroot[2];
            } else {
                molarVolLast = m_Vroot[0];
            }
        }
    } else if (m_NSolns == 1) {
        if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT
            || phaseRequested == FLUID_UNDEFINED)
        {
            molarVolLast = m_Vroot[0];
        } else {
            return -2.0;
        }
    } else if (m_NSolns == -1) {
        if (phaseRequested >= FLUID_LIQUID_0 || phaseRequested == FLUID_UNDEFINED
            || phaseRequested == FLUID_SUPERCRIT)
        {
            molarVolLast = m_Vroot[0];
        } else if (T > tcrit) {
            molarVolLast = m_Vroot[0];
        } else {
            return -2.0;
        }
    } else {
        molarVolLast = m_Vroot[0];
        return -1.0;
    }
    return mmw / molarVolLast;
}

CanteraDouble PengRobinson::densSpinodalLiquid() const
{
    CanteraDouble Vroot[3];
    CanteraDouble T = temperature();
    int nsol = solveCubic(T, pressure(), m_a, m_b, m_aAlpha_mix, Vroot);
    if (nsol != 3) {
        return critDensity();
    }

    auto resid = [this, T](CanteraDouble v) {
        CanteraDouble pp;
        return dpdVCalc(T, v, pp);
    };

    boost::uintmax_t maxiter = 100;
    auto [lower, upper] = bmt::toms748_solve(
        resid, Vroot[0], Vroot[1], bmt::eps_tolerance<CanteraDouble>(48), maxiter);

    CanteraDouble mmw = meanMolecularWeight();
    return mmw / (0.5 * (lower + upper));
}

CanteraDouble PengRobinson::densSpinodalGas() const
{
    CanteraDouble Vroot[3];
    CanteraDouble T = temperature();
    int nsol = solveCubic(T, pressure(), m_a, m_b, m_aAlpha_mix, Vroot);
    if (nsol != 3) {
        return critDensity();
    }

    auto resid = [this, T](CanteraDouble v) {
        CanteraDouble pp;
        return dpdVCalc(T, v, pp);
    };

    boost::uintmax_t maxiter = 100;
    auto [lower, upper] = bmt::toms748_solve(
        resid, Vroot[1], Vroot[2], bmt::eps_tolerance<CanteraDouble>(48), maxiter);

    CanteraDouble mmw = meanMolecularWeight();
    return mmw / (0.5 * (lower + upper));
}

CanteraDouble PengRobinson::dpdVCalc(CanteraDouble T, CanteraDouble molarVol, CanteraDouble& presCalc) const
{
    CanteraDouble denom = molarVol * molarVol + 2 * molarVol * m_b - m_b * m_b;
    CanteraDouble vpb = molarVol + m_b;
    CanteraDouble vmb = molarVol - m_b;
    return -GasConstant * T / (vmb * vmb) + 2 * m_aAlpha_mix * vpb / (denom*denom);
}

CanteraDouble PengRobinson::isothermalCompressibility() const
{
    calculatePressureDerivatives();
    return -1 / (molarVolume() * m_dpdV);
}

CanteraDouble PengRobinson::thermalExpansionCoeff() const
{
    calculatePressureDerivatives();
    return -m_dpdT / (molarVolume() * m_dpdV);
}

CanteraDouble PengRobinson::soundSpeed() const
{
    calculatePressureDerivatives();
    return molarVolume() * sqrt(-cp_mole() / cv_mole() * m_dpdV / meanMolecularWeight());
}

void PengRobinson::calculatePressureDerivatives() const
{
    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble pres;

    m_dpdV = dpdVCalc(T, mv, pres);
    CanteraDouble vmb = mv - m_b;
    CanteraDouble denom = mv * mv + 2 * mv * m_b - m_b * m_b;
    m_dpdT = (GasConstant / vmb - daAlpha_dT() / denom);
}

void PengRobinson::updateMixingExpressions()
{
    CanteraDouble temp = temperature();

    // Update individual alpha
    for (size_t j = 0; j < m_kk; j++) {
        CanteraDouble critTemp_j = speciesCritTemperature(m_a_coeffs(j,j), m_b_coeffs[j]);
        CanteraDouble sqt_alpha = 1 + m_kappa[j] * (1 - sqrt(temp / critTemp_j));
        m_alpha[j] = sqt_alpha*sqt_alpha;
    }

    // Update aAlpha_i, j
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            m_aAlpha_binary(i, j) = sqrt(m_alpha[i] * m_alpha[j]) * m_a_coeffs(i,j);
        }
    }
    calculateAB(m_a,m_b,m_aAlpha_mix);
}

void PengRobinson::calculateAB(CanteraDouble& aCalc, CanteraDouble& bCalc, CanteraDouble& aAlphaCalc) const
{
    bCalc = 0.0;
    aCalc = 0.0;
    aAlphaCalc = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        bCalc += moleFractions_[i] * m_b_coeffs[i];
        for (size_t j = 0; j < m_kk; j++) {
            aCalc += m_a_coeffs(i, j) * moleFractions_[i] * moleFractions_[j];
            aAlphaCalc += m_aAlpha_binary(i, j) * moleFractions_[i] * moleFractions_[j];
        }
    }
}

CanteraDouble PengRobinson::daAlpha_dT() const
{
    CanteraDouble daAlphadT = 0.0, k, Tc, sqtTr, coeff1, coeff2;
    for (size_t i = 0; i < m_kk; i++) {
        // Calculate first derivative of alpha for individual species
        Tc = speciesCritTemperature(m_a_coeffs(i,i), m_b_coeffs[i]);
        sqtTr = sqrt(temperature() / Tc);
        coeff1 = 1 / (Tc*sqtTr);
        coeff2 = sqtTr - 1;
        k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
    }
    // Calculate mixture derivative
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            daAlphadT += moleFractions_[i] * moleFractions_[j] * 0.5
                         * m_aAlpha_binary(i, j)
                         * (m_dalphadT[i] / m_alpha[i] + m_dalphadT[j] / m_alpha[j]);
        }
    }
    return daAlphadT;
}

CanteraDouble PengRobinson::d2aAlpha_dT2() const
{
    for (size_t i = 0; i < m_kk; i++) {
        CanteraDouble Tcrit_i = speciesCritTemperature(m_a_coeffs(i, i), m_b_coeffs[i]);
        CanteraDouble sqt_Tr = sqrt(temperature() / Tcrit_i);
        CanteraDouble coeff1 = 1 / (Tcrit_i*sqt_Tr);
        CanteraDouble coeff2 = sqt_Tr - 1;
        // Calculate first and second derivatives of alpha for individual species
        CanteraDouble k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
        m_d2alphadT2[i] = (k*k + k) * coeff1 / (2*sqt_Tr*sqt_Tr*Tcrit_i);
    }

    // Calculate mixture derivative
    CanteraDouble d2aAlphadT2 = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        CanteraDouble alphai = m_alpha[i];
        for (size_t j = 0; j < m_kk; j++) {
            CanteraDouble alphaj = m_alpha[j];
            CanteraDouble alphaij = alphai * alphaj;
            CanteraDouble term1 = m_d2alphadT2[i] / alphai + m_d2alphadT2[j] / alphaj;
            CanteraDouble term2 = 2 * m_dalphadT[i] * m_dalphadT[j] / alphaij;
            CanteraDouble term3 = m_dalphadT[i] / alphai + m_dalphadT[j] / alphaj;
            d2aAlphadT2 += 0.5 * moleFractions_[i] * moleFractions_[j]
                           * m_aAlpha_binary(i, j)
                           * (term1 + term2 - 0.5 * term3 * term3);
        }
    }
    return d2aAlphadT2;
}

void PengRobinson::calcCriticalConditions(CanteraDouble& pc, CanteraDouble& tc, CanteraDouble& vc) const
{
    if (m_b <= 0.0) {
        tc = 1000000.;
        pc = 1.0E13;
        vc = omega_vc * GasConstant * tc / pc;
        return;
    }
    if (m_a <= 0.0) {
        tc = 0.0;
        pc = 0.0;
        vc = 2.0 * m_b;
        return;
    }
    tc = m_a * omega_b / (m_b * omega_a * GasConstant);
    pc = omega_b * GasConstant * tc / m_b;
    vc = omega_vc * GasConstant * tc / pc;
}

int PengRobinson::solveCubic(CanteraDouble T, CanteraDouble pres, CanteraDouble a, CanteraDouble b, CanteraDouble aAlpha,
                             CanteraDouble Vroot[3]) const
{
    // Derive the coefficients of the cubic polynomial (in terms of molar volume v)
    CanteraDouble bsqr = b * b;
    CanteraDouble RT_p = GasConstant * T / pres;
    CanteraDouble aAlpha_p = aAlpha / pres;
    CanteraDouble an = 1.0;
    CanteraDouble bn = (b - RT_p);
    CanteraDouble cn = -(2 * RT_p * b - aAlpha_p + 3 * bsqr);
    CanteraDouble dn = (bsqr * RT_p + bsqr * b - aAlpha_p * b);

    CanteraDouble tc = a * omega_b / (b * omega_a * GasConstant);
    CanteraDouble pc = omega_b * GasConstant * tc / b;
    CanteraDouble vc = omega_vc * GasConstant * tc / pc;

    return MixtureFugacityTP::solveCubic(T, pres, a, b, aAlpha, Vroot,
                                         an, bn, cn, dn, tc, vc);
}

}
