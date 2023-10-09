/**
 * @file PDSS_SSVol.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/VPStandardStateTP.h"

namespace Cantera
{

PDSS_SSVol::PDSS_SSVol()
    : TCoeff_(4, 0.0)
{
}

void PDSS_SSVol::setVolumePolynomial(CanteraDouble* coeffs) {
    for (size_t i = 0; i < 4; i++) {
        TCoeff_[i] = coeffs[i];
    }
    volumeModel_ = SSVolume_Model::tpoly;
}

void PDSS_SSVol::setDensityPolynomial(CanteraDouble* coeffs) {
    for (size_t i = 0; i < 4; i++) {
        TCoeff_[i] = coeffs[i];
    }
    volumeModel_ = SSVolume_Model::density_tpoly;
}

void PDSS_SSVol::getParameters(AnyMap& eosNode) const
{
    PDSS::getParameters(eosNode);
    vector<AnyValue> data(4);
    if (volumeModel_ == SSVolume_Model::density_tpoly) {
        eosNode["model"] = "density-temperature-polynomial";
        data[0].setQuantity(TCoeff_[0], "kg/m^3");
        data[1].setQuantity(TCoeff_[1], "kg/m^3/K");
        data[2].setQuantity(TCoeff_[2], "kg/m^3/K^2");
        data[3].setQuantity(TCoeff_[3], "kg/m^3/K^3");
    } else {
        eosNode["model"] = "molar-volume-temperature-polynomial";
        data[0].setQuantity(TCoeff_[0], "m^3/kmol");
        data[1].setQuantity(TCoeff_[1], "m^3/kmol/K");
        data[2].setQuantity(TCoeff_[2], "m^3/kmol/K^2");
        data[3].setQuantity(TCoeff_[3], "m^3/kmol/K^3");
    }
    eosNode["data"] = std::move(data);
}

void PDSS_SSVol::initThermo()
{
    PDSS::initThermo();
    if (m_input.hasKey("model")) {
        const string& model = m_input["model"].asString();
        auto& data = m_input["data"].asVector<AnyValue>(4);
        if (model == "density-temperature-polynomial") {
            CanteraDouble coeffs[] {
                m_input.units().convert(data[0], "kg/m^3"),
                m_input.units().convert(data[1], "kg/m^3/K"),
                m_input.units().convert(data[2], "kg/m^3/K^2"),
                m_input.units().convert(data[3], "kg/m^3/K^3"),
            };
            setDensityPolynomial(coeffs);
        } else if (model == "molar-volume-temperature-polynomial") {
            CanteraDouble coeffs[] {
                m_input.units().convert(data[0], "m^3/kmol"),
                m_input.units().convert(data[1], "m^3/kmol/K"),
                m_input.units().convert(data[2], "m^3/kmol/K^2"),
                m_input.units().convert(data[3], "m^3/kmol/K^3"),
            };
            setVolumePolynomial(coeffs);
        }
    }
    m_minTemp = m_spthermo->minTemp();
    m_maxTemp = m_spthermo->maxTemp();
    m_p0 = m_spthermo->refPressure();
}

CanteraDouble PDSS_SSVol::intEnergy_mole() const
{
    CanteraDouble pV = m_pres * m_Vss;
    return m_h0_RT * GasConstant * m_temp - pV;
}

CanteraDouble PDSS_SSVol::cv_mole() const
{
    return (cp_mole() - m_V0);
}

void PDSS_SSVol::calcMolarVolume()
{
    if (volumeModel_ == SSVolume_Model::tpoly) {
        m_Vss = TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));
        m_V0 = m_Vss;
        dVdT_ = TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];
        d2VdT2_ = 2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];
    } else if (volumeModel_ == SSVolume_Model::density_tpoly) {
        CanteraDouble dens = TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));
        m_Vss = m_mw / dens;
        m_V0 = m_Vss;
        CanteraDouble dens2 = dens * dens;
        CanteraDouble ddensdT = TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];
        CanteraDouble d2densdT2 = 2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];
        dVdT_ = - m_mw / dens2 * ddensdT;
        d2VdT2_ = 2.0 * m_mw / (dens2 * dens) * ddensdT * ddensdT - m_mw / dens2 * d2densdT2;
    } else {
        throw NotImplementedError("PDSS_SSVol::calcMolarVolume");
    }
}

void PDSS_SSVol::setPressure(CanteraDouble p)
{
    m_pres = p;
    CanteraDouble deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        m_hss_RT = m_h0_RT;
        m_sss_R = m_s0_R;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R;
    } else {
        CanteraDouble del_pRT = deltaP / (GasConstant * m_temp);
        CanteraDouble sV_term = - deltaP / GasConstant * dVdT_;
        m_hss_RT = m_h0_RT + sV_term + del_pRT * m_Vss;
        m_sss_R = m_s0_R + sV_term;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R - m_temp * deltaP * d2VdT2_;
    }
}

void PDSS_SSVol::setTemperature(CanteraDouble temp)
{
    m_temp = temp;
    m_spthermo->updatePropertiesTemp(temp, &m_cp0_R, &m_h0_RT, &m_s0_R);
    calcMolarVolume();
    m_g0_RT = m_h0_RT - m_s0_R;
    CanteraDouble deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        m_hss_RT = m_h0_RT;
        m_sss_R = m_s0_R;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R;
    } else {
        CanteraDouble del_pRT = deltaP / (GasConstant * m_temp);
        CanteraDouble sV_term = - deltaP / GasConstant * dVdT_;
        m_hss_RT = m_h0_RT + sV_term + del_pRT * m_Vss;
        m_sss_R = m_s0_R + sV_term;
        m_gss_RT = m_hss_RT - m_sss_R;
        m_cpss_R = m_cp0_R - m_temp * deltaP * d2VdT2_;
    }
}

void PDSS_SSVol::setState_TP(CanteraDouble temp, CanteraDouble pres)
{
    m_pres = pres;
    setTemperature(temp);
}

CanteraDouble PDSS_SSVol::satPressure(CanteraDouble t)
{
    return 1.0E-200;
}

}
