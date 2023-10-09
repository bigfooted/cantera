//! @file Wall.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/numerics/Func1.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

bool WallBase::install(ReactorBase& rleft, ReactorBase& rright)
{
    // check if wall is already installed
    if (m_left || m_right) {
        return false;
    }
    m_left =  &rleft;
    m_right = &rright;
    m_left->addWall(*this, 0);
    m_right->addWall(*this, 1);
    return true;
}

void WallBase::setArea(CanteraDouble a) {
    m_area = a;
}

CanteraDouble Wall::velocity() const {
    if (m_vf) {
        return m_vf->eval(m_time);
    }
    return 0.;
}

CanteraDouble Wall::vdot(CanteraDouble t)
{
    warn_deprecated("Wall::vdot",
        "To be removed after Cantera 3.0; replaceable by 'expansionRate'.");
    CanteraDouble rate = m_k * m_area * (m_left->pressure() - m_right->pressure());

    if (m_vf) {
        rate += m_area * m_vf->eval(t);
    }
    return rate;
}

CanteraDouble Wall::expansionRate()
{
    CanteraDouble rate = m_k * m_area * (m_left->pressure() - m_right->pressure());

    if (m_vf) {
        rate += m_area * m_vf->eval(m_time);
    }
    return rate;
}

CanteraDouble Wall::heatFlux() const {
    if (m_qf) {
        return m_qf->eval(m_time);
    }
    return 0.;
}

CanteraDouble Wall::Q(CanteraDouble t)
{
    warn_deprecated("Wall::Q",
        "To be removed after Cantera 3.0; replaceable by 'heatRate'.");
    CanteraDouble q1 = (m_area * m_rrth) *
                (m_left->temperature() - m_right->temperature());
    if (m_emiss > 0.0) {
        CanteraDouble tl = m_left->temperature();
        CanteraDouble tr = m_right->temperature();
        q1 += m_emiss * m_area * StefanBoltz * (tl*tl*tl*tl - tr*tr*tr*tr);
    }

    if (m_qf) {
        q1 += m_area * m_qf->eval(t);
    }
    return q1;
}

CanteraDouble Wall::heatRate()
{
    CanteraDouble q1 = (m_area * m_rrth) *
                (m_left->temperature() - m_right->temperature());
    if (m_emiss > 0.0) {
        CanteraDouble tl = m_left->temperature();
        CanteraDouble tr = m_right->temperature();
        q1 += m_emiss * m_area * StefanBoltz * (tl*tl*tl*tl - tr*tr*tr*tr);
    }

    if (m_qf) {
        q1 += m_area * m_qf->eval(m_time);
    }
    return q1;
}

}
