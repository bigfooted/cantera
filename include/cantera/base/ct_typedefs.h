#pragma once
#ifndef AD_ENABLED
# define AD_ENABLED 0
#endif

#if AD_ENABLED
#include <codi.hpp>
using CanteraDoublePassive = double;

// Select AD type.
#if defined(AD_FORWARD)
using CanteraDouble = codi::RealForward;
#elif defined(AD_REVERSE)
using CanteraDouble = codi::RealReverse;
#else
#error AD mode configured with unknown type.
#endif

#else
using CanteraDoublePassive = double;
using CanteraDouble = double;
#endif
