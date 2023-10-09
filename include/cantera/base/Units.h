/**
 * @file Units.h
 * Header for unit conversion utilities, which are used to translate
 * user input from input files (See @ref inputGroup and
 * class @link Cantera::Units Units@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITS_H
#define CT_UNITS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class AnyValue;
class AnyMap;

//! @defgroup unitsGroup Unit Conversion
//! Unit conversion utilities.
//! @ingroup ioGroup

//! A representation of the units associated with a dimensional quantity.
/*!
 * Used for converting quantities between unit systems and checking for
 * dimensional consistency. Units objects are mainly used within UnitSystem
 * class to convert values from a user-specified Unit system to %Cantera's
 * base units (SI + kmol).
 * @ingroup unitsGroup
 */
class Units
{
public:
    //! Create a Units object with the specified dimensions.
    explicit Units(CanteraDouble factor=1.0, CanteraDouble mass=0, CanteraDouble length=0,
                   CanteraDouble time=0, CanteraDouble temperature=0, CanteraDouble current=0,
                   CanteraDouble quantity=0);

    //! Create an object with the specified dimensions
    //! @param units        A string representation of the units. See UnitSystem
    //!                     for a description of the formatting options.
    //! @param force_unity  ensure that conversion factor is equal to one
    explicit Units(const string& units, bool force_unity=false);

    //! Returns `true` if the specified Units are dimensionally consistent
    bool convertible(const Units& other) const;

    //! Return the factor for converting from this unit to Cantera's base
    //! units.
    CanteraDouble factor() const { return m_factor; }

    //! Multiply two Units objects, combining their conversion factors and
    //! dimensions
    Units& operator*=(const Units& other);

    //! Provide a string representation of these Units
    //! @param skip_unity  do not print '1' if conversion factor is equal to one
    string str(bool skip_unity=true) const;

    //! Raise these Units to a power, changing both the conversion factor and
    //! the dimensions of these Units.
    Units pow(CanteraDouble exponent) const;

    bool operator==(const Units& other) const;

    //! Return dimension of primary unit component
    //! ("mass", "length", "time", "temperature", "current", or "quantity")
    CanteraDouble dimension(const string& primary) const;

private:
    //! Scale the unit by the factor `k`
    void scale(CanteraDouble k) { m_factor *= k; }

    CanteraDouble m_factor = 1.0; //!< conversion factor to %Cantera base units
    CanteraDouble m_mass_dim = 0.0;
    CanteraDouble m_length_dim = 0.0;
    CanteraDouble m_time_dim = 0.0;
    CanteraDouble m_temperature_dim = 0.0;
    CanteraDouble m_current_dim = 0.0;
    CanteraDouble m_quantity_dim = 0.0;
    CanteraDouble m_pressure_dim = 0.0; //!< pseudo-dimension to track explicit pressure units
    CanteraDouble m_energy_dim = 0.0; //!< pseudo-dimension to track explicit energy units

    friend class UnitSystem;
};


//! Unit aggregation utility
/*!
 *  Provides functions for updating and calculating effective units from a stack
 *  of unit-exponent pairs. Matching units are aggregated, where a standard unit
 *  simplifies access when joining exponents. The utility is used in the context
 *  of effective reaction rate units.
 *
 * @note Helper utility class for internal use within Cantera.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 * @ingroup unitsGroup
 */
struct UnitStack
{
    UnitStack(const Units& standardUnits) {
        stack.reserve(2); // covers memory requirements for most applications
        stack.emplace_back(standardUnits, 0.);
    }

    //! Alternative constructor allows for direct assignment of vector
    UnitStack(std::initializer_list<pair<Units, CanteraDouble>> units)
        : stack(units) {}

    UnitStack() = default;

    //! Size of UnitStack
    size_t size() const { return stack.size(); }

    //! Get standard unit used by UnitStack
    Units standardUnits() const;

    //! Set standard units
    void setStandardUnits(Units& standardUnits);

    //! Effective exponent of standard unit
    CanteraDouble standardExponent() const;

    //! Join (update) exponent of standard units, where the updated exponent is
    //! the sum of the pre-existing exponent and the exponent passed as the argument.
    void join(CanteraDouble exponent);

    //! Update exponent of item with matching units; if it does not exist,
    //! add unit-exponent pair at end of stack
    void update(const Units& units, CanteraDouble exponent);

    //! Calculate product of units-exponent stack
    Units product() const;

    vector<pair<Units, CanteraDouble>> stack; //!< Stack uses vector of pairs
};


//! Unit conversion utility
/*!
 * Provides functions for converting dimensional values from a given unit system.
 * The main use is for converting values specified in input files to Cantera's
 * native unit system, which is SI units except for the use of kmol as the base
 * unit of quantity, that is, kilogram, meter, second, kelvin, ampere, and kmol.
 *
 * String representations of units can be written using multiplication,
 * division, and exponentiation. Spaces are ignored. Positive, negative, and
 * decimal exponents are permitted. Examples:
 *
 *     kg*m/s^2
 *     J/kmol
 *     m*s^-2
 *     J/kg/K
 *
 * Metric prefixes are recognized for all units, such as nm, hPa, mg, EJ, mL, kcal.
 *
 * Special functions for converting activation energies allow these values to be
 * expressed as either energy per quantity, energy (for example, eV), or temperature by
 * applying a factor of the Avogadro number or the gas constant where needed.
 *
 * @ingroup unitsGroup
 */
class UnitSystem
{
public:
    //! Create a unit system with the specified default units
    UnitSystem(std::initializer_list<string> units);

    //! Default constructor for unit system (needed as VS2019 does not
    //! recognize an optional argument with a default value)
    UnitSystem() : UnitSystem({}) {}

    //! Return default units used by the unit system
    map<string, string> defaults() const;

    //! Set the default units to convert from when explicit units are not
    //! provided. Defaults can be set for mass, length, time, quantity, energy,
    //! and pressure. Conversion using the pressure or energy units is done only
    //! when the target units explicitly contain pressure or energy units.
    //!
    //! * To use SI+kmol: `setDefaults({"kg", "m", "s", "Pa", "J", "kmol"});`
    //! * To use CGS+mol: `setDefaults({"cm", "g", "dyn/cm^2", "erg", "mol"});`
    void setDefaults(std::initializer_list<string> units);

    //! Set the default units using a map of dimension to unit pairs.
    //!
    //! Defaults for dimensions not specified will be left unchanged. To use
    //! Cantera's default units:
    //! ```
    //! UnitSystem system;
    //! map<string, string> defaults{
    //!     {"length", "m"}, {"mass", "kg"}, {"time", "s"},
    //!     {"quantity", "kmol"}, {"pressure", "Pa"}, {"energy", "J"},
    //!     {"activation-energy", "J/kmol"}
    //! };
    //! setDefaults(defaults);
    //! ```
    void setDefaults(const map<string, string>& units);

    //! Set the default units to convert from when using the
    //! `convertActivationEnergy` function.
    void setDefaultActivationEnergy(const string& e_units);

    //! Convert `value` from the units of `src` to the units of `dest`.
    CanteraDouble convert(CanteraDouble value, const string& src, const string& dest) const;
    CanteraDouble convert(CanteraDouble value, const Units& src, const Units& dest) const;

    //! Convert `value` to the specified `dest` units from the appropriate units
    //! for this unit system (defined by `setDefaults`)
    CanteraDouble convertTo(CanteraDouble value, const string& dest) const;
    CanteraDouble convertTo(CanteraDouble value, const Units& dest) const;

    //! Convert `value` from the specified `src` units to units appropriate for
    //! this unit system (defined by `setDefaults`)
    CanteraDouble convertFrom(CanteraDouble value, const string& src) const;
    CanteraDouble convertFrom(CanteraDouble value, const Units& src) const;

    //! Convert a generic AnyValue node to the units specified in `dest`. If the
    //! input is a CanteraDouble, convert it using the default units. If the input is a
    //! string, treat this as a dimensioned value, such as '988 kg/m^3' and convert
    //! from the specified units.
    CanteraDouble convert(const AnyValue& val, const string& dest) const;
    CanteraDouble convert(const AnyValue& val, const Units& dest) const;

    //! Convert a generic AnyValue node representing a reaction rate coefficient to the
    //! units specified in `dest`. Works like `convert(AnyValue&, Units&)` but with
    //! special handling for the case where the destination units are undefined.
    //!
    //! @since New in %Cantera 3.0
    CanteraDouble convertRateCoeff(const AnyValue& val, const Units& dest) const;

    //! Convert an array of AnyValue nodes to the units specified in `dest`. For
    //! each node, if the value is a CanteraDouble, convert it using the default units,
    //! and if it is a string, treat it as a value with the given dimensions.
    vector<CanteraDouble> convert(const vector<AnyValue>& vals, const string& dest) const;
    vector<CanteraDouble> convert(const vector<AnyValue>& vals, const Units& dest) const;

    //! Convert `value` from the units of `src` to the units of `dest`, allowing
    //! for the different dimensions that can be used for activation energies
    CanteraDouble convertActivationEnergy(CanteraDouble value, const string& src,
                                   const string& dest) const;

    //! Convert `value` to the units specified by `dest` from the default
    //! activation energy units
    CanteraDouble convertActivationEnergyTo(CanteraDouble value, const string& dest) const;
    CanteraDouble convertActivationEnergyTo(CanteraDouble value, const Units& dest) const;

    //! Convert `value` from the units specified by `src` to the default
    //! activation energy units
    CanteraDouble convertActivationEnergyFrom(CanteraDouble value, const string& src) const;

    //! Convert a generic AnyValue node to the units specified in `dest`. If the
    //! input is a CanteraDouble, convert it using the default units. If the input is a
    //! string, treat this as a dimensioned value, such as '2.7e4 J/kmol', and
    //! convert from the specified units.
    CanteraDouble convertActivationEnergy(const AnyValue& val, const string& dest) const;

    //! Get the changes to the defaults from `other` to this UnitSystem
    AnyMap getDelta(const UnitSystem& other) const;

private:
    //! Factor to convert mass from this unit system to kg
    CanteraDouble m_mass_factor = 1.0;

    //! Factor to convert length from this unit system to meters
    CanteraDouble m_length_factor = 1.0;

    //! Factor to convert time from this unit system to seconds
    CanteraDouble m_time_factor = 1.0;

    //! Factor to convert pressure from this unit system to Pa
    CanteraDouble m_pressure_factor = 1.0;

    //! Factor to convert energy from this unit system to J
    CanteraDouble m_energy_factor = 1.0;

    //! Factor to convert activation energy from this unit system to J/kmol
    CanteraDouble m_activation_energy_factor = 1.0;

    //! Factor to convert quantity from this unit system to kmol
    CanteraDouble m_quantity_factor = 1.0;

    //! True if activation energy units are set explicitly, rather than as a
    //! combination of energy and quantity units
    bool m_explicit_activation_energy = false;

    //! Map of dimensions (mass, length, etc.) to names of specified default
    //! units
    map<string, string> m_defaults;
};

}

#endif
