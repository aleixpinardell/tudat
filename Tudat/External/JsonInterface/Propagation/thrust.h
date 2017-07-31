/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_THRUST_H
#define TUDAT_JSONINTERFACE_THRUST_H

#include <Tudat/SimulationSetup/PropagationSetup/thrustSettings.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

/// DIRECTION

static std::map< ThrustDirectionGuidanceTypes, std::string > thrustDirectionTypes =
{
    { colinear_with_state_segment_thrust_direction, "colinearWithStateSegment" },
    { thrust_direction_from_existing_body_orientation, "fromExistingBodyOrientation" },
    { custom_thrust_direction, "customDirection" },
    { custom_thrust_orientation, "customOrientation" }
};

//! `ThrustDirectionGuidanceTypes` not supported by `json_interface`.
static std::vector< ThrustDirectionGuidanceTypes > unsupportedThrustDirectionTypes =
{
    custom_thrust_direction,
    custom_thrust_orientation
};

//! Convert `ThrustDirectionGuidanceTypes` to `json`.
inline void to_json( json& jsonObject, const ThrustDirectionGuidanceTypes& directionType )
{
    jsonObject = json_interface::stringFromEnum( directionType, thrustDirectionTypes );
}

//! Convert `json` to `AvailableAcceleration`.
inline void from_json( const json& jsonObject, ThrustDirectionGuidanceTypes& directionType )
{
    directionType = json_interface::enumFromString( jsonObject.get< std::string >( ), thrustDirectionTypes );
}

//! Create a `json` object from a shared pointer to a `ThrustDirectionGuidanceSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings );

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings );


/// MAGNITUDE

//! Map of `ThrustMagnitudeTypes` string representations.
static std::map< ThrustMagnitudeTypes, std::string > thrustMagnitudeTypes =
{
    { constant_thrust_magnitude, "constant" },
    { from_engine_properties_thrust_magnitude, "fromEngineProperties" },
    { thrust_magnitude_from_time_function, "timeDependent" },
    { thrust_magnitude_from_dependent_variables, "variableDependent" }
};

//! `ThrustMagnitudeTypes` not supported by `json_interface`.
static std::vector< ThrustMagnitudeTypes > unsupportedThrustMagnitudeTypes =
{
    thrust_magnitude_from_time_function,
    thrust_magnitude_from_dependent_variables
};

//! Convert `ThrustMagnitudeTypes` to `json`.
inline void to_json( json& jsonObject, const ThrustMagnitudeTypes& magnitudeType )
{
    jsonObject = json_interface::stringFromEnum( magnitudeType, thrustMagnitudeTypes );
}

//! Convert `json` to `AvailableAcceleration`.
inline void from_json( const json& jsonObject, ThrustMagnitudeTypes& magnitudeType )
{
    magnitudeType = json_interface::enumFromString( jsonObject.get< std::string >( ), thrustMagnitudeTypes );
}

//! Create a `json` object from a shared pointer to a `ThrustEngineSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ThrustEngineSettings >& magnitudeSettings );

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ThrustEngineSettings >& magnitudeSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_THRUST_H
