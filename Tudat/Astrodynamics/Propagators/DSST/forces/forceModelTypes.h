#ifndef TUDAT_DSST_FORCEMODELTYPES_H
#define TUDAT_DSST_FORCEMODELTYPES_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

#include "zonalSphericHarmoicGravity.h"
#include "thirdBodyCentralGravity.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

/*
//! Function to identify the derived class type of an acceleration model.
tudat::basic_astrodynamics::AvailableAcceleration getForceModelType(
        const boost::shared_ptr< ForceModel > forceModel )
{
    using namespace tudat::basic_astrodynamics;

    // Nominal type is undefined
    AvailableAcceleration accelerationType = undefined_acceleration;

    // Check for each accelerarion mdoel type implemented as AvailableAcceleration.
    if( boost::dynamic_pointer_cast< ThirdBody >( forceModel ) != NULL )
    {
        accelerationType = third_body_central_gravity;
    }
    else
    {
        throw std::runtime_error(
                    "Error, acceleration model not identified when getting acceleration type." );
    }

    // Return identified type.
    return accelerationType;

}
*/

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_FORCEMODELTYPES_H
