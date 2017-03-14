/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "atmosphericDrag.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void AtmosphericDrag::updateMembers( )
{
    CentralBodyPerturbed::updateMembers();

    ballistic = aux.CD * aux.area / ( 2 * aux.mass );
    omega = aux.centralBody->getRotationalVelocity();

    NonConservative::updateMembers();
}


//! Update the values of the minimum and maximum true longitude for the averaging integral.
void AtmosphericDrag::determineIntegrationLimits( ) {
    using namespace mathematical_constants;

    if ( hmax <= 0.0 ) { // No altitude limit specified, thus consider drag thorughout all the revolution
        L1 = -PI;
        L2 =  PI;
        return;
    }

    // Eq. 3.4-(2b)
    const double cosflim = ( a * ( 1 - std::pow( ecc, 2 ) ) / ( R + hmax ) - 1 ) / ecc;
    if ( cosflim > 1 ) { // Drag can be completely neglected
        L1 = 0;
        L2 = 0;
        return;
    } else if ( cosflim < -1 ) { // Drag has to be considered throughout all the revolution
        L1 = -PI;
        L2 =  PI;
        return;
    }
    // Otherwise, get the limits:

    const double flim = std::acos( cosflim );

    // longitude of perigee = omega + I * Omega
    // Eq. 2.1.2-(1)
    const double longitudePerigeeSin = h / ecc;
    const double longitudePerigeeCos = k / ecc;
    const double longitudePerigee = std::atan2( longitudePerigeeSin, longitudePerigeeCos );

    // Eq. 3.4-(2b)
    L1 = -flim + longitudePerigee;
    L2 =  flim + longitudePerigee;
}

//! Set the values of the minimum and maximum true longitude for the averaging integral.
Eigen::Vector3d AtmosphericDrag::getDisturbingAcceleration() {
    // From r = [x, y, z] to [h, lon, lat]
    const double h = r.norm() - R;
    const double lon = 0.0;  // FIXME
    const double lat = 0.0;  // FIXME

    Eigen::Vector3d vatm = omega.cross( r );
    Eigen::Vector3d vrel = v - vatm;
    std::cout << std::round(h/1e3) << " km: " << atmosphereModel->getDensity( h, lon, lat, epoch ) << std::endl;
    return -ballistic * atmosphereModel->getDensity( h, lon, lat, epoch ) * vrel.norm() * vrel;
}


} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat
