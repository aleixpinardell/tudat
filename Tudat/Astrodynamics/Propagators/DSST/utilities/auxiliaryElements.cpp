/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "auxiliaryElements.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

//! Update the instance's members according to the current equinoctialElements.
void AuxiliaryElements::updateMembers()
{
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::mathematical_constants;
    using std::sqrt;
    using std::pow;

    // Retrograde factor
    I = retrogradeElements ? -1 : 1;

    // Get equinoctial elements
    sma = equinoctialElements[ semiMajorAxisIndex ];
    h   = equinoctialElements[ hIndex ];
    k   = equinoctialElements[ kIndex ];
    p   = equinoctialElements[ pIndex ];
    q   = equinoctialElements[ qIndex ];
    lm  = equinoctialElements[ meanLongitudeIndex ];

    // Define repeated values
    const double h2 = pow( h, 2 );
    const double k2 = pow( k, 2 );
    const double p2 = pow( p, 2 );
    const double q2 = pow( q, 2 );

    // Eccentricity
    ecc = sqrt( h2 + k2 );

    // Keplerian mean motion
    n = sqrt( mu() / pow( sma, 3 ) );

    // Keplerian period
    period = 2 * PI * sma * sqrt( sma / mu() );

    // Get eccentric and true longitudes
    le = getEccentricLongitude( equinoctialElements );
    lv = getTrueLongitude( equinoctialElements );

    // Get A, B and C
    A = sqrt( mu() * sma );
    B = sqrt( 1 - h2 - k2 );
    C = 1 + p2 + q2;

    // Get equinoctial reference frame vectors
    Eigen::Matrix< double, 3, 3 > fgw = getEquinoctialReferenceFrameBasisVectors( p, q, retrogradeElements );
    f = fgw.col( 0 );
    g = fgw.col( 1 );
    w = fgw.col( 2 );

    // Get direction cosines (for central body, not to be used for third body calculations)
    alpha = f( 2 );
    beta  = g( 2 );
    gamma = w( 2 );
}


} // namespace dsst

} // namespace propagators

} // namespace tudat
