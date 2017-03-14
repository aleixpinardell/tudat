/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "nonConservative.h"

#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"

#include "boost/bind.hpp"


namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void NonConservative::updateMembers( )
{
    // b = 1 / ( 1 + B );

    determineIntegrationLimits();
}


//! Update r, v, X, Y, Xdot and Ydot and partial derivatvies of the equinoctial elements wrt to v.
Eigen::Vector6d NonConservative::integrand( const double trueLongitude ) {
    using namespace Eigen;
    using namespace mathematical_constants;
    using namespace orbital_element_conversions;

    // True longitude
    L = trueLongitude;

    // Get epoch
    Vector6d components = aux.equinoctialElements.getComponents( meanType );
    const double auxMeanLongitude = components( fastVariableIndex );
    components( fastVariableIndex ) = L;
    const double meanLongitude = getMeanLongitudeFromTrue( components );
    epoch = aux.epoch + ( meanLongitude - auxMeanLongitude ) / meanMotion;

    // Change of variable: https://pomax.github.io/bezierinfo/legendre-gauss.html [Eq. 2]
    // L = 0.5 * ( ( L2 - L1 ) * L + L2 + L1 );

    // Compute other things...
    const double sinL = std::sin( L );
    const double cosL = std::cos( L );

    // R  [ Eq. 2.1.4-(6) ]
    R = a * B * B / ( 1 + h * sinL + k * cosL );

    // X, Y  [ Eq. 2.1.4-(7) ]
    X = R * sinL;
    Y = R * cosL;

    // Xdot, Ydot  [ Eq. 2.1.4-(8) ]
    Xdot = -meanMotion * a / B * ( h + sinL );
    Ydot =  meanMotion * a / B * ( k + cosL );

    // Get f, g, w, µ
    Vector3d &f = aux.f;
    Vector3d &g = aux.g;
    Vector3d &w = aux.w;
    const double mu = aux.mu;

    // r, v  [ Eq. 2.1.4-(9) ]
    r = X    * f + Y    * g;
    v = Xdot * f + Ydot * g;

    // V
    V = v.norm();

    // Partial derivatives of the equinoctial elements wrt to v  [ Eq. 2.1.7-(3) ]
    const Vector3d W = ( I * q * Y - p * X ) / B * w;
    partials.col( 0 ) = 2 / ( std::pow( meanMotion, 2 ) * a ) * v;                                      // da/dv
    partials.col( 1 ) = ( ( 2 * Xdot * Y - X * Ydot ) * f - X * Xdot * g ) / mu + k / B * W;            // dh/dv
    partials.col( 2 ) = ( ( 2 * X * Ydot - Xdot * Y ) * g - Y * Ydot * f ) / mu - h / B * W;            // dk/dv
    partials.col( 3 ) =     C * Y / ( 2 * A * B ) * w;                                                  // dp/dv
    partials.col( 4 ) = I * C * X / ( 2 * A * B ) * w;                                                  // dq/dv
    partials.col( 5 ) = -2 / A * r + ( k * partials.col( 1 ) - h * partials.col( 2 ) ) / ( 1 + B ) + W; // dlambda/dv

    // Disturbing acceleration
    perturbingAcceleration = getDisturbingAcceleration();

    // Compute vector of integrands
    return std::pow( R / a, 2 ) * partials.transpose() * perturbingAcceleration;
}


//! Get the mean element rates for the current auxiliary elements [ Eq. 3.4-(1b) ]
Eigen::Vector6d NonConservative::computeMeanElementRates( )
{
    using namespace Eigen;
    using namespace numerical_quadrature;
    using namespace mathematical_constants;

    if ( L1 >= L2 ) { // Return immediately if drag can be completely neglected
        return Vector6d::Zero();
    }

    GaussianQuadrature< double, Vector6d > integral( boost::bind( &NonConservative::integrand, this, _1 ),
                                                     L1, L2, N() );
    return 1 / ( 2 * PI * B ) * integral.getQuadrature();
}


//! Get the short period terms for the current auxiliary elements [ Section 4.4 ]
Eigen::Vector6d NonConservative::computeShortPeriodTerms( )
{
    // Numerical quadrature
    // GaussianQuadrature gaussSin( integrandSin, L1, L2, N );
    // GaussianQuadrature gaussCos( integrandCos, L1, L2, N );

    return Eigen::Vector6d::Zero();
}


} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat
