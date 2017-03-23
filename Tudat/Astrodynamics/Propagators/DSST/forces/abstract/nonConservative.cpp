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

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void NonConservative::updateMembers( )
{
    using namespace mathematical_constants;

    determineIntegrationLimits();

    // Determine number of nodes for the Gaussian quadrature
    N = std::max( std::ceil( maximumScalableNumberOfQuadratureNodes * ( L2 - L1 ) / ( 2 * PI ) ),
                  double( fixedNumberOfQuadratureNodes ) );
    // std::cout << N << std::endl;
}


//! Update the acceleration model to the current epoch and get the perturbing acceleration
Eigen::Vector3d NonConservative::getPerturbingAcceleration()
{
    // Create current global Cartesian state from r and v
    Eigen::Matrix< double, Eigen::Dynamic, 1 > cartesianState( 6 );
    cartesianState << r, v;
    // std::cout << r.norm() << std::endl;
    // std::cout << v.norm() << std::endl;
    // std::cout << "e = " << orbital_element_conversions::convertCartesianToKeplerianElements( ( Eigen::Vector6d() << r, v ).finished(), aux.mu )( 1 ) << std::endl;
    // std::cout << EquinoctialElements::fromCartesian( cartesianState, aux.mu, orbital_element_conversions::meanType ).toKeplerian().transpose() << std::endl;
    cartesianState += aux.centralBody->getStateInBaseFrameFromEphemeris( epoch );
    // std::cout << cartesianState.transpose() << std::endl;
    // std::cout << "My altitude (substep) = " << ( r.norm() - aux.centralBody->getShapeModel()->getAverageRadius() ) / 1e3 << " km" << std::endl;
    // std::cout << aux.equinoctialElements.toKeplerian().transpose() << std::endl;

    // std::cout << std::setprecision( 12 ) << aux.centralGravityAM->getCurrentPositionOfBodyExertingAcceleration().transpose() << std::endl;
    // std::cout << std::setprecision( 12 ) << aux.centralGravityAM->getCurrentPositionOfBodySubjectToAcceleration().transpose() << std::endl;

    // Update environment
    environmentUpdater->updateEnvironment( epoch, { { transational_state, cartesianState } } );
    accelerationModel->updateMembers( epoch );
    aux.centralGravityAM->updateBaseMembers();

    // std::cout << std::setprecision( 12 ) << aux.centralGravityAM->getCurrentPositionOfBodyExertingAcceleration().transpose() << std::endl;
    // std::cout << std::setprecision( 12 ) << aux.centralGravityAM->getCurrentPositionOfBodySubjectToAcceleration().transpose() << std::endl;

    // Get acceleration
    // std::cout << "a = " << accelerationModel->getAcceleration().norm() << " m/s^2\n" << std::endl;
    return accelerationModel->getAcceleration();
}


//! Update r, v, X, Y, Xdot and Ydot and partial derivatvies of the equinoctial elements wrt to v.
Eigen::Vector6d NonConservative::integrand( const double trueLongitude )
{
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
    // std::cout << std::setprecision( 10 ) << epoch << std::endl;

    // Compute other things...
    const double sinL = std::sin( L );
    const double cosL = std::cos( L );

    // R  [ Eq. 2.1.4-(6) ]
    R = a * B * B / ( 1 + h * sinL + k * cosL );

    // X, Y  [ Eq. 2.1.4-(7) ]
    X = R * cosL;
    Y = R * sinL;

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

    // std::cout << aux.epoch << ": r = \n" << r.norm() << "\n" << std::endl;

    // V
    V = v.norm();

    // Partial derivatives of the equinoctial elements wrt to v  [ Eq. 2.1.7-(3) ]
    const Vector3d W = ( I * q * Y - p * X ) / A * w;
    partials.col( 0 ) = 2 / ( std::pow( meanMotion, 2 ) * a ) * v;                                      // da/dv
    partials.col( 1 ) = ( ( 2 * Xdot * Y - X * Ydot ) * f - X * Xdot * g ) / mu + k / B * W;            // dh/dv
    partials.col( 2 ) = ( ( 2 * X * Ydot - Xdot * Y ) * g - Y * Ydot * f ) / mu - h / B * W;            // dk/dv
    partials.col( 3 ) =     C * Y / ( 2 * A * B ) * w;                                                  // dp/dv
    partials.col( 4 ) = I * C * X / ( 2 * A * B ) * w;                                                  // dq/dv
    partials.col( 5 ) = -2 / A * r + ( k * partials.col( 1 ) - h * partials.col( 2 ) ) / ( 1 + B ) + W; // dlambda/dv

    // Disturbing acceleration
    perturbingAcceleration = getPerturbingAcceleration();

    // Compute vector of integrands:               6x3 * 3x1                    ->  6x1  ( Vector6d )
    return std::pow( R / a, 2 ) * partials.transpose() * perturbingAcceleration;
}


//! Get the mean element rates for the current auxiliary elements [ Eq. 3.4-(1b) ]
Eigen::Vector6d NonConservative::computeMeanElementRates( )
{
    // std::cout << "My altitude = " << ( aux.equinoctialElements.toCartesian( aux.mu ).segment( 0, 3 ).norm() - aux.centralBody->getShapeModel()->getAverageRadius() ) / 1e3 << " km" << std::endl;

    using namespace Eigen;
    using namespace numerical_quadrature;
    using namespace mathematical_constants;

    if ( L1 >= L2 ) { // Return immediately if drag can be completely neglected
        return Vector6d::Zero();
    }

    GaussianQuadrature< double, Vector6d > integral( boost::bind( &NonConservative::integrand, this, _1 ), L1, L2, N );
    Vector6d quadratureResult = integral.getQuadrature();
    Vector6d integralResult = quadratureResult / ( 2 * PI * B );
    // std::cout << "Mean element rates due to drag: " << integralReulst.transpose() << std::endl;


    // FIXME: is this necessary?
    /*
    // Reset central body and environment to current actual state

    // Create current global Cartesian state from r and v
    Eigen::Matrix< double, Eigen::Dynamic, 1 > cartesianState( 6 );
    cartesianState << aux.equinoctialElements.toCartesian( aux.mu );
    cartesianState += aux.centralBody->getStateInBaseFrameFromEphemeris( aux.epoch );

    // Update environment
    environmentUpdater->updateEnvironment( aux.epoch, { { transational_state, cartesianState } } );
    // accelerationModel->updateMembers( aux.epoch );
    */
    return integralResult;
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

} // namespace sst

} // namespace propagators

} // namespace tudat
