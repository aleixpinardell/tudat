/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "zonalSphericalHarmonicGravity.h"
// #include "Tudat/Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

// #include "tudat/Tudat/Astrodynamics/Propagators/DSST/utilities/jacobiPolynomials.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Function to compute the U derivatives for the current state [ Eq. 3.1-(6) ]
Eigen::Vector6d ZonalSphericalHarmonicGravity::computeMeanDisturbingFunctionPartialDerivatives()
{
    using namespace coefficients_factories;

    // R / a up to power maxAR3Pow
    const double R_a = R / a;
    powR_a = Vectord( N + 1 );
    powR_a( 0 ) = 1.0;
    for ( unsigned int i = 1; i <= N; i++ ) {
        powR_a( i ) = R_a * powR_a( i - 1 );
    }

    // Initialise mean potential function
    U = 0.0;

    // Potential derivatives
    Eigen::Vector6d dU = Eigen::Vector6d::Zero();

    // Vns coefficients
    NestedVectord Vns = VnsFactory.getCoefficientsUpTo( N, S );

    // K_0{-n-1,s} coefficients and derivatives with respect to \Chi
    NegativeK0nsCoefficientsFactory mK0nsFactory( Chi, N + 1, S );
    NestedVectord  mK0ns = mK0nsFactory.getCoefficients();
    NestedVectord dmK0ns = mK0nsFactory.getDerivatives();

    // Qns coefficients and derivatives with respect to \gamma
    QnsCoefficientsFactory QnsFactory( gamma, N, S );
    NestedVectord Qns  = QnsFactory.getCoefficients();
    NestedVectord dQns = QnsFactory.getDerivatives();

    // Gs coefficients and partial derivatives with respect to h, k, \alpha and \beta
    GsHsCoefficientsFactory GsFactory( h, k, alpha, beta, S );
    Vectord Gs        = GsFactory.getGsCoefficients();
    Vectord Hs        = GsFactory.getHsCoefficients();
    Vectord dGsdh     = GsFactory.getGsHDerivatives();
    Vectord dGsdk     = GsFactory.getGsKDerivatives();
    Vectord dGsdalpha = GsFactory.getGsAlphaDerivatives();
    Vectord dGsdbeta  = GsFactory.getGsBetaDerivatives();

    for ( unsigned int s = 0; s <= S; s++ ) {
        // Get the current Gs coefficient
        const double G = Gs( s );

        // Kronecker symbol ( 2 - delta_{0,s} )
        const double delta0s = ( s == 0 ) ? 1.0 : 2.0;

        for ( unsigned int n = s + 2; n <= N; n++ ) {
            // Vns is zero for odd ( n - s ), thus only consider when ( n - s ) is even
            if ( ( n - s ) % 2 == 0 ) {
                // Coefficients for current n and s
                const double  mK0 =  mK0ns( n + 1, s );
                const double dmK0 = dmK0ns( n + 1, s );

                // Other terms
                const double R_an  = powR_a( n );
                const double coef0 = delta0s * R_an * J( n ) * Vns( n, s );
                const double coef1 = coef0 * Qns( n, s );
                const double coef2 = coef1 * mK0;
                const double coef3 = coef2 * G;

                //Compute U:
                U += coef3;

                Eigen::Vector6d dUterm;
                dUterm << coef3 * ( n + 1 ),                                // dU / da
                        coef1 * ( mK0 * dGsdh( s ) + hChi3 * G * dmK0 ),    // dU / dh
                        coef1 * ( mK0 * dGsdk( s ) + kChi3 * G * dmK0 ),    // dU / dk
                        coef2 * dGsdalpha( s ),                             // dU / dalpha
                        coef2 * dGsdbeta( s ),                              // dU / dbeta
                        coef0 * mK0 * dQns( n, s ) * G;                     // dU / dgamma

                dU += dUterm;
            }
        }
    }

    // multiply by -µ/a
    const double mu_a = -mu / a;
    U *= mu_a;
    dU *= mu_a;

    dU( 0 ) /= -a;
    return dU;
}


//! Get the short period terms for the current auxiliary elements
Eigen::Vector6d ZonalSphericalHarmonicGravity::computeShortPeriodTerms( )
{
    return Eigen::Vector6d::Zero();
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
