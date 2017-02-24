/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "coefficientsfactory.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

//! Initializer.
CoefficientsFactory::CoefficientsFactory() {
    VNS[ NSKey(0, 0) ] = 1.0;
    VNS[ NSKey(1, 0) ] = 0.0;
    VNS[ NSKey(1, 1) ] = 0.5;
}

//! Function to compute the dK_0^{n,s}/d\chi coefficients for a given \chi from the recurrence formula 3.2-(3)
NSSetPair CoefficientsFactory::computeK0nsAndDerivatives( const double chi, const int nMax, const int sMax )
{
    // Initialization
    vv_double dK0ns( nMax + 1 );
    for ( int i = 0; i <= nMax; i++ )
    {
        dK0ns[ i ] = std::vector< double >( std::min( i + 2, sMax + 1 ) );
    }

    // first elements
    // dK0ns[0][0] = 1;
    // dK0ns[0][1] = -1;

    // Get Kns coefficients
    vv_double K0ns = computeK0ns( chi, nMax, sMax );

    // Compute recursively
    for ( double n = 1; n <= nMax; n++ )
    {
        for ( double s = 0; s < dK0ns[n].size(); s++ )
        {
            if ( n == s || n == s - 1 )
            {
                dK0ns[n][s] = 0;
            }
            else if ( n > s && n > 2 )
            {
                dK0ns[n][s] = (2*n + 1) / (n + 1) * dK0ns[n - 1][s] -
                        (n + s) * (n - s) / ( n * (n + 1) * chi*chi ) * dK0ns[n - 2][s] +
                        2*(n + 2)*(n - s) / ( n * (n + 1) * chi*chi*chi ) * K0ns[n - 2][s];
            }
        }
    }

    return NSSetPair( K0ns , dK0ns );
}


//! Function to compute the K_0^{n,s} coefficients for a given \chi from the recurrence formula 2.3.7-(7)
vv_double CoefficientsFactory::computeK0ns( const double chi, const int nMax, const int sMax )
{
    // Initialization
    vv_double K0ns( nMax + 1 );
    for ( int i = 0; i <= nMax; i++ )
    {
        K0ns[ i ] = new_v_double( std::min( i + 2, sMax + 1 ) );
    }

    // first elements
    K0ns[0][0] = 1;
    K0ns[0][1] = -1;

    // Compute recursively
    for ( double n = 1; n <= nMax; n++ )
    {
        for ( double s = 0; s < K0ns[n].size(); s++ )
        {
            if ( n == s && n >= 1 )
            {
                K0ns[n][s] = (2*s + 1) / (s + 1) * K0ns[s - 1][s];
            }
            else if ( n == s - 1 && n >= 1 )
            {
                K0ns[n][s] = (1 - 2*s) / s * K0ns[s - 2][s - 1];
            }
            else if ( n >= s + 1 && s + 1 >= 2 )
            {
                K0ns[n][s] = (2*n + 1) / (n + 1) * K0ns[n - 1][s] -
                        (n + s) * (n - s) / ( n * (n + 1) * chi*chi ) * K0ns[n - 2][s];
            }
        }
    }

    return K0ns;
}


//! Function to compute the Q_{n,s} coefficients evaluated at \gamma from the recurrence formula 2.8.3-(2)
vv_double CoefficientsFactory::computeQns( const double gamma, const int nMax, const int sMax )
{
    // Initialization
    const int sDim = std::min( sMax + 1, nMax ) + 1;
    const int rows = nMax + 1;
    vv_double Qns( rows );
    for ( int i = 0; i <= nMax; i++ )
    {
        const int snDim = std::min( i + 1, sDim );
        Qns[ i ] = new_v_double( snDim );
    }

    // first element
    Qns[0][0] = 1;

    for ( int n = 1; n <= nMax; n++ )
    {
        const int snDim = std::min( n + 1, sDim );
        for ( int s = 0; s < snDim; s++ )
        {
            if ( n == s )
            {
                Qns[n][s] = (2.0 * s - 1.0) * Qns[s - 1][s - 1];
            }
            else if ( n == (s + 1) )
            {
                Qns[n][s] = (2.0 * s + 1.0) * gamma * Qns[s][s];
            }
            else
            {
                Qns[n][s] = (2.0 * n - 1.0) * gamma * Qns[n - 1][s] - (n + s - 1.) * Qns[n - 2][s];
                Qns[n][s] /= n - s;
            }
        }
    }

    return Qns;
}


//! Function to compute the dQ_{n,s}/d\gamma coefficients evaluated at \gamma from the recurrence formula 3.1-(8)
vv_double CoefficientsFactory::computeQnsDerivatives( const double gamma, const int nMax, const int sMax )
{
    // Initialization
    const int sDim = std::min( sMax + 1, nMax ) + 1;
    const int rows = nMax + 1;
    vv_double dQns( rows );
    for ( int i = 0; i <= nMax; i++ )
    {
        const int snDim = std::min( i + 1, sDim );
        dQns[ i ] = new_v_double( snDim );
    }

    // Get Qns until (n,s+1)
    vv_double Qns = computeQns( gamma, nMax, sMax + 1 );

    for ( int n = 0; n <= nMax; n++ )
    {
        const int snDim = std::min( n + 1, sDim );
        for ( int s = 0; s < snDim; s++ )
        {
            if ( n == s )
            {
                dQns[n][s] = 0;
            }
            else
            {
                dQns[n][s] = Qns[n][s + 1];
            }
        }
    }

    return dQns;
}



//! Function to compute recursively G_s and H_s polynomials from equation 3.1-(5)
Eigen::Matrix2Xd CoefficientsFactory::computeGsHs(
        const double k, const double h, const double alpha, const double beta, const int order )
{
    // Constant terms
    const double hamkb = h * alpha - k * beta;
    const double kaphb = k * alpha + h * beta;

    // Initialization
    Eigen::Matrix2Xd GsHs( 2, order + 1 );
    GsHs(0,0) = 1.0;
    GsHs(1,0) = 0.0;

    for ( int s = 1; s <= order; s++ ) {
        // Gs coefficient
        GsHs(0,s) = kaphb * GsHs(0, s - 1) - hamkb * GsHs(1, s - 1);
        // Hs coefficient
        GsHs(1,s) = hamkb * GsHs(0, s - 1) + kaphb * GsHs(1, s - 1);
    }

    return GsHs;
}


//! Function to compute derivatives of G_s with respect to k, h, alpha and beta from equation 3.1-(9)
std::map< std::string, Eigen::VectorXd > CoefficientsFactory::computeGsDerivatives(
        const double k, const double h, const double alpha, const double beta, const int order )
{
    // Initialization
    std::map< std::string, Eigen::VectorXd > dGs;
    std::vector< std::string > keys = { "h", "k", "alpha", "beta" };
    for ( unsigned int i = 0; i < keys.size(); i++ ) {
        dGs[ keys[ i ] ] = Eigen::VectorXd( order + 1 );
        dGs[ keys[ i ] ]( 0 ) = 0.0;
    }

    Eigen::Matrix2Xd GsHs = computeGsHs( k, h, alpha, beta, order );
    Eigen::VectorXd  Gs   = GsHs.row( 0 );
    Eigen::VectorXd  Hs   = GsHs.row( 1 );

    for ( double s = 1; s <= order; s++ ) {
        dGs[ "h" ]( s )      = s * beta  * Gs(s - 1) - s * alpha * Hs(s - 1);
        dGs[ "k" ]( s )      = s * alpha * Gs(s - 1) + s * beta  * Hs(s - 1);
        dGs[ "alpha" ]( s )  = s * k     * Gs(s - 1) - s * h     * Hs(s - 1);
        dGs[ "beta" ]( s )   = s * h     * Gs(s - 1) + s * k     * Hs(s - 1);
    }

    return dGs;
}


//! Function to compute the V_{n,s} coefficients from Eq. 2.8.2-(1)(2).
NSMap CoefficientsFactory::computeVns( const int order )
{
    if ( order > LAST_VNS_ORDER )
    {
        const int min = (LAST_VNS_ORDER - 2 < 0) ? 0 : (LAST_VNS_ORDER - 2);
        for ( int n = min; n < order; n++ )
        {
            for ( int s = 0; s < n + 1; s++ )
            {
                if ( (n - s) % 2 != 0 )
                {
                    VNS[ NSKey(n, s) ] = 0.0;
                }
                else
                {
                    // s = n
                    if ( n == s && (s + 1) < order )
                    {
                        VNS[ NSKey(s + 1, s + 1) ] = VNS[ NSKey(s, s) ] / (2 * s + 2);
                    }
                    // otherwise
                    if ( (n + 2) < order )
                    {
                        VNS[ NSKey(n + 2, s)] = VNS[ NSKey(n, s) ] * (-n + s - 1) / (n + s + 2);
                    }
                }
            }
        }
        LAST_VNS_ORDER = order;
    }
    return VNS;
}

//! Function to get a V_{n,s}^m coefficient from Eq. 2.7.1-(6) and Sec. 2.8.2
double CoefficientsFactory::getVmns( const int m, const int n, const int s )
{
    if ( m > n )
    {
        std::cerr << "Cannot compute the Vmns coefficient with m > n (" << m << " > " << n << ")" << std::endl;
    }
    const double fns = factorial( n + std::abs( s ) );
    const double fnm = factorial( n  - m );

    double result = 0;
    // If (n - s) is odd, the Vmsn coefficient is null
    if ( ( n - s ) % 2 == 0 ) {
        // Update the Vns coefficients
        if ( ( n + 1 ) > LAST_VNS_ORDER )
        {
            computeVns( n + 1 );
        }
        if ( s >= 0 )
        {
            result = fns * VNS[ NSKey(n, s) ] / fnm;
        }
        else
        {
            // If s < 0 : Vmn-s = (-1)^(-s) Vmns
            const int mops = ( s % 2 == 0 ) ? 1 : -1;
            result = mops * fns * VNS[ NSKey(n, -s) ] / fnm;
        }
    }
    return result;
}


} // namespace dsst

} // namespace propagators

} // namespace tudat
