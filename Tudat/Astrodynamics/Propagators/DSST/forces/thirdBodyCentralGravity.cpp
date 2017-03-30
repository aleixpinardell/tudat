/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "thirdBodyCentralGravity.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{

//! Update instance's members that are computed from the current auxiliary elements.
void ThirdBodyCentralGravity::updateMembers( )
{
    // Update vector to third body
    CentralGravityAM thirdBodyCentralGravity = thirdBodyAM->getAccelerationModelForBodyUndergoingAcceleration();

    mu3 = thirdBodyCentralGravity->getGravitationalParameterFunction()();

    r3 = thirdBodyCentralGravity->getCurrentPositionOfBodyExertingAcceleration() -
            thirdBodyCentralGravity->getCurrentPositionOfBodySubjectToAcceleration();

    // std::cout << "SUN: " << r3.transpose() << std::endl;

    ConservativeThirdBodyPerturbed::updateMembers();
}


/*

//! Get the short period terms for the current auxiliary elements [ Sections 2.5.2 and 4.2 ]
Eigen::Vector6d ThirdBodyCentralGravity::computeShortPeriodTerms( )
{
    /// update order, suborder of Vns, Qns
    /// Store Gns up to max n
    /// Cs, Ss up to max s, both for (alpha,beta) and (k,h) [ not (h,k) !!! ]
    /// getwCoefficient
    /// constant terms, such as e^...
    /// Summation -> Cj, Sj
    /// Cj, Sj -> Cij, Sij. Not in the book for third body, see Orekit.
    /// Apparently, will need derivatives of Cj, Sj wrt a, alpha, beta, gamma, h, k, lambda. Not in the book.
    /// Finally, the short period terms are given by Eq. 2.5.2-(14) as a function of Cij, Sij, eccentric longitude
    ///
    /// Expressions for Cj, Sj, Cij and Sij are unique to third body. Not necessary to create separate factories.
    /// Same as with getwCoefficient. Moreover, they need to be recalculated for every integration step, so no need
    /// to store them. Can be created and destructed inside computeShortPeriodTerms().

    const unsigned int eccPow = std::max( maxEccPow, maxEccPowShort );

    // Vns coefficients
    NestedVectord Vns = VnsFactory.getCoefficientsUpTo( maxAR3Pow, eccPow );

    // Qns coefficients
    QnsFactory.setSuborder( eccPow );

    // C^j, S^j = getFourierCjSjCoefficients()

    return Eigen::Vector6d::Zero();
}



//! Function to obtain the coefficients e^{-|j-s|} * w^{n,s}_j(e) and e^{-(j+s)} * w^{n,s}_{-j}(e) [ Eq. 4.2-(16,10) ]
std::pair< double, double > ThirdBodyCentralGravity::getewjCoefficients(
        const int n, const int s, const int j )
{
    using namespace coefficients_factories;
    using std::fabs;
    using std::pow;

    double ewpj = 0.0;
    double ewmj = 0.0;

    // Construct JacobiKey
    const double absjms = fabs( j - s );
    const double absjps = fabs( j + s );
    const JacobiKey key( absjms, absjps );

    // Create JacobiPolynomialsFactory (with initial order 0) if it does not exist
    if ( jacobiFactories.count( key ) == 0 ) {
        jacobiFactories[ key ] = JacobiPolynomialsFactory( Chi, key, 0 );
    }

    // Shared terms
    // c = e / ( 1 + sqrt(1 - e²) ) = e / ( 1 + B ) = e * b
    const double eterm = pow( ecc, -absjms );
    const double c = ecc * b;
    const double c2 = c * c;
    const double omc2 = 1 - c2;
    const double term1 = pow( omc2 / ( 1 + c2 ), n );
    const double term2 = pow( -c, absjms );
    const double term = term1 * term2;
    JacobiPolynomialsFactory &jacobiFactory = jacobiFactories[ key ];

    // Compute and return
    if( fabs( s ) > fabs( j ) ) {
        ewpj = eterm * term * factorial( n + s ) * factorial( n - s ) / ( factorial( n + j ) * factorial( n - j ) )
                * pow( omc2, -fabs( s ) ) * jacobiFactory.getPolynomial( n - fabs(s) );
    } else {
        ewpj = eterm * term * pow( omc2, -fabs( j ) ) * jacobiFactory.getPolynomial( n - fabs(j) );
    }

    // e terms

    return { ewpj, ewmj };
}

*/


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
