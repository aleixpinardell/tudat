/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN


#include <limits>
#include <boost/test/unit_test.hpp>
#include <eigen/Eigen/Eigen>

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_dsst_propagator )

BOOST_AUTO_TEST_CASE( coefficientsFactory )
{
    using namespace tudat::propagators::dsst;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::mathematical_constants;
    auto tolerance = 100 * std::numeric_limits< double >::epsilon( );

    ////////////////////////////
    /// EQUINOCTIAL ELEMENTS ///
    ////////////////////////////

    // Convert TO equinoctial elements

    const double mu = 4e14;
    Eigen::Vector6d keplerian;
    keplerian << 7e6, 0.1, 1, 2, 3, 4;
    Eigen::Vector6d cartesian = convertKeplerianToCartesianElements( keplerian, mu );

    Eigen::Vector6d equinoctialFromKep = convertKeplerianToEquinoctialElements( keplerian, mu );
    Eigen::Vector6d equinoctialFromCar = convertCartesianToEquinoctialElements( cartesian, mu );

    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( equinoctialFromKep[ i ], equinoctialFromCar[ i ], tolerance );
    }
    // std::cout << equinoctialFromKep << std::endl;
    // std::cout << equinoctialFromCar << std::endl;


    // Convert FROM equinoctial elements

    Eigen::Vector6d cartesianFromEqui = convertEquinoctialToCartesianElements( equinoctialFromCar, mu );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( cartesian[ i ], cartesianFromEqui[ i ], tolerance );
    }
//    std::cout << cartesian << std::endl;
//    std::cout << cartesianFromEqui << std::endl;

    Eigen::Vector6d keplerianFromEqui = convertEquinoctialToKeplerianElements( equinoctialFromKep, mu );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( keplerian[ i ], keplerianFromEqui[ i ], tolerance );
    }
//    std::cout << keplerian << std::endl;
//    std::cout << keplerianFromEqui << std::endl;


    // Modified equinoctial elements

    Eigen::Vector6d equinoctial = equinoctialFromKep;
    Eigen::Vector6d modEquinoctialFromEqui = convertEquinoctialToModifiedEquinoctialElements( equinoctial );
    Eigen::Vector6d modEquinoctialFromKep = convertKeplerianToModifiedEquinoctialElements( keplerian, false );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( modEquinoctialFromEqui[ i ], modEquinoctialFromKep[ i ], tolerance );
    }

    Eigen::Vector6d equinoctialFromModEqui = convertModifiedEquinoctialToEquinoctialElements( equinoctial );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( equinoctial[ i ], equinoctialFromModEqui[ i ], tolerance );
    }


    // Get eccentric and true longitudes

    const double ectrLongitude = getEccentricLongitude( equinoctial );
    const double trueLongitude = getTrueLongitude( equinoctial );
    const double meanLongitude = equinoctial[ meanLongitudeIndex ];

    const double trueAnomaly = keplerian [ trueAnomalyIndex ];
    const double e = keplerian [ eccentricityIndex ];
    const double ectrAnomaly = convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly, e );
    const double meanAnomaly = convertEllipticalEccentricAnomalyToMeanAnomaly( ectrAnomaly, e );

    BOOST_CHECK_CLOSE_FRACTION( computeModulo( trueLongitude - trueAnomaly, 2 * PI ),
                                computeModulo( ectrLongitude - ectrAnomaly, 2 * PI ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION( computeModulo( trueLongitude - trueAnomaly, 2 * PI ),
                                computeModulo( meanLongitude - meanAnomaly, 2 * PI ), tolerance );

//    std::cout << std::fmod( trueLongitude - trueAnomaly + 2 * M_PI, 2 * M_PI ) << std::endl;
//    std::cout << std::fmod( ectrLongitude - ectrAnomaly + 2 * M_PI, 2 * M_PI ) << std::endl;
//    std::cout << std::fmod( meanLongitude - meanAnomaly + 2 * M_PI, 2 * M_PI ) << std::endl;



    ////////////////////////////
    /// COEFFICIENTS FACTORY ///
    ////////////////////////////

    CoefficientsFactory coefficientsFactory;

    // K_0^{n,s} and derivatives
    NSSetPair K0nsAndDerivatives = coefficientsFactory.computeK0nsAndDerivatives( 2.1, 3, 3 );
    vv_double K0ns = K0nsAndDerivatives.first;
    vv_double dK0ns = K0nsAndDerivatives.second;
    // BOOST_CHECK_CLOSE_FRACTION( Qns[2][1], 4.5, epsilon );
//    for ( auto const &ent : K0ns ) {
//        for ( auto const &entent : ent ) {
//            std::cout << entent << ", ";
//        }
//        std::cout << "; ";
//    }
//    std::cout << std::endl;
//    for ( auto const &ent : dK0ns ) {
//        for ( auto const &entent : ent ) {
//            std::cout << entent << ", ";
//        }
//        std::cout << "; ";
//    }
//    std::cout << std::endl;

    // Q_{n,s}
    vv_double Qns = coefficientsFactory.computeQns( 1.5, 4, 4 );
    BOOST_CHECK_CLOSE_FRACTION( Qns[2][1], 4.5, tolerance );
//    for ( auto const &ent : Qns ) {
//        for ( auto const &entent : ent ) {
//            std::cout << entent << ", ";
//        }
//        std::cout << "; ";
//    }
//    std::cout << std::endl;

    // d(Q_{n,s}/d(\gamma)
    vv_double dQns = coefficientsFactory.computeQnsDerivatives( 1.5, 4, 4 );
    BOOST_CHECK_CLOSE_FRACTION( dQns[3][2], Qns[3][3], tolerance );
//    for ( auto const &ent : dQns ) {
//        for ( auto const &entent : ent ) {
//            std::cout << entent << ", ";
//        }
//        std::cout << "; ";
//    }
//    std::cout << std::endl;

    // G_s, H_s polynomials
    Eigen::Matrix2Xd GsHs = coefficientsFactory.computeGsHs( 1.0, 1.0, 1.5, 0.5, 2 );
    BOOST_CHECK_CLOSE_FRACTION( GsHs(0,2), 3.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( GsHs(1,2), 4.0, tolerance );
    // std::cout << GsHs << std::endl;

    // G_s derivatives
    std::map< std::string, Eigen::VectorXd > dGs = coefficientsFactory.computeGsDerivatives( 1.0, 1.0, 1.5, 0.5, 2 );
    BOOST_CHECK_CLOSE_FRACTION( dGs["h"](2),    -1.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( dGs["k"](2),     7.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( dGs["alpha"](2), 2.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( dGs["beta"](2),  6.0, tolerance );

    // V_{n,s} coefficients
    NSMap Vns = coefficientsFactory.computeVns(5);
    BOOST_CHECK_CLOSE_FRACTION( Vns[ NSKey(0,0) ], 1, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( Vns[ NSKey(3,1) ], -0.125, tolerance );
//    for ( auto const &ent : Vns ) {
//        std::cout << ent.first.first << "," << ent.first.second << ": " << ent.second << std::endl;
//    }

    // V_{n,s}^m coefficient
    double Vmns = coefficientsFactory.getVmns( 3, 4, 2 );
    BOOST_CHECK_CLOSE_FRACTION( Vmns, -15, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
