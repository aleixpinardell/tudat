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

#include <chrono>
#include <limits>
#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>
#include <eigen/Eigen/Eigen>

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

#include "Tudat/Astrodynamics/Propagators/DSST/forces/availableForceModels.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/osculatingToMean.h"


namespace tudat
{

namespace unit_tests
{

double bodyMass() {
    return 3000;
}

double bodyArea() {
    return 15;
}

double bodyCD() {
    return 2.2;
}

double bodyCR() {
    return 1.5;
}


Eigen::Vector3d earthPosition() {
    return ( Eigen::Vector3d() << 150e9, 0, 0 ).finished();
}

double earthGM() {
    return 3.9860044e14;
}

double earthMeanRadius() {
    return 6371e3;
}

Eigen::Vector3d earthRotationalVelocity() {
    return ( Eigen::Vector3d() << 0, 0, 2 * mathematical_constants::PI / physical_constants::JULIAN_DAY ).finished();
}

Eigen::Vector3d sunPosition() {
    return ( Eigen::Vector3d() << 1e4, -1e3, 1e2 ).finished();
}

double sunGM() {
    return 1.327e20;
}


BOOST_AUTO_TEST_SUITE( test_dsst_propagator )

BOOST_AUTO_TEST_CASE( dsst_propagation )
{
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace mathematical_constants;

    using namespace aerodynamics;
    using namespace electro_magnetism;
    using namespace simulation_setup;

    using namespace propagators::dsst;
    using namespace propagators::dsst::force_models;

    // auto tolerance = 100 * std::numeric_limits< double >::epsilon( );


    double initialEpoch = 7.9025e7;  // 4 Jul 2002

    // Transform to equinoctial
    Eigen::Vector6d keplerian;
    keplerian << 24288e3, 0.7293, 0.1744, 0.0226, 3.3657, 2.9189;
    EquinoctialElements equinoctial = EquinoctialElements::fromKeplerian( keplerian, eccentricType );


    // Create bodies

    // SPACECRAFT
    boost::shared_ptr< Body > propagatedBody = boost::make_shared< Body >( "Spacecraft" );
    propagatedBody->defineMass( bodyMass );
    propagatedBody->defineCrossSectionalArea( bodyArea );
    propagatedBody->defineDragCoefficient( bodyCD );
    // propagatedBody->defineRadiationPressureCoefficient( bodyCR );

    // EARTH
    boost::shared_ptr< SphericalCelestialBody > centralBody = boost::make_shared< SphericalCelestialBody >( "Earth" );
    centralBody->definePosition( earthPosition );
    centralBody->defineGravitationalParameter( earthGM );
    centralBody->defineRadius( earthMeanRadius );
    centralBody->defineRotationalVelocity( earthRotationalVelocity );

    // SUN
    boost::shared_ptr< CelestialBody > thirdBody = boost::make_shared< CelestialBody >( "Sun" );
    thirdBody->definePosition( sunPosition );
    thirdBody->defineGravitationalParameter( sunGM );


    // Model GGM02C
    const double R = 6378136.30;                          // reference radius for the geopotential expansion [m]
    std::vector< double > Js = { -4.8416938905481E-04,    // J2
                                  9.5718508415439E-07,    // J3
                                  5.3999143526074E-07,    // J4
                                  6.8715981001179E-08,    // J5
                                 -1.4994011973125E-07,    // J6
                                  9.0504630229132E-08 };  // J7

    // Create auxiliary elements (used throughout the propagation by all force models)
    AuxiliaryElements aux( propagatedBody, centralBody );
    aux.updateState( initialEpoch, equinoctial );

    // Acceleration models from Tudat
    std::vector< AvailableAcceleration > accelerationModels
            = { spherical_harmonic_gravity, third_body_central_gravity, aerodynamic, cannon_ball_radiation_pressure };

    // Create all force models
    std::vector< boost::shared_ptr< ForceModel > > forceModels;
    for ( auto accelerationModel : accelerationModels ) { // for each perturbation
        switch ( accelerationModel ) {
        case spherical_harmonic_gravity: {
            ZonalSphericalHarmonicGravity forceModel( aux, R, Js );
            forceModels.push_back( boost::make_shared< ZonalSphericalHarmonicGravity >( forceModel ) );
            break;
        }
        case third_body_central_gravity: {
            ThirdBodyCentralGravity forceModel( aux, thirdBody );
            forceModels.push_back( boost::make_shared< ThirdBodyCentralGravity >( forceModel ) );
            break;
        }
        case aerodynamic: {
            const double altitudeLimit = 800e3;
            auto atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );
            auto atmosphereModel = createAtmosphereModel( atmosphereSettings, "Earth" );
            AtmosphericDrag forceModel( aux, atmosphereModel, altitudeLimit );
            forceModels.push_back( boost::make_shared< AtmosphericDrag >( forceModel ) );
            break;
        }
        case cannon_ball_radiation_pressure: {
            auto radiationPressureInterface = boost::make_shared< RadiationPressureInterface >(
                        boost::lambda::constant( 3.839E26 ), sunPosition, earthPosition, bodyCR(), bodyArea() );
            ConservativeRadiationPressure forceModel( aux, thirdBody, radiationPressureInterface );
            forceModels.push_back( boost::make_shared< ConservativeRadiationPressure >( forceModel ) );
            break;
        }
        default:
            throw std::runtime_error( "The acceleration model is not supported by the  propagator." );
            break;
        }
    }

    // Transform aux to mean elements before starting with the propagation
    element_conversions::transformOsculatingToMeanElements( aux, forceModels );



    // Start integration
    auto t_ini = std::chrono::steady_clock::now();

    std::vector< double > integrationEpochs = { initialEpoch };
    std::map< double, Eigen::Vector6d > meanElementRates;
    for ( double epoch : integrationEpochs ) { // for each integration step
        Eigen::Vector6d meanElementRate = Eigen::Vector6d::Zero();
        for ( auto forceModel : forceModels ) {
            Eigen::Vector6d Ai = forceModel->getMeanElementRates();
            meanElementRate += Ai;
            std::cout << Ai << "\n" << std::endl;
            // std::cout << Ai.norm() << "\n" << std::endl;
        }

        // Add mean motion to fast variable derivative
        meanElementRate( fastVariableIndex ) += aux.n;

        // Store mean element rate
        meanElementRates[ epoch ] = meanElementRate;

        // Now we have k1 = A(t) for the RK4 integrator (step-size ∆)
        // We need to compute k2, k3 and k4, for different a's and t's
        // E.g. for k2:
        // // update aux elements -> a(t) + ∆/2*k1
        // // update aux epoch    -> t+∆/2
        // // compute meanElements only (no short period terms needed)
        // // k2 = A( a(t) + ∆/2*k1, t+∆/2 )
        // Same for k3, k4
        // Then we can integrate: a(t+∆) = f(a(t),k1,k2,k3,k4,∆)
        // Update aux with new integrated state
        // Compute and store shortPeriodTerms
    }
    std::cout << meanElementRates[ initialEpoch ] << "\n" << std::endl;

    auto t_end = std::chrono::steady_clock::now();

    double computationTime =
            std::chrono::duration_cast< std::chrono::microseconds >( t_end - t_ini ).count() / 1e3;

    std::cout << "Took " << computationTime << " ms." << std::endl;

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat







