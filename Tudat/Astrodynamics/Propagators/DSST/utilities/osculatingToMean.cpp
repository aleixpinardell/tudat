/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "osculatingToMean.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

//! Transform osculating to mean elements.
void transformOsculatingToMeanElements( AuxiliaryElements &auxiliaryElements,
                    std::vector< boost::shared_ptr< ForceModel > > &forceModels )
{
    using namespace tudat::mathematical_constants;
    using namespace tudat::orbital_element_conversions;

    // Store initial osculating elements
    Eigen::Vector6d osculating = auxiliaryElements.equinoctialElements;

    // Reference to mean elements
    Eigen::Vector6d &mean = auxiliaryElements.equinoctialElements;

    // threshold for each parameter
    Eigen::Vector6d keplerian = convertEquinoctialToKeplerianElements(
                osculating, auxiliaryElements.mu(), auxiliaryElements.retrogradeElements );
    const double epsilon    = 1.0e-13;
    const double thresholdA = epsilon * ( 1 + std::fabs( keplerian( semiMajorAxisIndex ) ) );
    const double thresholdE = epsilon * ( 1 + keplerian( eccentricityIndex ) );
    const double thresholdI = epsilon * ( 1 + keplerian( inclinationIndex ) );
    const double thresholdL = epsilon * PI;
    Eigen::Vector6d thresholds;
    thresholds << thresholdA, thresholdE, thresholdE, thresholdI, thresholdI, thresholdL;

    unsigned int i = 0;
    const unsigned int maxIterations = 200;
    while ( i++ < maxIterations ) {
        // Initialize rebuilt to be equal to the current mean elements
        Eigen::Vector6d rebuilt = mean;

        // Add short period terms to rebuilt to make it osculating
        for ( auto forceModel : forceModels ) {
            rebuilt += forceModel->getShortPeriodTerms();
        }

        // Get difference between osculating and rebuilt
        Eigen::Vector6d delta = osculating - rebuilt;

        // Convert last element to range [-π, π]
        delta( 5 ) -= 2*PI * std::floor( 0.5 * ( delta( 5 ) / PI + 1 ) );

        // Check convergence
        if ( ( delta.array().abs() < thresholds.array() ).all() ) {
            return;
        }

        // Update mean elements and auxiliaryElements' members
        mean += delta;
        auxiliaryElements.updateMembers();
    }

    throw std::runtime_error( "Impossible to transform osculating to mean elements: "
                              "maximum number of iterations (200) exceeded." );
}




} // namespace dsst

} // namespace propagators

} // namespace tudat
