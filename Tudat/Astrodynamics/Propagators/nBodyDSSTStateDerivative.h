/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NBODYCDSSTSTATEDERIVATIVE_H
#define TUDAT_NBODYCDSSTSTATEDERIVATIVE_H

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

#include "Tudat/Astrodynamics/Propagators/DSST/forces/availableForceModels.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/osculatingToMean.h"

namespace tudat
{

namespace propagators
{

//! Class for computing the state derivative of translational motion of N bodies, using a DSST propagator.
template< typename StateScalarType = double, typename TimeType = double >
class NBodyDSSTStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyDSSTStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                              const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                              const std::vector< std::string >& bodiesToIntegrate,
                              const Eigen::Matrix< StateScalarType, 6, 1 >& initialCartesianElements,
                              const TimeType& initialEpoch ) :
        NBodyStateDerivative< StateScalarType, TimeType >( accelerationModelsPerBody, centralBodyData,
                                                           dsst, bodiesToIntegrate )
    {
        using namespace basic_astrodynamics;
        using namespace orbital_element_conversions;

        using namespace aerodynamics;
        using namespace electro_magnetism;
        using namespace gravitation;

        using namespace sst;
        using namespace sst::force_models;

        // Check only one propagated body
        if ( bodiesToIntegrate.size() > 1 )
        {
            throw std::runtime_error( "DSST propagator does not support simultaneous propagation of multiple bodies." );
        }

        // Get acceleration models for propagated body
        const std::string propagatedBody = *bodiesToIntegrate.begin();
        const SingleBodyAccelerationMap accelerationModelsMap = accelerationModelsPerBody.at( propagatedBody );

        // First iteration through acceleration models, only to find central gravity model
        CentralGravityAM centralGravityAM;
        for ( auto ent: accelerationModelsMap )
        {
            const std::string bodyExertingAccelerations = ent.first;
            const std::vector< AccelerationModel3dPointer > accelerationModels = ent.second;
            for ( auto accelerationModel: accelerationModels )
            {
                centralGravityAM = boost::dynamic_pointer_cast<
                        SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > >( accelerationModel );
                if ( centralGravityAM != NULL ) { break; }
            }
            if ( centralGravityAM != NULL ) { break; }
        }

        if ( centralGravityAM == NULL ) {
            throw std::runtime_error( "Could not find a central or spherical harmonics gravity model required "
                                      "for setting up the DSST propagator." );
        }

        // Convert from Cartesian to (eccentric) equinoctial elements
        EquinoctialElements equinoctialElements = EquinoctialElements::fromCartesian(
                    initialCartesianElements.template cast< double >(),
                    centralGravityAM->getGravitationalParameterFunction()(), eccentricType );

        // Create auxiliary elements
        auxiliaryElements = AuxiliaryElements( initialEpoch, equinoctialElements, centralGravityAM );

        // Iterate through all acceleration models to create the DSST force models
        for ( auto ent: accelerationModelsMap )
        {
            const std::string bodyExertingAccelerations = ent.first;
            const std::vector< AccelerationModel3dPointer > accelerationModels = ent.second;
            for ( AccelerationModel3dPointer accelerationModel: accelerationModels )
            {
                AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );
                switch ( accelerationType ) {
                case central_gravity:
                {
                    // No DSST force model required
                    break;
                }
                case spherical_harmonic_gravity:
                {
                    // std::cout << "spherical_harmonic_gravity" << std::endl;
                    SphericalHarmonicsAM sphericalHarmonicsAM = boost::dynamic_pointer_cast<
                            SphericalHarmonicsGravitationalAccelerationModel >( accelerationModel );
                    forceModels.push_back( boost::make_shared< ZonalSphericalHarmonicGravity >(
                                               auxiliaryElements, sphericalHarmonicsAM ) );
                    break;
                }
                case third_body_central_gravity:
                {
                    // std::cout << "third_body_central_gravity" << std::endl;
                    ThirdBodyAM thirdBodyAM = boost::dynamic_pointer_cast<
                            ThirdBodyCentralGravityAcceleration >( accelerationModel );
                    forceModels.push_back( boost::make_shared< ThirdBodyCentralGravity >(
                                               auxiliaryElements, thirdBodyAM ) );
                    break;
                }
                case aerodynamic:
                {
                    // std::cout << "aerodynamic" << std::endl;
                    AerodynamicAM aerodynamicAM = boost::dynamic_pointer_cast<
                            AerodynamicAcceleration >( accelerationModel );
                    forceModels.push_back( boost::make_shared< AtmosphericDrag >(
                                               auxiliaryElements, aerodynamicAM ) );
                    break;
                }
                case cannon_ball_radiation_pressure:
                {
                    // std::cout << "cannon_ball_radiation_pressure" << std::endl;
                    RadiationPressureAM radiationPressureAM = boost::dynamic_pointer_cast<
                            CannonBallRadiationPressureAcceleration >( accelerationModel );
                    if ( radiationPressureAM->getRadiationPressureInterface()->getOccultingBodyRadii().size() > 0 ) {
                        throw std::runtime_error( "DSST propagator does not support consideration of occultations "
                                                  "during the evaluation of solar radiation pressure." );
                    } else {
                        forceModels.push_back( boost::make_shared< ConservativeRadiationPressure >(
                                                   auxiliaryElements, radiationPressureAM ) );
                    }
                    break;
                }
                default:
                {
                    throw std::runtime_error( "Error, DSST propagator does not support accelerations of type: " +
                                              boost::lexical_cast< std::string >( accelerationType ) );
                    break;
                }
                }
            }
        }

        // Convert initial state from osculating to mean elements
        // FIXME: Currently, short-period terms are implemented as vectors of zeros, so this does virtually nothing
        // element_conversions::transformOsculatingToMeanElements( auxiliaryElements, forceModels );
    }


    //! Destructor
    ~NBodyDSSTStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system.
    /*!
     * Calculates the state derivative (velocity+acceleration of each body) of the translational motion of the system
     * at the given time and position/velocity of bodies.
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing Caartesian
     *  position/velocity of the bodies being integrated. The order of the values is defined by the order of bodies in
     *  bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current state derivative (velocity+acceleration) of system of bodies integrated numerically
     *  (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        using namespace Eigen;
        using namespace orbital_element_conversions;
        using namespace sst;

        Eigen::Vector6d equinoctialComponents = stateOfSystemToBeIntegrated.template cast< double >();
        auxiliaryElements.updateState( time, equinoctialComponents );

        Vector6d meanElementRates = Vector6d::Zero();

        for ( auto forceModel: forceModels ) {
            Eigen::Vector6d Ai = forceModel->getMeanElementRates();
            meanElementRates += Ai;
            // std::cout << Ai << "\n" << std::endl;
        }
        // std::cout << "meanElementRates :\n" << meanElementRates << "\n" << std::endl;

        // Add mean motion to the mean longitude rate
        meanElementRates( fastVariableIndex ) += auxiliaryElements.n;

        // Update state derivative (only one body being propagated, thus only first column needs to be updated)
        stateDerivative.col( 0 ) = meanElementRates.template cast< StateScalarType >();
    }

    //! Function to convert the state in the conventional form to the DSST propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the Cowell propagator,
     * the two are equivalent, and this function returns the input state.
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form' (which is equal to outputSolution).
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        using namespace Eigen;
        using namespace orbital_element_conversions;
        using namespace sst;

        Matrix< StateScalarType, Dynamic, Dynamic > dsstSolution;
        dsstSolution.resizeLike( cartesianSolution );

        // std::cout << "cartesianSolution: " << cartesianSolution << std::endl;

        for ( unsigned int i = 0; i < cartesianSolution.cols(); i++ )
        {
            Vector6d cartesianComponents = cartesianSolution.col( i ).template cast< double >();

            EquinoctialElements equinoctialElements =
                    EquinoctialElements::fromCartesian( cartesianComponents, auxiliaryElements.mu, meanType );

            dsstSolution.col( i ) = equinoctialElements.getComponents().template cast< StateScalarType >();
        }

        // std::cout << "dsstSolution: " << dsstSolution << std::endl;

        return dsstSolution;
    }

    //! Function to convert the DSST propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. For the Cowell propagator,
     * the two are equivalent, and this function returns the input state.
     * In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration, equal to conventional form for this class).
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentCartesianLocalSoluton State (internalSolution), converted to the 'conventional form',
     *  which is equal to outputSolution for this class (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
            const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        using namespace Eigen;
        using namespace orbital_element_conversions;
        using namespace sst;

        // std::cout << "internalSolution: " << internalSolution << std::endl;

        Vector6d equinoctialComponents = internalSolution.template cast< double >();

        EquinoctialElements equinoctialElements( equinoctialComponents, meanType,
                                                 auxiliaryElements.equinoctialElements.isRetrograde() );

        Vector6d cartesian = equinoctialElements.toCartesian( auxiliaryElements.mu );

        currentCartesianLocalSoluton = cartesian.template cast< StateScalarType >();
        // std::cout << "currentCartesianLocalSoluton: \n" << cartesian.head(3).norm() << "\n\n" << std::endl;

    }


private:

    //! Auxiliary elements shared amongst all the force models
    sst::AuxiliaryElements auxiliaryElements;

    //! List of DSST force models
    std::vector< boost::shared_ptr< sst::force_models::ForceModel > > forceModels;

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NBODYCDSSTSTATEDERIVATIVE_H
