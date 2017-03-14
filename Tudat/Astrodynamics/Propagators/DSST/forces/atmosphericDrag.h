#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

#include "Tudat/Astrodynamics/Propagators/DSST/forces/abstract/nonConservative.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/abstract/centralBodyPertubed.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Abstract class for perturbations that can be expressed as a disturbing potential
class AtmosphericDrag : public CentralBodyPerturbed, public NonConservative {
public:

    AtmosphericDrag( AuxiliaryElements &auxiliaryElements, boost::shared_ptr< aerodynamics::AtmosphereModel >
                     atmosphereModel, double maximumAltitude = 0.0,
                     const unsigned int maximumScalableNumberOfQuadratureAbscissae = 48,
                     const unsigned int fixedNumberOfQuadratureAbscissae = 0 ) :
      ForceModel( auxiliaryElements ),
      CentralBodyPerturbed( auxiliaryElements ),
      NonConservative( auxiliaryElements, maximumScalableNumberOfQuadratureAbscissae, fixedNumberOfQuadratureAbscissae ),
      atmosphereModel( atmosphereModel ), hmax( maximumAltitude )
    {
        boost::shared_ptr< SphericalCelestialBody > planet =
                boost::dynamic_pointer_cast< SphericalCelestialBody >( auxiliaryElements.centralBody );
        if ( planet == NULL ) {
            throw std::runtime_error( "Atmospheric drag can only be modelled around spherical central bodies." );
        }
        R = planet->getRadius();
    }


private:

    //! Set up the force model.
    void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );

    // Common factors for integral computation
    //! CD * A / ( 2 * m )
    double ballistic;

    //! Central's body rotational velocity
    Eigen::Vector3d omega;


    //! Reference to the atmosphere model
    boost::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel;

    //! Radius of the central body [ m ]
    double R = 0;

    //! Distance limit for consideration of drag [ m ]
    const double hmax;


    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    void determineIntegrationLimits( );

    //! Get the value of the current disturbing acceleration, to be called after base class has updated r, v and epoch.
    Eigen::Vector3d getDisturbingAcceleration();


};



} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H
