#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"

// #include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

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


typedef boost::shared_ptr< aerodynamics::AerodynamicAcceleration > AerodynamicAM;


//! Abstract class for perturbations that can be expressed as a disturbing potential
class AtmosphericDrag final : public CentralBodyPerturbed, public NonConservative {
public:

    AtmosphericDrag( AuxiliaryElements &auxiliaryElements, AerodynamicAM aerodynamicAM, double maximumAltitude = 0.0,
                     const unsigned int maximumScalableNumberOfQuadratureAbscissae = 48,
                     const unsigned int fixedNumberOfQuadratureAbscissae = 0 ) :
      ForceModel( auxiliaryElements ),
      CentralBodyPerturbed( auxiliaryElements ),
      NonConservative( auxiliaryElements, aerodynamicAM,
                       maximumScalableNumberOfQuadratureAbscissae, fixedNumberOfQuadratureAbscissae ),
      hmax( maximumAltitude )
    {
        // R = auxiliaryElements.centralGravityAM->get
    }


private:

    //! Pointer to aerodynamic acceleration model
    AerodynamicAM aerodynamicAM;

    //! Set up the force model.
    void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    //! Radius of the central body [ m ]
    double R = 6371e3;  // FIXME

    //! Distance limit for consideration of drag [ m ]
    const double hmax;


    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    void determineIntegrationLimits( );


};



} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H
