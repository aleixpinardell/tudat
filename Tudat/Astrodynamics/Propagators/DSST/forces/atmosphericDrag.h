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

namespace sst
{

namespace force_models
{


typedef boost::shared_ptr< aerodynamics::AerodynamicAcceleration > AerodynamicAM;


//! Abstract class for perturbations that can be expressed as a disturbing potential
class AtmosphericDrag final : public CentralBodyPerturbed, public NonConservative {
public:

    AtmosphericDrag( AuxiliaryElements &auxiliaryElements, const std::string &perturbingBody,
                     AerodynamicAM aerodynamicAM, const double maximumAltitude = TUDAT_NAN,
                     const unsigned int maximumScalableNumberOfQuadratureNodes = 48,
                     const unsigned int fixedNumberOfQuadratureNodes = 0 ) :
        ForceModel( auxiliaryElements ),
        CentralBodyPerturbed( auxiliaryElements ),
        NonConservative( auxiliaryElements, perturbingBody, aerodynamicAM,
                         maximumScalableNumberOfQuadratureNodes, fixedNumberOfQuadratureNodes )
    {
        // FIXME: user should be able to specify this altitude limit somehow during simulation setup
        hmax = maximumAltitude;
        if ( hmax != hmax ) {  // hmax is not a number, i.e. user has not defined it
            if ( perturbingBody == "Earth" ) {
                hmax = 600e3;   // drag is only considered below 600 km altitude
            } else {
                hmax = 0.0;     // 0 -> drag is considered for any altitude
            }
        }
    }


private:

    //! Set up the force model.
    void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    //! Distance limit for consideration of drag [ m ]
    double hmax;


    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    void determineIntegrationLimits( );


};



} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H
