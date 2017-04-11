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


//! Struct for storing settings for AtmosphericDrag ForceModel
struct AtmosphericDragSettings: NonConservativeSettings {
    //! Default constructor
    AtmosphericDragSettings( double altitudeLimit                   = 0.0,
                             unsigned int numberOfQuadratureNodes   = 40,
                             bool numberOfQuadratureNodesIsScalable = true,
                             bool alwaysIncludeCentralNode          = true )
        : NonConservativeSettings( numberOfQuadratureNodes, numberOfQuadratureNodesIsScalable,
                                   alwaysIncludeCentralNode ),
          altitudeLimit( altitudeLimit ) { }

    //! Altitude limit above which atmospheric drag is completely neglected.
    //! If set to 0, drag is considered for all altitudes.
    double altitudeLimit;
};


//! Abstract class for perturbations that can be expressed as a disturbing potential
class AtmosphericDrag final : public CentralBodyPerturbed, public NonConservative {
public:

    AtmosphericDrag( AuxiliaryElements &auxiliaryElements, const std::string &perturbingBody,
                     AerodynamicAM aerodynamicAM, boost::shared_ptr< AtmosphericDragSettings > settings = NULL ) :
        ForceModel( auxiliaryElements, settings ),
        CentralBodyPerturbed( auxiliaryElements ),
        NonConservative( auxiliaryElements, perturbingBody, aerodynamicAM, settings ) { }

    boost::shared_ptr< AtmosphericDragSettings > getSettings( ) {
        boost::shared_ptr< AtmosphericDragSettings > castedSettings =
                boost::dynamic_pointer_cast< AtmosphericDragSettings >( settings );
        return castedSettings != NULL ? castedSettings : boost::make_shared< AtmosphericDragSettings >( );
    }


private:

    //! Set up the force model.
    void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );

    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    void determineIntegrationLimits( );


};



} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_ATMOSPHERICDRAG_H
