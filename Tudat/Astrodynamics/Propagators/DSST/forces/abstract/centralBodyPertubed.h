#ifndef TUDAT_DSST_CENTRALBODYPERTURBED_H
#define TUDAT_DSST_CENTRALBODYPERTURBED_H

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{


//! Abstract class for perturbations cuased by central bodies
class CentralBodyPerturbedForceModel : public virtual ForceModel {
protected:

    CentralBodyPerturbedForceModel( AuxiliaryElements &auxiliaryElements ) :
    ForceModel( auxiliaryElements             ),
    centralBody   ( auxiliaryElements.centralBody ) { }


    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );

    //! The central body exerting the acceleration
    const Body &centralBody;

    //! Standard gravitational parameter μ for the central body in m³/s²
    double mu() {
        return centralBody.getGravitationalParameter();
    }

};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_CENTRALBODYPERTURBED_H
