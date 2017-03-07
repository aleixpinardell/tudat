#ifndef TUDAT_DSST_THIRDBODYPERTURBED_H
#define TUDAT_DSST_THIRDBODYPERTURBED_H

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{


//! Abstract class for perturbations cuased by third bodies
class ThirdBodyPerturbedForceModel : public virtual ForceModel {
protected:

    ThirdBodyPerturbedForceModel( AuxiliaryElements &auxiliaryElements, Body &thirdBody ) :
    ForceModel( auxiliaryElements ),
    thirdBody     ( thirdBody	  ) { }


    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );

    //! The third body exerting the acceleration
    Body &thirdBody;

    //! Position of the third body from the centre of mass of the central body
    Eigen::Vector3d r3;

    //! Distance from center of mass of the central body to the 3rd body
    double R3;

    //! Standard gravitational parameter μ for the third body in m³/s²
    double mu3() {
        return thirdBody.getGravitationalParameter();
    }

};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_THIRDBODYPERTURBED_H
