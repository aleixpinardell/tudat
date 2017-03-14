#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYPERTURBED_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYPERTURBED_H

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Abstract class for perturbations cuased by third bodies
class ThirdBodyPerturbed : public virtual ForceModel {
protected:

    ThirdBodyPerturbed( AuxiliaryElements &auxiliaryElements, boost::shared_ptr< CelestialBody > thirdBody ) :
    ForceModel( auxiliaryElements ), thirdBody( thirdBody ) { }


    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );

    //! The third body exerting the acceleration
    boost::shared_ptr< CelestialBody > thirdBody;

    //! Standard gravitational parameter μ for the third body in m³/s²
    double mu3;

    //! Position of the third body from the centre of mass of the central body
    Eigen::Vector3d r3;

    //! Distance from center of mass of the central body to the 3rd body
    double R3;

};



} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYPERTURBED_H
