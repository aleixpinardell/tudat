#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_CENTRALBODYPERTURBED_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_CENTRALBODYPERTURBED_H

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Abstract class for perturbations cuased by central bodies
class CentralBodyPerturbed : public virtual ForceModel {
protected:

    CentralBodyPerturbed( AuxiliaryElements &auxiliaryElements ) :
    ForceModel ( auxiliaryElements             ),
    centralBody( auxiliaryElements.centralBody ),
    mu(          auxiliaryElements.mu          ) { }


    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );

    //! The central body exerting the acceleration
    boost::shared_ptr< CelestialBody > centralBody;

    //! Standard gravitational parameter μ for the central body in m³/s²
    double &mu;

};



} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_CENTRALBODYPERTURBED_H
