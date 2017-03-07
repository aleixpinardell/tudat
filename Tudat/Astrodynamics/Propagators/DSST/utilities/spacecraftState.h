#ifndef TUDAT_DSST_SPACECRAFTSTATE_H
#define TUDAT_DSST_SPACECRAFTSTATE_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/forceModel.h"
#include "orbit.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


class SpacecraftState {
public:
    Orbit orbit;

    Eigen::Vector6d equinoctialComponents;

    SpacecraftState(const Orbit orbit) : orbit(orbit) { }

    AbsoluteDate getDate() const {
        return orbit.getDate();
    }

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_SPACECRAFTSTATE_H
