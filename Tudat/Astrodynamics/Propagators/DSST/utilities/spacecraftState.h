#ifndef TUDAT_DSST_SPACECRAFTSTATE_H
#define TUDAT_DSST_SPACECRAFTSTATE_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/orbit.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

class SpacecraftState {
public:
    const Orbit orbit;

    Vector6 equinoctialComponents;

    SpacecraftState(const Orbit orbit) : orbit(orbit) { }

    AbsoluteDate getDate() const {
        return orbit.getDate();
    }

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_SPACECRAFTSTATE_H
