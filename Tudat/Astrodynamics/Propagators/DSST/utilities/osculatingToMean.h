#ifndef TUDAT_DSST_OSCULATINGTOMEAN_H
#define TUDAT_DSST_OSCULATINGTOMEAN_H

#include "Tudat/Astrodynamics/Propagators/DSST/forces/abstract/forceModel.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

//! Transform osculating to mean elements.
/*!
 * Transforms the equinoctial elements contained by osculatingAuxiliaryElements from osculating to mean.
 * Based on a fixed-point interative process from Section 7.4-2.
 * \param osculatingAuxiliaryElements AuxiliaryElements containing the osculating equinoctial elements.
 * \param forceModels The force models generating the short period terms needed to do the transformation.
 */
void transformOsculatingToMeanElements( AuxiliaryElements &osculatingAuxiliaryElements,
                    std::vector< boost::shared_ptr< ForceModel > > &forceModels );


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_OSCULATINGTOMEAN_H
