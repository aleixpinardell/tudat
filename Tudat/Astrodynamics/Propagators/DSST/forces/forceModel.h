#ifndef TUDAT_DSST_FORCEMODEL_H
#define TUDAT_DSST_FORCEMODEL_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/auxiliaryElements.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/shortPeriodTerms.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{


class DSSTForceModel {
    // Destructor
    // virtual ~DSSTForceModel() {}

    virtual std::vector< ShortPeriodTerms > initialize( AuxiliaryElements auxiliaryElements, const bool meanOnly )
    {
        return std::vector< ShortPeriodTerms >();
    }

    virtual void initializeStep( AuxiliaryElements auxiliaryElements ) { }

    virtual v_double getMeanElementRates( ) { return v_double(); }

    virtual void updateShortPeriodTerms( std::vector< SpacecraftState > meanStates );
};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_FORCEMODEL_H
