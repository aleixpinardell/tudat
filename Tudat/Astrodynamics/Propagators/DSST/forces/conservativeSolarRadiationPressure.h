#ifndef TUDAT_DSST_CONSERVATIVESOLARRADIATIONPRESSURE_H
#define TUDAT_DSST_CONSERVATIVESOLARRADIATIONPRESSURE_H

// #include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "abstract/conservativeThirdBodyPertubed.h"

// #include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


//! Final class for the contribution of a third body's central gravity
class DSSTConservativeSolarRadiationPressureForceModel final : public DSSTConservativeThirdBodyPerturbedForceModel
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     * \param thirdBody The third body exerting the acceleration.
     */
    DSSTConservativeSolarRadiationPressureForceModel( AuxiliaryElements &auxiliaryElements, Body &thirdBody ) :
        DSSTForceModel( auxiliaryElements ),
        DSSTConservativeThirdBodyPerturbedForceModel( auxiliaryElements, thirdBody ) { }


private:

    double computeUfactor() {
	return ;
    }

    //! Minimum value for n in the expansion
    unsigned int nMin( unsigned int s ) {
	return std::max( 1, (int) s );
    }

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_CONSERVATIVERADIATIONPRESSURE_H
