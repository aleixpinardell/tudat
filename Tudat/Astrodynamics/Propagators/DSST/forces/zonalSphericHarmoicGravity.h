#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_ZONALSPHERICHARMONICGRAVITY_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_ZONALSPHERICHARMONICGRAVITY_H

#include "abstract/sphericHarmoicGravity.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsFactories.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Final class for the contribution of the zonal terms of the spherical harmonic expansion of a gravity
class ZonalSphericalHarmonicGravity final : public SphericalHarmonicGravity
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     */
    ZonalSphericalHarmonicGravity(
            AuxiliaryElements &auxiliaryElements, SphericalHarmonicsAM sphericalHarmonicsAM )
        :
        ForceModel( auxiliaryElements ),
        SphericalHarmonicGravity( auxiliaryElements, sphericalHarmonicsAM )
    {
        if ( shOrder > 0 ) {
            throw std::runtime_error( "DSST propagator does not support tesseral terms." );
        }
    }


private:

    //! Function to compute the U derivatives for the current state [ Eq. 3.2-(2) ]
    Eigen::Vector6d computeMeanDisturbingFunctionPartialDerivatives( );

    //! Get the short period terms for the current auxiliary elements [ Sections 2.5.2 and 4.2 ]
    Eigen::Vector6d computeShortPeriodTerms( );


    //! Qns coefficients generator
    coefficients_factories::QnsCoefficientsFactory QnsFactory;

    //! Gs coefficients generator
    coefficients_factories::GsHsCoefficientsFactory GsFactory;

    //! R/a up to the Nth power
    Vectord powR_a;

};


} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_ZONALSPHERICHARMONICGRAVITY_H
