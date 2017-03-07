#ifndef TUDAT_DSST_ZONALSPHERICHARMONICGRAVITY_H
#define TUDAT_DSST_ZONALSPHERICHARMONICGRAVITY_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "abstract/sphericHarmoicGravity.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


//! Final class for the contribution of the zonal terms of the spherical harmonic expansion of a gravity
class ZonalSphericalHarmonicGravityForceModel final : public SphericalHarmonicGravityForceModel
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     */
    ZonalSphericalHarmonicGravityForceModel( AuxiliaryElements &auxiliaryElements, double referenceRadius,
                                       std::vector< double > &Jterms ) :
        ForceModel( auxiliaryElements ),
        SphericalHarmonicGravityForceModel( auxiliaryElements, referenceRadius, Jterms ) { }


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


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_ZONALSPHERICHARMONICGRAVITY_H
