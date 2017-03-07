#ifndef TUDAT_DSST_THIRDBODYCENTRALGRAVITY_H
#define TUDAT_DSST_THIRDBODYCENTRALGRAVITY_H

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
class ThirdBodyCentralGravityForceModel final : public ConservativeThirdBodyPerturbedForceModel
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     * \param thirdBody The third body exerting the acceleration.
     */
    ThirdBodyCentralGravityForceModel( AuxiliaryElements &auxiliaryElements, Body &thirdBody ) :
        ForceModel( auxiliaryElements ),
        ConservativeThirdBodyPerturbedForceModel( auxiliaryElements, thirdBody ) { }


private:

    double computeUfactor() {
        return mu3() / R3;
    }

    //! Minimum value for n in the expansion
    unsigned int nMin( unsigned int s ) {
        return std::max( 2, (int) s );
    }


    /*

    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// Short period terms /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    //! e / ( 1 + sqrt(1 - e²) )
     double c;

    //! Map of Jacobi polynomials factories used by getwCoefficient
    std::map< coefficients_factories::JacobiKey, coefficients_factories::JacobiPolynomialsFactory > jacobiFactories;

    //! Function to obtain the coefficients e^{-|j-s|} * w^{n,s}_j(e) and e^{-(j+s)} * w^{n,s}_{-j}(e) [Eq. 4.2-(16,10)]
    std::pair< double, double > getewjCoefficients( const int n, const int s, const int j );

    */
};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_THIRDBODYCENTRALGRAVITY_H
