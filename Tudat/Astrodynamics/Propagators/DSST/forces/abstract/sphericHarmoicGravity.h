#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_SPHERICHARMONICGRAVITY_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_SPHERICHARMONICGRAVITY_H

#include "centralBodyPertubed.h"
#include "conservative.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactories.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/upperBounds.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{

//! Abstract class for gravity perturbations expressed as a spheric harmonic expansion.
//! Zonal and tesseral terms are treated differently, so derived classes will implement these terms separately.
class SphericalHarmonicGravity :
        public CentralBodyPerturbed, public Conservative
{
protected:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     */
    SphericalHarmonicGravity( AuxiliaryElements &auxiliaryElements, double referenceRadius, NestedVectord Jterms ) :
          ForceModel( auxiliaryElements	),
          CentralBodyPerturbed( auxiliaryElements ),
          Conservative( auxiliaryElements ),
          R( referenceRadius ),
          Jterms( Jterms )
    {
        shDegree = Jterms.size() + 1;
        shOrder = 0;
        for ( unsigned int n = 2; n <= shDegree; n++ ) {
            std::vector< double > terms = Jterms[ n - 2 ];
            unsigned int m = terms.size() - 1;
            if ( m > n ) {
                throw std::runtime_error( "The order cannot exceed the degree of the spheric harmonic expansion.");
            }
            if ( m > shOrder ) {
                shOrder = m;
            }
        }
    }


private:

    void setUp() {
        Conservative::setUp();
    }

    //! Update instance's members based on current auxiliary elements.
    void updateMembers( );

    //! Set the values of N and S for the series expansion of the disturbing potential.
    void determineTruncationValues( );


protected:

    // FIXME: do we also need reference µ?
    // µ = 398600.44150E+09 for GGM02C, GGM02S
    // µ = 398600.44189E+09 according to Wikipedia

    //! Reference radius used for the geopotential expansion
    const double R;

    //! Jterms from degree 2 up to N
    const NestedVectord Jterms;

    //! Returns J_{n,m}
    double J( unsigned int n, unsigned int m ) {
        return Jterms( n - 2, m );
    }

    //! Returns J_n = J_{n,0}
    double J( unsigned int n ) {
        return J( n, 0 );
    }

    //! Degree and order of the spherical harmonic expansion
    unsigned int shDegree;
    unsigned int shOrder;


    //! Truncation tolerance
    const double TRUNCATION_TOLERANCE = 1e-4;

    //! Number of points for interpolation
    // const int INTERPOLATION_POINTS = 3;

    //! Maximum power for eccentricity used in short periodic computation
    // const int MAX_ECCPOWER_SP = 4;

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

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_SPHERICHARMONICGRAVITY_H
