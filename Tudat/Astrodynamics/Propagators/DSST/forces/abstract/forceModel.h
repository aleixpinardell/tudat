#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_FORCEMODEL_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_FORCEMODEL_H

// #include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
// #include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/vectors.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/auxiliaryElements.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{

//! Base abstract class for all force models in T
class ForceModel {
public:

    Eigen::Vector6d getMeanElementRates( ) {
        prepareForComputating();
        meanElementRatesUpToDate = true;
        return computeMeanElementRates();
    }

    Eigen::Vector6d getShortPeriodTerms( ) {
        if ( computingShortPeriodTermsRequiresMeanElementRatesToBeComputedFirst() && ! meanElementRatesUpToDate ) {
            getMeanElementRates();
        } else {
            prepareForComputating();
        }
        return computeShortPeriodTerms();
    }


protected:

    ForceModel( AuxiliaryElements &auxiliaryElements ) :
        aux               ( auxiliaryElements                     ),
        I                 ( auxiliaryElements.I                   ),
        a                 ( auxiliaryElements.sma                 ),
        h                 ( auxiliaryElements.h                   ),
        k                 ( auxiliaryElements.k                   ),
        p                 ( auxiliaryElements.p                   ),
        q                 ( auxiliaryElements.q                   ),
        ecc               ( auxiliaryElements.ecc                 ),
        A                 ( auxiliaryElements.A                   ),
        B                 ( auxiliaryElements.B                   ),
        C                 ( auxiliaryElements.C                   ),
        meanMotion        ( auxiliaryElements.n                   ),
        lastUpdateElements( auxiliaryElements.equinoctialElements ) { }


    AuxiliaryElements &aux;

    //! Retrograde factor I
    const int &I;

    // Equinoctial elements (according to T notation)
    //! a
    const double &a;
    //! e_y
    const double &h;
    //! e_x
    const double &k;
    //! h_y
    const double &p;
    //! h_x
    const double &q;

    //! Eccentricity
    const double &ecc;

    // Common factors for potential computation
    //! A = n * a²
    const double &A;
    //! B = sqrt(1 - e²)
    const double &B;
    //! C = 1 + p² + q²
    const double &C;

    //! The mean motion (n)
    const double &meanMotion;

    // Direction cosines of the symmetry axis.
    // Calculation depends on whether the perturbation is caused by a central body or a third body.
    //! α
    double alpha;
    //! β
    double beta;
    //! γ
    double gamma;


private:

    //! Whether the force model has been set up.
    bool setupCompleted = false;

    //! Set up the force model.
    //! To be implemented by derived classes. Must call updateMembers() at least once.
    //! Can include any code that needs to be run only once during the instance's lifetime.
    virtual void setUp() = 0;

    EquinoctialElements lastUpdateElements;

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( ) = 0;

    //! Set up the force model if it has not been set up yet, or update to the current auxiliary elements if needed.
    void prepareForComputating() {
        if ( ! setupCompleted ) {
            setUp();
            setupCompleted = true;
        } else if ( aux.equinoctialElements != lastUpdateElements ) {
            updateMembers();
            lastUpdateElements = aux.equinoctialElements;
        }
    }

    //! Whether the mean element rates have been computed for the current auxiliary elements.
    bool meanElementRatesUpToDate = false;

    //! Returns the mean element rates and, optionally, updates properties used later by computeShortPeriodTerms().
    virtual Eigen::Vector6d computeMeanElementRates( ) = 0;

    //! Whether computing the short period terms requires the mean element rates to be computed first.
    virtual bool computingShortPeriodTermsRequiresMeanElementRatesToBeComputedFirst() {
        return true;
    }

    //! Returns the short period terms.
    virtual Eigen::Vector6d computeShortPeriodTerms( ) = 0;


};



} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_FORCEMODEL_H
