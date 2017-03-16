#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_NONCONSERVATIVE_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_NONCONSERVATIVE_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


typedef boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > AccelerationModel;


//! Abstract class for perturbations that can be expressed as a disturbing potential
class NonConservative : public virtual ForceModel {
protected:

    //! Derived (non-abstract) classes must call this constructor with either
    //! maximumScalableNumberOfQuadratureAbscissae or fixedNumberOfQuadratureAbscissae ≥ 2
    NonConservative( AuxiliaryElements &auxiliaryElements, AccelerationModel accelerationModel,
                     const unsigned int maximumScalableNumberOfQuadratureAbscissae,
                     const unsigned int fixedNumberOfQuadratureAbscissae = 0 ) :
        ForceModel( auxiliaryElements ),
        accelerationModel( accelerationModel ),
        maximumScalableNumberOfQuadratureAbscissae( maximumScalableNumberOfQuadratureAbscissae ),
        fixedNumberOfQuadratureAbscissae( fixedNumberOfQuadratureAbscissae ) { }


    //! Pointer to the dynamic simulator
    // FIXME

    //! Pointer to the corresponding acceleration model
    AccelerationModel accelerationModel;

    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    // Common factors for integral computation
    //! b = 1 / ( 1 + B )
    // double b;


    //! Lower limit for the true longitude in the averaging integral
    double L1 = -mathematical_constants::PI;

    //! Upper limit for the true longitude in the averaging integral ( L2 - L1 ≤ 2π )
    double L2 = mathematical_constants::PI;


    //! The value of the true longitude during the current integration step
    double L;

    //! The epoch during the current integration step
    double epoch;

    //! The propagated body distance |x, y, z| during the current integration step
    double R;

    //! The propagated body speed |vx, vy, vz| during the current integration step
    double V;

    //! The x components in the equinoctial reference frame during the current integration step
    double X;

    //! The y components in the equinoctial reference frame during the current integration step
    double Y;

    //! The vx components in the equinoctial reference frame during the current integration step
    double Xdot;

    //! The vy components in the equinoctial reference frame during the current integration step
    double Ydot;

    //! The propagated body position (x, y, z) during the current integration step
    Eigen::Vector3d r;

    //! The propagated body velocity (vx, vy, vz) during the current integration step
    Eigen::Vector3d v;

    //! Partial derivatives of the equinoctial elements wrt v (/dvx, /dvy, /dvz) during the current integration step
    Eigen::Matrix< double, 3, 6 > partials;

    //! The perturbing acceleration q (qx, qy, qz) during the current integration step
    Eigen::Vector3d perturbingAcceleration;


private:

    //! Number of steps for the numerical quadrature of the averaging integral when the limits of the integral
    //! differ by 2π.
    const unsigned int maximumScalableNumberOfQuadratureAbscissae;

    //! Fixed number of steps for the numerical quadrature of the averaging integral, regardless of integral limits.
    const unsigned int fixedNumberOfQuadratureAbscissae;

    //! The largerst of fixedNumberOfQuadratureAbscissae and maximumScalableNumberOfQuadratureAbscissae*√[(L2-L1)/2π]
    unsigned int N() const {
        using namespace mathematical_constants;
        return std::max( std::ceil( maximumScalableNumberOfQuadratureAbscissae * std::sqrt( ( L2 - L1 ) / ( 2*PI ) ) ),
                         double( fixedNumberOfQuadratureAbscissae ) );
    }

    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    virtual void determineIntegrationLimits( ) = 0;

    //! Update r, v, X, Y, Xdot and Ydot and partial derivatvies of the equinoctial elements wrt to v.
    void updateCartesianRelatedMembers( );

    //! Get the value of the integrand for a given `trueLongitude` (for each of the equinoctial elements)
    Eigen::Vector6d integrand( const double trueLongitude );

    //! Set the values of the minimum and maximum true longitude for the averaging integral.
    virtual Eigen::Vector3d getDisturbingAcceleration( ) {
        // FIXME: update dynamic simulator to epoch -> update environment...
        return accelerationModel->getAcceleration();
    }

    //! Get the mean element rates for the current auxiliary elements [ Eq. 3.1-(1) ]
    Eigen::Vector6d computeMeanElementRates( );

    //! Returns the short period terms.
    Eigen::Vector6d computeShortPeriodTerms( );

};



} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_NONCONSERVATIVE_H
