#ifndef TUDAT_DSST_AUXILIARYELEMENTS_H
#define TUDAT_DSST_AUXILIARYELEMENTS_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "coefficientsFactory.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

class AuxiliaryElements
{
public:
    //! Reference to the central body
    const Body &centralBody;

    //! Date, in seconds since J2000
    double epoch;

    //! Mean equinoctial elements
    Eigen::Vector6d equinoctialElements;

    //! Whether retrograde elements are being used
    bool retrogradeElements;

    //! Update the instance's members according to the current equinoctialElements.
    void updateMembers();

    //! Update to a new state.
    /*!
     * Updates instance's members using a new state.
     * \param epoch The new time in seconds since J2000.
     * \param equinoctialElements A vector containing the new equinoctial elements.
     * \param useRetrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
     * with inclinations close to 180 degrees). Default is false, which avoids singularities for equatorial orbits.
     */
    void updateState( const double epoch, const Eigen::Vector6d &equinoctialElements,
                      const bool useRetrogradeElements = false ) {
        this->epoch               = epoch;
        this->equinoctialElements = equinoctialElements;
        this->retrogradeElements  = useRetrogradeElements;
        updateMembers();
    }

    //! Constructor.
    /*!
     * Simple constructor.
     * \param centralBody A reference to the central body.
     * \param epoch The current time in seconds since J2000.
     * \param equinoctialElements A vector containing the current equinoctial elements.
     * \param useRetrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
     * with inclinations close to 180 degrees). Default is false, which avoids singularities for equatorial orbits.
     */
    AuxiliaryElements( Body &centralBody, const double epoch, const Eigen::Vector6d &equinoctialElements,
                       const bool useRetrogradeElements = false ) :
        centralBody( centralBody ) {
        updateState( epoch, equinoctialElements, useRetrogradeElements );
    }


    // Properties that depend on a referenced object (i.e. thirdBody) are not stored but accessed every time.

    //! Central body's gravitational parameter
    double mu() {
        return centralBody.getGravitationalParameter();
    }


    // Properties that depend on a stored object (i.e. equinoctialElements) are stored and updated by updateState().
    // In this way, they are only computed once for each integration step, even though they are used many times.

    //! Retrograde factor I
    int I;

    //! Eccentricity
    double ecc;

    //! Keplerian mean motion
    double n;

    //! Keplerian period
    double period;

    //! Semi-major axis
    double sma;

    //! y component of eccentricity vector
    double h;

    //! x component of eccentricity vector
    double k;

    //! y component of inclination vector
    double p;

    //! x component of inclination vector
    double q;

    //! Mean longitude
    double lm;

    //! True longitude
    double lv;

    //! Eccentric longitude
    double le;

    //! A = sqrt(μ * a)
    double A;

    //! B = sqrt(1 - h² - k²)
    double B;

    //! C = 1 + p² + q²
    double C;

    //! Equinoctial frame f vector
    Eigen::Vector3d f;

    //! Equinoctial frame g vector
    Eigen::Vector3d g;

    //! Equinoctial frame w vector
    Eigen::Vector3d w;

    //! Direction cosine α (for central body, not to be used for third body calculations)
    double alpha;

    //! Direction cosine β (for central body, not to be used for third body calculations)
    double beta;

    //! Direction cosine γ (for central body, not to be used for third body calculations)
    double gamma;

    //! Vns coefficients generator, used in several perturbations
    coefficients_factories::VnsCoefficientsFactory VnsFactory;


    //! Default assignment operator
    AuxiliaryElements & operator= ( const AuxiliaryElements & ) = default;

};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_AUXILIARYELEMENTS_H
