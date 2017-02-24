#ifndef TUDAT_DSST_AUXILIARYELEMENTS_H
#define TUDAT_DSST_AUXILIARYELEMENTS_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "absoluteDate.h"
#include "spacecraftState.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

class AuxiliaryElements
{
public:
    //! Constructor.
    /** Simple constructor.
     * \param equinoctialElements A vector containing the current equinoctial elements.
     * \param centralBodyGravitationalParameter G*M of the central body [m^3/s^2].
     * \param retrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
     * with inclinations of 180 degrees). Default is false, which avoids singularities for equatorial orbits.
     */
    AuxiliaryElements( const SpacecraftState state, const bool retrogradeElements );

    //! Orbit date
    AbsoluteDate date;

    //! Central body gravitational parameter
    double mu;

    //! Eccentricity
    double ecc;

    //! Keplerian mean motion
    double n;

    //! Keplerian period
    double period;

    //! Semi-major axis
    double sma;

    //! x component of eccentricity vector
    double k;

    //! y component of eccentricity vector
    double h;

    //! x component of inclination vector
    double q;

    //! y component of inclination vector
    double p;

    //! Mean longitude
    double lm;

    //! True longitude
    double lv;

    //! Eccentric longitude
    double le;

    //! Retrograde factor I
    int    I;

    //! A = sqrt(μ * a)
    double A;

    //! B = sqrt(1 - h² - k²)
    double B;

    //! C = 1 + p² + q²
    double C;

    //! Equinoctial frame f vector
    Vector3 f;

    //! Equinoctial frame g vector
    Vector3 g;

    //! Equinoctial frame w vector
    Vector3 w;

    //! Direction cosine α
    double alpha;

    //! Direction cosine β
    double beta;

    //! Direction cosine γ
    double gamma;
};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_AUXILIARYELEMENTS_H
