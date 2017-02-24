#ifndef TUDAT_DSST_ORBIT_H
#define TUDAT_DSST_ORBIT_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/absoluteDate.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

class Orbit {

private:
    /** Frame in which are defined the orbital parameters. */
    /* FIXME
    const Frame frame;
    */

    /** Date of the orbital parameters. */
    const AbsoluteDate date;

    /** Value of mu used to compute position and velocity (m³/s²). */
    const double mu;

protected:
    Orbit( /*const Frame frame,*/ const AbsoluteDate date, const double gravitationalParameter ) :
        /*frame(frame),*/ date(date), mu(gravitationalParameter) { }

public:
    AbsoluteDate getDate() const {
        return date;
    }

    double getCentralBodyGravitatinalParameter() const {
        return mu;
    }

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_ORBIT_H
