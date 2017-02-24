#ifndef TUDAT_DSST_HERMITEINTERPOLATOR_H
#define TUDAT_DSST_HERMITEINTERPOLATOR_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


/** Polynomial interpolator using both sample values and sample derivatives.
 * <p>
 * The interpolation polynomials match all sample points, including both values
 * and provided derivatives. There is one polynomial for each component of
 * the values vector. All polynomials have the same degree. The degree of the
 * polynomials depends on the number of points and number of derivatives at each
 * point. For example the interpolation polynomials for n sample points without
 * any derivatives all have degree n-1. The interpolation polynomials for n
 * sample points with the two extreme points having value and first derivative
 * and the remaining points having value only all have degree n+1. The
 * interpolation polynomial for n sample points with value, first and second
 * derivative for all points all have degree 3n-1.
 * </p>
 *
 */
class HermiteInterpolator {
private:
    /** Sample abscissae. */
    v_double abscissae;

    /** Top diagonal of the divided differences array. */
    vv_double topDiagonal;

    /** Bottom diagonal of the divided differences array. */
    vv_double bottomDiagonal;

public:
    /** Create an empty interpolator.
     */
    // HermiteInterpolator() { }

    /** Add a sample point.
     * <p>
     * This method must be called once for each sample point. It is allowed to
     * mix some calls with values only with calls with values and first
     * derivatives.
     * </p>
     * <p>
     * The point abscissae for all calls <em>must</em> be different.
     * </p>
     * @param x abscissa of the sample point
     * @param value value and derivatives of the sample point
     * (if only one row is passed, it is the value, if two rows are
     * passed the first one is the value and the second the derivative
     * and so on)
     * @exception MathIllegalArgumentException if the abscissa difference between added point
     * and a previous point is zero (i.e. the two points are at same abscissa)
     * @exception MathRuntimeException if the number of derivatives is larger
     * than 20, which prevents computation of a factorial
     */
    void addSamplePoint(const double x, const vv_double value)
    {

        for (unsigned int i = 0; i < value.size(); ++i) {

            v_double y = value[i];
            if (i > 1) {
                double inv = 1.0 / factorial(i);
                for (unsigned int j = 0; j < y.size(); ++j) {
                    y[j] *= inv;
                }
            }

            // update the bottom diagonal of the divided differences array
            const unsigned int n = abscissae.size();
            bottomDiagonal.insert(bottomDiagonal.begin() + n - i, y);
            v_double bottom0 = y;
            for (unsigned int j = i; j < n; ++j) {
                v_double bottom1 = bottomDiagonal[n - (j + 1)];
                const double denom = (x - abscissae[n - (j + 1)]);
                if ( denom == 0.0 ) {
                    throw std::runtime_error( "Dividing by 0 during addition of sample ploint to interpolator." );
                }
                const double inv = 1.0 / denom;
                for (unsigned int k = 0; k < y.size(); ++k) {
                    bottom1[k] = inv * (bottom0[k] - bottom1[k]);
                }
                bottom0 = bottom1;
            }

            // update the top diagonal of the divided differences array
            topDiagonal.push_back(bottom0);

            // update the abscissae array
            abscissae.push_back(x);

        }

    }


    /** Interpolate value at a specified abscissa.
     * <p>
     * Calling this method is equivalent to call the {@link PolynomialFunction#value(double)
     * value} methods of all polynomials returned by {@link #getPolynomials() getPolynomials},
     * except it does not build the intermediate polynomials, so this method is faster and
     * numerically more stable.
     * </p>
     * @param x interpolation abscissa
     * @return interpolated value
     * @exception MathIllegalArgumentException if sample is empty
     */
    v_double value(double x) {

        v_double value(topDiagonal[0].size());
        double valueCoeff = 1;
        for (unsigned int i = 0; i < topDiagonal.size(); ++i) {
            v_double dividedDifference = topDiagonal[i];
            for (unsigned int k = 0; k < value.size(); ++k) {
                value[k] += dividedDifference[k] * valueCoeff;
            }
            const double deltaX = x - abscissae[i];
            valueCoeff *= deltaX;
        }

        return value;

    }


};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_JACOBIPOLYNOMIALS_H
