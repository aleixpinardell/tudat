#ifndef TUDAT_DSST_JACOBIPOLYNOMIALS_H
#define TUDAT_DSST_JACOBIPOLYNOMIALS_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "polynomialFunction.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

/** Inner class for Jacobi polynomials keys.
 * <p>
 * Please note that this class is not original content but is a copy from the
 * Hipparchus library. This library is published under the
 * Apache License, version 2.0.
 * </p>
 *
 * @see org.hipparchus.analysis.polynomials.PolynomialsUtils
 */
class JacobiKey {

private:
    /** First exponent. */
    const int v;

    /** Second exponent. */
    const int w;

public:
    /** Simple constructor.
     * @param v first exponent
     * @param w second exponent
     */
    JacobiKey(const int v, const int w) : v(v), w(w) { }

    /** Get hash code.
     * @return hash code
     */
    int hashCode() {
        return (v << 16) ^ w;
    }

    /** Check if the instance represent the same key as another instance.
     * @param key other key
     * @return true if the instance and the other key refer to the same polynomial
     */
    bool equals(const JacobiKey key) const {
        return (v == key.v) && (w == key.w);
    }


    bool operator < (const JacobiKey& key) const
    {
        return v == key.v ? w < key.w : v < key.v;
    }

};


/** Provider of the Jacobi polynomials P<sub>l</sub><sup>v,w</sup>.
 * <p>
 * This class is used for DSSTTesseral contribution computation.
 * </p>
 *
 */
class JacobiPolynomials {

private:

    /** Storage map. */
    static std::map<JacobiKey, std::vector<PolynomialFunction>> MAP;

    /** Private constructor as class is a utility. */
    JacobiPolynomials() { }


public:
    /** Returns the value and derivatives of the Jacobi polynomial P<sub>l</sub><sup>v,w</sup> evaluated at γ.
     * <p>
     * This method is guaranteed to be thread-safe
     * </p>
     * @param l degree of the polynomial
     * @param v v value
     * @param w w value
     * @param gamma γ value
     * @return value and derivatives of the Jacobi polynomial P<sub>l</sub><sup>v,w</sup>(γ)
     */
    static DerivativeStructure getValue(const int l, const int v, const int w, const DerivativeStructure gamma);

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_JACOBIPOLYNOMIALS_H
