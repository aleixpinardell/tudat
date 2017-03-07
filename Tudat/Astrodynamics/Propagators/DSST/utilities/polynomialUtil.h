#ifndef TUDAT_DSST_POLYNOMIALUTIL_H
#define TUDAT_DSST_POLYNOMIALUTIL_H

#include "boost/boost/rational.hpp"

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "jacobiPolynomials.h"
#include "polynomialFunction.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

typedef boost::rational< int > BigFraction;


/*
class BigFraction : public boost::rational< int >
{
public:
    BigFraction( int n, int d ) : boost::rational< int >( n, d ) { }

    double doubleValue() const
    {
        return double( numerator() ) / double( denominator() );
    }
};
*/

/*
{
private:
    const int numerator;
    const int denominator;

public:
    BigFraction( const int numerator, const int denominator ) : numerator(numerator), denominator(denominator) { }

    double doubleValue() const
    {
        return double( numerator ) / double( denominator );
    }
};
*/

/** Interface for recurrence coefficients generation. */
class RecurrenceCoefficientsGenerator {
public:
    /**
     * Generate recurrence coefficients.
     * @param k highest degree of the polynomials used in the recurrence
     * @return an array of three coefficients such that
     * \( P_{k+1}(x) = (a[0] + a[1] x) P_k(x) - a[2] P_{k-1}(x) \)
     */
    virtual std::vector<BigFraction> generate(int k);
};


class JacobiRecurrenceCoefficientsGenerator : public RecurrenceCoefficientsGenerator
{
private:
    const int v;
    const int w;
public:
    JacobiRecurrenceCoefficientsGenerator( const int v, const int w ) : v(v), w(w) { }

    std::vector<BigFraction> generate(int k)
    {
        k++;
        const int kvw      = k + v + w;
        const int twoKvw   = kvw + k;
        const int twoKvwM1 = twoKvw - 1;
        const int twoKvwM2 = twoKvw - 2;
        const int den      = 2 * k *  kvw * twoKvwM2;

        return std::vector<BigFraction> ( {   BigFraction(twoKvwM1 * (v * v - w * w), den),
                                              BigFraction(twoKvwM1 * twoKvw * twoKvwM2, den),
                                              BigFraction(2 * (k + v - 1) * (k + w - 1) * twoKvw, den)
                                          } );
    }
};


class PolynomialsUtils
{

private:
    static std::map<JacobiKey, std::vector<BigFraction>> JACOBI_COEFFICIENTS;

public:


    /** Compute polynomial coefficients up to a given degree.
      * @param degree maximal degree
      * @param maxDegree current maximal degree
      * @param generator recurrence coefficients generator
      * @param coefficients list where the computed coefficients should be appended
      */
     static void computeUpToDegree( const int degree, const int maxDegree,
             RecurrenceCoefficientsGenerator generator, std::vector<BigFraction> coefficients) {

         int startK = (maxDegree - 1) * maxDegree / 2;
         for (int k = maxDegree; k < degree; ++k) {

             // start indices of two previous polynomials Pk(X) and Pk-1(X)
             int startKm1 = startK;
             startK += k;

             // Pk+1(X) = (a[0] + a[1] X) Pk(X) - a[2] Pk-1(X)
             std::vector<BigFraction> ai = generator.generate(k);

             BigFraction ck     = coefficients[startK];
             BigFraction ckm1   = coefficients[startKm1];

             // degree 0 coefficient
             coefficients.push_back( ck * ai[0] - ckm1 * ai[2] );

             // degree 1 to degree k-1 coefficients
             for (int i = 1; i < k; ++i) {
                 const BigFraction ckPrev = ck;
                 ck     = coefficients[ startK + i ];
                 ckm1   = coefficients[ startKm1 + i ];
                 coefficients.push_back( ck * ai[0] + ckPrev * ai[1] - ckm1 * ai[2] );
             }

             // degree k coefficient
             const BigFraction ckPrev = ck;
             ck = coefficients[ startK + k ];
             coefficients.push_back( ck * ai[0] + ckPrev * ai[1] );

             // degree k+1 coefficient
             coefficients.push_back( ck * ai[1] );

         }

     }


    /** Get the coefficients array for a given degree.
      * @param degree degree of the polynomial
      * @param coefficients list where the computed coefficients are stored
      * @param generator recurrence coefficients generator
      * @return coefficients array
      */
     static PolynomialFunction buildPolynomial( const int degree, const std::vector<BigFraction> coefficients,
                                                const RecurrenceCoefficientsGenerator generator ) {

             const int maxDegree = (int) std::floor(std::sqrt(2 * coefficients.size())) - 1;
             if (degree > maxDegree) {
                 computeUpToDegree(degree, maxDegree, generator, coefficients);
             }

         // coefficient  for polynomial 0 is  l [0]
         // coefficients for polynomial 1 are l [1] ... l [2] (degrees 0 ... 1)
         // coefficients for polynomial 2 are l [3] ... l [5] (degrees 0 ... 2)
         // coefficients for polynomial 3 are l [6] ... l [9] (degrees 0 ... 3)
         // coefficients for polynomial 4 are l[10] ... l[14] (degrees 0 ... 4)
         // coefficients for polynomial 5 are l[15] ... l[20] (degrees 0 ... 5)
         // coefficients for polynomial 6 are l[21] ... l[27] (degrees 0 ... 6)
         // ...
         const int start = degree * (degree + 1) / 2;

         SingleVectord a(degree + 1);
         for (int i = 0; i <= degree; ++i) {
             a[i] = boost::rational_cast<double>( coefficients[start + i] );
         }

         // build the polynomial
         return PolynomialFunction(a);

     }


    /**
     * Create a Jacobi polynomial.
     * <p><a href="http://mathworld.wolfram.com/JacobiPolynomial.html">Jacobi
     * polynomials</a> are orthogonal polynomials.
     * They can be defined by the following recurrence relations:</p><p>
     * \(
     *    P_0^{vw}(x) = 1 \\
     *    P_{-1}^{vw}(x) = 0 \\
     *    2k(k + v + w)(2k + v + w - 2) P_k^{vw}(x) = \\
     *    (2k + v + w - 1)[(2k + v + w)(2k + v + w - 2) x + v^2 - w^2] P_{k-1}^{vw}(x) \\
     *  - 2(k + v - 1)(k + w - 1)(2k + v + w) P_{k-2}^{vw}(x)
     * \)
     * </p>
     * @param degree degree of the polynomial
     * @param v first exponent
     * @param w second exponent
     * @return Jacobi polynomial of specified degree
     */
    static PolynomialFunction createJacobiPolynomial(const int degree, const int v, const int w)
    {
        // select the appropriate list
        const JacobiKey key(v, w);

        if (JACOBI_COEFFICIENTS.count(key) == 0) {

            // allocate a new list for v, w
            std::vector<BigFraction> list;
            JACOBI_COEFFICIENTS[key] = list;

            // Pv,w,0(x) = 1;
            list.push_back(BigFraction(1, 1));

            // P1(x) = (v - w) / 2 + (2 + v + w) * X / 2
            list.push_back(BigFraction(v - w, 2));
            list.push_back(BigFraction(2 + v + w, 2));

        }

        return buildPolynomial( degree, JACOBI_COEFFICIENTS[key], JacobiRecurrenceCoefficientsGenerator(v,w) );

    }


};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_POLYNOMIALUTIL_H
