#ifndef TUDAT_DSST_POLYNOMIALFUNCTION_H
#define TUDAT_DSST_POLYNOMIALFUNCTION_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "derivativeStructure.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

class PolynomialFunction
{
private:
    Vectord coefficients;

public:
    PolynomialFunction(Vectord c) {
        unsigned int n = c.size();
        if (n == 0) {
            throw std::runtime_error( "Empty polynomials coefficients array" );
        }
        while ((n > 1) && (c[n - 1] == 0)) {
            --n;
        }
        coefficients = Vectord(n);
        for ( unsigned int i = 0; i < n; i++ )
        {
            coefficients[i] = c[i];
        }
    }

    static double evaluate(Vectord coefficients, double argument) {
        int n = coefficients.size();
        if (n == 0) {
            throw std::runtime_error( "Empty polynomials coefficients array" );
        }
        double result = coefficients[n - 1];
        for (int j = n - 2; j >= 0; j--) {
            result = argument * result + coefficients[j];
        }
        return result;
    }

    double value(double x) {
        return evaluate(coefficients, x);
    }


    DerivativeStructure value(const DerivativeStructure t) {
        int n = coefficients.size();
        if (n == 0) {
            throw std::runtime_error( "Empty polynomials coefficients array" );
        }
        DerivativeStructure result(t.getFreeParameters(), t.getOrder(), coefficients[n - 1]);
        for (int j = n - 2; j >= 0; j--) {
            result = result.multiply(t).add(coefficients[j]);
        }
        return result;
    }


    /**
     * Add a polynomial to the instance.
     *
     * @param p Polynomial to add.
     * @return a new polynomial which is the sum of the instance and {@code p}.
     */
    PolynomialFunction add(const PolynomialFunction p) {
        // identify the lowest degree polynomial
        const unsigned int lowLength  = std::min(coefficients.size(), p.coefficients.size());
        const unsigned int highLength = std::max(coefficients.size(), p.coefficients.size());

        // build the coefficients array
        Vectord newCoefficients(highLength);
        for (unsigned int i = 0; i < lowLength; ++i) {
            newCoefficients[i] = coefficients[i] + p.coefficients[i];
        }
        Vectord longCoefficients = (coefficients.size() < p.coefficients.size()) ? p.coefficients : coefficients;
        for (unsigned int i = lowLength; i < highLength; ++i) {
            newCoefficients[i] = longCoefficients[i];
        }

        return PolynomialFunction(newCoefficients);
    }


    /**
     * Multiply the instance by a polynomial.
     *
     * @param p Polynomial to multiply by.
     * @return a new polynomial equal to this times {@code p}
     */
    PolynomialFunction multiply(const PolynomialFunction p) {
        Vectord newCoefficients = Vectord(coefficients.size() + p.coefficients.size() - 1);

        for (int i = 0; i < (int) newCoefficients.size(); ++i) {
            newCoefficients[i] = 0.0;
            for (int j = std::max(0, i + 1 - ((int) p.coefficients.size()));
                 j < std::min((int) coefficients.size(), i + 1); ++j ) {
                newCoefficients[i] += coefficients[j] * p.coefficients[i-j];
            }
        }

        return PolynomialFunction(newCoefficients);
    }

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_POLYNOMIALFUNCTION_H
