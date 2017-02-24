#ifndef TUDAT_DSST_COEFFICIENTSFACTORY_H
#define TUDAT_DSST_COEFFICIENTSFACTORY_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{


//! Key (N,S) formed by a pair of integers.
class NSKey : public std::pair< int, int > {
public:
    NSKey( int n, int s ) : std::pair< int, int >( n, s ) { }

    int getN() {
        return first;
    }

    int getS() {
        return second;
    }
};

typedef std::map< NSKey, double > NSMap;
typedef std::pair< vv_double, vv_double > NSSetPair;


class CoefficientsFactory
{
public:
    //! Initializer.
    CoefficientsFactory();

    //! Function to compute the dK_0^{n,s}/d\chi coefficients for a given \chi from the recurrence formula 3.2-(3)
    /*!
     * \brief Function to compute the derivatives of the K_0^{n,s} coefficients with respect to \chi.
     * \param chi 1/(sqrt(1-e^2)).
     * \param nMax Maximum value for n.
     * \param sMax Maximum value for s.
     * \return A pair containing the K_0^{n,s} coefficients and the derivatives.
     */
    NSSetPair computeK0nsAndDerivatives( const double chi, const int nMax, const int sMax );

    //! Function to compute the K_0^{n,s} coefficients for a given \chi from the recurrence formula 2.3.7-(7)
    vv_double computeK0ns( const double chi, const int nMax, const int sMax );

    //! Function to compute the Q_{n,s} coefficients evaluated at \gamma from the recurrence formula 2.8.3-(2)
    vv_double computeQns( const double gamma, const int nMax, const int sMax );

    //! Function to compute the dQ_{n,s}/d\gamma coefficients evaluated at \gamma from the recurrence formula 3.1-(8)
    vv_double computeQnsDerivatives( const double gamma, const int nMax, const int sMax );

    //! Function to compute recursively G_s and H_s polynomials from equation 3.1-(5)
    Eigen::Matrix2Xd computeGsHs(
            const double k, const double h, const double alpha, const double beta, const int order );

    //! Function to compute derivatives of G_s with respect to k, h, alpha and beta from equation 3.1-(9)
    std::map< std::string, Eigen::VectorXd > computeGsDerivatives(
            const double k, const double h, const double alpha, const double beta, const int order );

    //! Function to compute the V_{n,s} coefficients from Eq. 2.8.2-(1/2)
    /*!
     * Returns a map with the coefficients for n = 0...order-1, s = 0...n
     * \param order Order of the computation.
     * \return map  Map with the V_{n,s} coefficients.
     */
    NSMap computeVns( const int order );

    //! Function to get a V_{n,s}^m coefficient from Eq. 2.7.1-(6) and Sec. 2.8.2
    /*!
     * Get a V_{n,s}^m coefficient.
     * \param m m
     * \param n n
     * \param s s
     * \return The V_{n,s}^m coefficient
     * \throws std::runtime_error when m > n
     */
    double getVmns( const int m, const int n, const int s );

private:
    //! Internal storage of the polynomial values. Reused for further computation.
    NSMap VNS;

    //! Last computed order for V<sub>ns</sub> coefficients.
    int LAST_VNS_ORDER = 2;

};





} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_COEFFICIENTSFACTORY_H
