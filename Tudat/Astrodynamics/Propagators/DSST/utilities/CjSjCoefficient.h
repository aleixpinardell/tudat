#ifndef TUDAT_DSST_CJSJCOEFFICIENT_H
#define TUDAT_DSST_CJSJCOEFFICIENT_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "complex"


namespace tudat
{

namespace propagators
{

namespace dsst
{

typedef std::complex< double > Complex;

class CjSjCoefficient {

private:
    /** Last computed order j. */
    int jLast;

    /** Complex base (k + ih) of the C<sub>j</sub>, S<sub>j</sub> series. */
    Complex kih;

    /** List of computed elements. */
    std::vector< Complex > cjsj;

public:
    /** C<sub>j</sub>(k, h) and S<sub>j</sub>(k, h) constructor.
     * @param k k value
     * @param h h value
     */
    CjSjCoefficient(const double k, const double h) {
        kih  = Complex(k, h);
        cjsj.push_back(Complex(1, 0));
        cjsj.push_back(kih);
        jLast = 1;
    }

    /** Get the C<sub>j</sub> coefficient.
     * @param j order
     * @return C<sub>j</sub>
     */
    double getCj(const int j) {
        if (j > jLast) {
            // Update to order j
            updateCjSj(j);
        }
        return cjsj[j].real();
    }

    /** Get the S<sub>j</sub> coefficient.
     * @param j order
     * @return S<sub>j</sub>
     */
    double getSj(const int j) {
        if (j > jLast) {
            // Update to order j
            updateCjSj(j);
        }
        return cjsj[j].imag();
    }

    /** Get the dC<sub>j</sub> / dk coefficient.
     * @param j order
     * @return dC<sub>j</sub> / d<sub>k</sub>
     */
    double getDcjDk(const int j) {
        return j == 0 ? 0 : j * getCj(j - 1);
    }

    /** Get the dS<sub>j</sub> / dk coefficient.
     * @param j order
     * @return dS<sub>j</sub> / d<sub>k</sub>
     */
    double getDsjDk(const int j) {
        return j == 0 ? 0 : j * getSj(j - 1);
    }

    /** Get the dC<sub>j</sub> / dh coefficient.
     * @param j order
     * @return dC<sub>i</sub> / d<sub>k</sub>
     */
    double getDcjDh(const int j) {
        return j == 0 ? 0 : -j * getSj(j - 1);
    }

    /** Get the dS<sub>j</sub> / dh coefficient.
     * @param j order
     * @return dS<sub>j</sub> / d<sub>h</sub>
     */
    double getDsjDh(const int j) {
        return j == 0 ? 0 : j * getCj(j - 1);
    }

private:
    /** Update the cjsj up to order j.
     * @param j order
     */
    void updateCjSj(const int j) {
        Complex last = *cjsj.end();
        for (int i = jLast; i < j; i++) {
            const Complex next = last * kih;
            cjsj.push_back(next);
            last = next;
        }
        jLast = j;
    }
};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_CJSJCOEFFICIENT_H
