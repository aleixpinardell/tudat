#ifndef TUDAT_DSST_THIRDBODY_H
#define TUDAT_DSST_THIRDBODY_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/forceModel.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/CjSjCoefficient.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/shortPeriodicsInterpolatedCoefficient.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/shortPeriodTerms.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/spacecraftState.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/timeSpanMap.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/upperBounds.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{


class DSSTThirdBody : public DSSTForceModel
{
private:
    Body & thirdBody;

public:
    //! Constructor.
    /** Simple constructor.
     * \param gravitationalParameter Gravitational parameter of the third body [m^3/s^2].
     * \param distance Distance to the thrid body [m]
     */
    DSSTThirdBody( const Body &thirdBody, AuxiliaryElements &auxiliaryElements );

    std::vector< ShortPeriodTerms > initialize( AuxiliaryElements auxiliaryElements , const bool meanOnly );

    void update( AuxiliaryElements &auxiliaryElements );

    Eigen::Vector6d getMeanElementRates( );

    Eigen::Vector6d computeUDerivatives( );

    //! Get the position of the third body with respect to the central body [ m, m, m ]
    Eigen::Vector3d getPositionFromCentralBody() const {
        return thirdBody.getPosition() - auxiliaryElements.centralBody.getPosition();
    }

    //! Get the distance between the third body and the central body [ m ]
    double getDistanceFromCentralBody() const {
        return getPositionFromCentralBody().norm();
    }


private:
    //! Max power for summation
    static const int    MAX_POWER = 22;

    //! Truncation tolerance for big, eccentric  orbits
    static constexpr double BIG_TRUNCATION_TOLERANCE = 1.e-1;

    //! Truncation tolerance for small orbits
    static constexpr double SMALL_TRUNCATION_TOLERANCE = 1.9e-6;

    //! Number of points for interpolation
    static const int    INTERPOLATION_POINTS = 3;

    //! Maximum power for eccentricity used in short periodic computation
    static const int    MAX_ECCPOWER_SP = 4;

    //! Retrograde factor I
    static const int    I = 1;

    /* FIXME
    //! The 3rd body to consider
    const CelestialBody    body;
    */

    //! Name of the third body
    std::string bodyName;

    //! Position of the third body
    Eigen::Vector3d bodyPosition;

    //! Standard gravitational parameter μ for the body in m³/s²
    double gm;

    //! Distance from center of mass of the central body to the 3rd body
    double R3;

    // Equinoctial elements (according to DSST notation)
    //! a
    double a;
    //! e_y
    double h;
    //! e_x
    double k;
    //! h_y
    double p;
    //! h_x
    double q;

    //! Eccentricity
    double ecc;

    // Direction cosines of the symmetry axis
    //! α
    double alpha;
    //! β
    double beta;
    //! γ
    double gamma;

    // Common factors for potential computation
    //! A = n * a²
    double A;
    //! B = sqrt(1 - e²)
    double B;
    //! C = 1 + p² + q²
    double C;
    //! B²
    double B2;
    //! B³
    double B3;

    //! The mean motion (n)
    double meanMotion;

    //! \Chi = 1 / sqrt(1 - e²) = 1 / B
    double Chi;
    //! \Chi²
    double Chi2;
    //! \Chi³
    double Chi3;
    //! -2 * a / A
    double m2aoA;
    //! B / A
    double BoA;
    //! 1 / (A * B)
    double ooAB;
    //! -C / (2 * A * B)
    double mCo2AB;
    //! B / A(1 + B)
    double BoABpo;

    //! Max power for a/R3 in the serie expansion
    unsigned int maxAR3Pow = INT_MIN;

    //! Max power for e in the serie expansion
    unsigned int maxEccPow = INT_MIN;

    //! Max power for e in the serie expansion (for short periodics)
    unsigned int maxEccPowShort;

    //! Max frequency of F
    unsigned int maxFreqF;

    /* FIXME
    //! An array that contains the objects needed to build the Hansen coefficients. The index is s
    const HansenThirdBodyLinear[] hansenObjects;
    */

    //! The current value of the mean disturbing function. Needed for the short periodic contribution
    double U;

    //! Vns coefficients
    coefficients_factories::VnsCoefficientsFactory VnsFactory;

    /*
    //! K0ns coefficients
    coefficients_factories::K0nsCoefficientsFactory K0nsFactory;
    */

    //! Qns coefficients
    coefficients_factories::QnsCoefficientsFactory QnsFactory;
    // DoubleVectord Qns;

    //! Gs coefficients
    coefficients_factories::GsHsCoefficientsFactory GsFactory;

    //! a / R3 up to power maxAR3Pow
    SingleVectord aoR3Pow;

    //! mu3 / R3
    double muoR3;

    //! b = 1 / (1 + sqrt(1 - e²)) = 1 / (1 + B)
    double b;

    //! h * \Chi³
    double hChi3;

    //! k * \Chi³
    double kChi3;




    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// GnsCoefficients ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////


    /** The G<sub>n,s</sub> coefficients and their derivatives.
      * <p>
      * See Danielson, 4.2-17
      *
      */
    class GnsCoefficients {

    private:
        DSSTThirdBody * const owner;

        /** Maximum value for n index. */
        int nMax;

        /** Maximum value for s index. */
        int sMax;

        /** The coefficients G<sub>n,s</sub>. */
        DoubleVectord gns;

        /** The derivatives of the coefficients G<sub>n,s</sub> by a. */
        DoubleVectord dgnsda;

        /** The derivatives of the coefficients G<sub>n,s</sub> by γ. */
        DoubleVectord dgnsdgamma;

        /**
         * Compute the coefficient G<sub>n,s</sub> and its derivatives.
         * <p>
         * Only the derivatives by a and γ are computed as all others are 0
         * </p>
         */
        void generateCoefficients();

    public:
        /** Standard constructor.
         *
         * @param nMax maximum value for n indes
         * @param sMax maximum value for s index
         */
        GnsCoefficients( DSSTThirdBody &owner, const int nMax, const int sMax);

        /** Get the coefficient G<sub>n,s</sub>.
         *
         * @param n n index
         * @param s s index
         * @return the coefficient G<sub>n,s</sub>
         */
        double getGns(const int n, const int s) {
            return this->gns[n][s];
        }

        /** Get the derivative dG<sub>n,s</sub> / da.
         *
         * @param n n index
         * @param s s index
         * @return the derivative dG<sub>n,s</sub> / da
         */
        double getdGnsda(const int n, const int s) {
            return this->dgnsda[n][s];
        }

        /** Get the derivative dG<sub>n,s</sub> / dγ.
         *
         * @param n n index
         * @param s s index
         * @return the derivative dG<sub>n,s</sub> / dγ
         */
        double getdGnsdgamma(const int n, const int s) {
            return this->dgnsdgamma[n][s];
        }
    };




    //////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// WnsjEtomjmsCoefficient ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    /** This class covers the coefficients e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup>
     * and their derivatives by h and k.
     *
     * <p>
     * Starting from Danielson 4.2-9,10,11 and taking into account that fact that: <br />
     * c = e / (1 + (1 - e²)<sup>1/2</sup>) = e / (1 + B) = e * b <br/>
     * the expression e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup>
     * can be written as: <br/ >
     * - for |s| > |j| <br />
     * e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup> =
     *          (((n + s)!(n - s)!)/((n + j)!(n - j)!)) *
     *          (-b)<sup>|j-s|</sup> *
     *          ((1 - c²)<sup>n-|s|</sup>/(1 + c²)<sup>n</sup>) *
     *          P<sub>n-|s|</sub><sup>|j-s|, |j+s|</sup>(χ) <br />
     * <br />
     * - for |s| <= |j| <br />
     * e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup> =
     *          (-b)<sup>|j-s|</sup> *
     *          ((1 - c²)<sup>n-|j|</sup>/(1 + c²)<sup>n</sup>) *
     *          P<sub>n-|j|</sub><sup>|j-s|, |j+s|</sup>(χ)
     * </p>
     *
     */
    class WnsjEtomjmsCoefficient {

    private:
        DSSTThirdBody * const owner;

        /** The value c.
         * <p>
         *  c = e / (1 + (1 - e²)<sup>1/2</sup>) = e / (1 + B) = e * b <br/>
         * </p>
         *  */
        double c;

        /** c².*/
        double c2;

        /** db / dh. */
        double dbdh;

        /** db / dk. */
        double dbdk;

        /** de / dh. */
        double dedh;

        /** de / dk. */
        double dedk;

        /** dc / dh = e * db/dh + b * de/dh. */
        double dcdh;

        /** dc / dk = e * db/dk + b * de/dk. */
        double dcdk;

        /** The values (1 - c²)<sup>n</sup>. <br />
         * The maximum possible value for the power is N + 1 */
        SingleVectord omc2tn;

        /** The values (1 + c²)<sup>n</sup>. <br />
         * The maximum possible value for the power is N + 1 */
        SingleVectord opc2tn;

        /** The values b<sup>|j-s|</sup>. */
        SingleVectord btjms;

    public:
        /**
         * Standard constructor.
         */
        WnsjEtomjmsCoefficient( DSSTThirdBody & owner );

        /** Compute the value of the coefficient e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup>
         * and its derivatives by h and k. <br />
         *
         * @param j j index
         * @param s s index
         * @param n n index
         * @return an array containing the value of the coefficient at index 0, the derivative by k
         * at index 1 and the derivative by h at index 2
         */
        SingleVectord computeWjnsEmjmsAndDeriv(const int j, const int s, const int n);

    };





    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// CjSjAlphaBetaKH ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    /** This class computes the terms containing the coefficients C<sub>j</sub> and S<sub>j</sub> of (α, β) or (k, h).
     *
     * <p>
     * The following terms and their derivatives by k, h, alpha and beta are considered: <br/ >
     * - sign(j-s) * C<sub>s</sub>(α, β) * S<sub>|j-s|</sub>(k, h)
     * + S<sub>s</sub>(α, β) * C<sub>|j-s|</sub>(k, h) <br />
     * - C<sub>s</sub>(α, β) * S<sub>j+s</sub>(k, h) - S<sub>s</sub>(α, β) * C<sub>j+s</sub>(k, h) <br />
     * - C<sub>s</sub>(α, β) * C<sub>|j-s|</sub>(k, h)
     * - sign(j-s) * S<sub>s</sub>(α, β) * S<sub>|j-s|</sub>(k, h) <br />
     * - C<sub>s</sub>(α, β) * C<sub>j+s</sub>(k, h) + S<sub>s</sub>(α, β) * S<sub>j+s</sub>(k, h) <br />
     * For the ease of usage the above terms are renamed A<sub>js</sub>, B<sub>js</sub>, D<sub>js</sub>
     * and E<sub>js</sub> respectively <br />
     * See the CS Mathematical report $3.5.3.2 for more details
     * </p>
     */
    class CjSjAlphaBetaKH {

    private:
        DSSTThirdBody * const owner;

        /** The C<sub>j</sub>(k, h) and the S<sub>j</sub>(k, h) series. */
        CjSjCoefficient cjsjkh;

        /** The C<sub>j</sub>(α, β) and the S<sub>j</sub>(α, β) series. */
        CjSjCoefficient cjsjalbe;

        /** The coeficient sign(j-s) * C<sub>s</sub>(α, β) * S<sub>|j-s|</sub>(k, h)
     * + S<sub>s</sub>(α, β) * C<sub>|j-s|</sub>(k, h)
     * and its derivative by k, h, α and β. */
        SingleVectord coefAandDeriv = SingleVectord(5);

        /** The coeficient C<sub>s</sub>(α, β) * S<sub>j+s</sub>(k, h)
     * - S<sub>s</sub>(α, β) * C<sub>j+s</sub>(k, h)
     * and its derivative by k, h, α and β. */
        SingleVectord coefBandDeriv = SingleVectord(5);

        /** The coeficient C<sub>s</sub>(α, β) * C<sub>|j-s|</sub>(k, h)
     * - sign(j-s) * S<sub>s</sub>(α, β) * S<sub>|j-s|</sub>(k, h)
     * and its derivative by k, h, α and β. */
        SingleVectord coefDandDeriv = SingleVectord(5);

        /** The coeficient C<sub>s</sub>(α, β) * C<sub>j+s</sub>(k, h) + S<sub>s</sub>(α, β) * S<sub>j+s</sub>(k, h)
     * and its derivative by k, h, α and β. */
        SingleVectord coefEandDeriv = SingleVectord(5);

    public:
        /**
     * Standard constructor.
     */
        CjSjAlphaBetaKH( DSSTThirdBody &owner) : owner(&owner), cjsjkh(CjSjCoefficient(owner.k, owner.h)),
            cjsjalbe(CjSjCoefficient(owner.alpha, owner.beta)) {}

        /** Compute the coefficients and their derivatives for a given (j,s) pair.
     * @param j j index
     * @param s s index
     */
        void computeCoefficients(const int j, const int s);

        /** Get the value of coefficient A<sub>j,s</sub>.
     *
     * @return the coefficient A<sub>j,s</sub>
     */
        double getCoefA() {
            return coefAandDeriv[0];
        }

        /** Get the value of coefficient dA<sub>j,s</sub>/dk.
     *
     * @return the coefficient dA<sub>j,s</sub>/dk
     */
        double getdCoefAdk() {
            return coefAandDeriv[1];
        }

        /** Get the value of coefficient dA<sub>j,s</sub>/dh.
     *
     * @return the coefficient dA<sub>j,s</sub>/dh
     */
        double getdCoefAdh() {
            return coefAandDeriv[2];
        }

        /** Get the value of coefficient dA<sub>j,s</sub>/dα.
     *
     * @return the coefficient dA<sub>j,s</sub>/dα
     */
        double getdCoefAdalpha() {
            return coefAandDeriv[3];
        }

        /** Get the value of coefficient dA<sub>j,s</sub>/dβ.
     *
     * @return the coefficient dA<sub>j,s</sub>/dβ
     */
        double getdCoefAdbeta() {
            return coefAandDeriv[4];
        }

        /** Get the value of coefficient B<sub>j,s</sub>.
     *
     * @return the coefficient B<sub>j,s</sub>
     */
        double getCoefB() {
            return coefBandDeriv[0];
        }

        /** Get the value of coefficient dB<sub>j,s</sub>/dk.
     *
     * @return the coefficient dB<sub>j,s</sub>/dk
     */
        double getdCoefBdk() {
            return coefBandDeriv[1];
        }

        /** Get the value of coefficient dB<sub>j,s</sub>/dh.
     *
     * @return the coefficient dB<sub>j,s</sub>/dh
     */
        double getdCoefBdh() {
            return coefBandDeriv[2];
        }

        /** Get the value of coefficient dB<sub>j,s</sub>/dα.
     *
     * @return the coefficient dB<sub>j,s</sub>/dα
     */
        double getdCoefBdalpha() {
            return coefBandDeriv[3];
        }

        /** Get the value of coefficient dB<sub>j,s</sub>/dβ.
     *
     * @return the coefficient dB<sub>j,s</sub>/dβ
     */
        double getdCoefBdbeta() {
            return coefBandDeriv[4];
        }

        /** Get the value of coefficient D<sub>j,s</sub>.
     *
     * @return the coefficient D<sub>j,s</sub>
     */
        double getCoefD() {
            return coefDandDeriv[0];
        }

        /** Get the value of coefficient dD<sub>j,s</sub>/dk.
     *
     * @return the coefficient dD<sub>j,s</sub>/dk
     */
        double getdCoefDdk() {
            return coefDandDeriv[1];
        }

        /** Get the value of coefficient dD<sub>j,s</sub>/dh.
     *
     * @return the coefficient dD<sub>j,s</sub>/dh
     */
        double getdCoefDdh() {
            return coefDandDeriv[2];
        }

        /** Get the value of coefficient dD<sub>j,s</sub>/dα.
     *
     * @return the coefficient dD<sub>j,s</sub>/dα
     */
        double getdCoefDdalpha() {
            return coefDandDeriv[3];
        }

        /** Get the value of coefficient dD<sub>j,s</sub>/dβ.
     *
     * @return the coefficient dD<sub>j,s</sub>/dβ
     */
        double getdCoefDdbeta() {
            return coefDandDeriv[4];
        }

        /** Get the value of coefficient E<sub>j,s</sub>.
     *
     * @return the coefficient E<sub>j,s</sub>
     */
        double getCoefE() {
            return coefEandDeriv[0];
        }

        /** Get the value of coefficient dE<sub>j,s</sub>/dk.
     *
     * @return the coefficient dE<sub>j,s</sub>/dk
     */
        double getdCoefEdk() {
            return coefEandDeriv[1];
        }

        /** Get the value of coefficient dE<sub>j,s</sub>/dh.
     *
     * @return the coefficient dE<sub>j,s</sub>/dh
     */
        double getdCoefEdh() {
            return coefEandDeriv[2];
        }

        /** Get the value of coefficient dE<sub>j,s</sub>/dα.
     *
     * @return the coefficient dE<sub>j,s</sub>/dα
     */
        double getdCoefEdalpha() {
            return coefEandDeriv[3];
        }

        /** Get the value of coefficient dE<sub>j,s</sub>/dβ.
     *
     * @return the coefficient dE<sub>j,s</sub>/dβ
     */
        double getdCoefEdbeta() {
            return coefEandDeriv[4];
        }
    };




    ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// FourierCjSjCoefficients ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    /** Computes the C<sup>j</sup> and S<sup>j</sup> coefficients Danielson 4.2-(15,16)
     * and their derivatives.
     *  <p>
     *  CS Mathematical Report $3.5.3.2
     *  </p>
     */
    class FourierCjSjCoefficients {

    private:
        DSSTThirdBody * const owner;

        /** The coefficients G<sub>n, s</sub> and their derivatives. */
        GnsCoefficients gns;

        /** the coefficients e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup> and their derivatives by h and k. */
        WnsjEtomjmsCoefficient wnsjEtomjmsCoefficient;

        /** The terms containing the coefficients C<sub>j</sub> and S<sub>j</sub> of (α, β) or (k, h). */
        CjSjAlphaBetaKH ABDECoefficients;

        /** The Fourier coefficients C<sup>j</sup> and their derivatives.
         * <p>
         * Each column of the matrix contains the following values: <br/>
         * - C<sup>j</sup> <br/>
         * - dC<sup>j</sup> / da <br/>
         * - dC<sup>j</sup> / dk <br/>
         * - dC<sup>j</sup> / dh <br/>
         * - dC<sup>j</sup> / dα <br/>
         * - dC<sup>j</sup> / dβ <br/>
         * - dC<sup>j</sup> / dγ <br/>
         * </p>
         */
        DoubleVectord cj;

        /** The S<sup>j</sup> coefficients and their derivatives.
         * <p>
         * Each column of the matrix contains the following values: <br/>
         * - S<sup>j</sup> <br/>
         * - dS<sup>j</sup> / da <br/>
         * - dS<sup>j</sup> / dk <br/>
         * - dS<sup>j</sup> / dh <br/>
         * - dS<sup>j</sup> / dα <br/>
         * - dS<sup>j</sup> / dβ <br/>
         * - dS<sup>j</sup> / dγ <br/>
         * </p>
         */
        DoubleVectord sj;

        /** The Coefficients C<sup>j</sup><sub>,λ</sub>.
         * <p>
         * See Danielson 4.2-21
         * </p>
         */
        SingleVectord cjlambda;

        /** The Coefficients S<sup>j</sup><sub>,λ</sub>.
         * <p>
         * See Danielson 4.2-21
         * </p>
         */
        SingleVectord sjlambda;

        /** Maximum value for n. */
        int nMax;

        /** Maximum value for s. */
        int sMax;

        /** Maximum value for j. */
        int jMax;

        /**
         * Compute all coefficients.
         */
        void computeCoefficients();

    public:
        /**
         * Constructor.
         *
         * \param nMax maximum value for n index
         * \param sMax maximum value for s index
         * \param jMax maximum value for j index
         */
        FourierCjSjCoefficients( DSSTThirdBody & owner, const int nMax, const int sMax, const int jMax );

        /** Get the Fourier coefficient C<sup>j</sup>.
         * \param j j index
         * \return C<sup>j</sup>
         */
        double getCj(const int j) {
            return cj[0][j];
        }

        /** Get the derivative dC<sup>j</sup>/da.
         * \param j j index
         * \return dC<sup>j</sup>/da
         */
        double getdCjda(const int j) {
            return cj[1][j];
        }

        /** Get the derivative dC<sup>j</sup>/dk.
         * \param j j index
         * \return dC<sup>j</sup>/dk
         */
        double getdCjdk(const int j) {
            return cj[2][j];
        }

        /** Get the derivative dC<sup>j</sup>/dh.
         * \param j j index
         * \return dC<sup>j</sup>/dh
         */
        double getdCjdh(const int j) {
            return cj[3][j];
        }

        /** Get the derivative dC<sup>j</sup>/dα.
         * \param j j index
         * \return dC<sup>j</sup>/dα
         */
        double getdCjdalpha(const int j) {
            return cj[4][j];
        }

        /** Get the derivative dC<sup>j</sup>/dβ.
         * \param j j index
         * \return dC<sup>j</sup>/dβ
         */
        double getdCjdbeta(const int j) {
            return cj[5][j];
        }

        /** Get the derivative dC<sup>j</sup>/dγ.
         * \param j j index
         * \return dC<sup>j</sup>/dγ
         */
        double getdCjdgamma(const int j) {
            return cj[6][j];
        }

        /** Get the Fourier coefficient S<sup>j</sup>.
         * \param j j index
         * \return S<sup>j</sup>
         */
        double getSj(const int j) {
            return sj[0][j];
        }

        /** Get the derivative dS<sup>j</sup>/da.
         * \param j j index
         * \return dS<sup>j</sup>/da
         */
        double getdSjda(const int j) {
            return sj[1][j];
        }

        /** Get the derivative dS<sup>j</sup>/dk.
         * \param j j index
         * \return dS<sup>j</sup>/dk
         */
        double getdSjdk(const int j) {
            return sj[2][j];
        }

        /** Get the derivative dS<sup>j</sup>/dh.
         * \param j j index
         * \return dS<sup>j</sup>/dh
         */
        double getdSjdh(const int j) {
            return sj[3][j];
        }

        /** Get the derivative dS<sup>j</sup>/dα.
         * \param j j index
         * \return dS<sup>j</sup>/dα
         */
        double getdSjdalpha(const int j) {
            return sj[4][j];
        }

        /** Get the derivative dS<sup>j</sup>/dβ.
         * \param j j index
         * \return dS<sup>j</sup>/dβ
         */
        double getdSjdbeta(const int j) {
            return sj[5][j];
        }

        /** Get the derivative dS<sup>j</sup>/dγ.
         * \param j j index
         * \return dS<sup>j</sup>/dγ
         */
        double getdSjdgamma(const int j) {
            return sj[6][j];
        }

        /** Get the coefficient C⁰<sub>,λ</sub>.
         * \return C⁰<sub>,λ</sub>
         */
        double getC0Lambda() {
            return cjlambda[0];
        }

        /** Get the coefficient C<sup>j</sup><sub>,λ</sub>.
         * \param j j index
         * \return C<sup>j</sup><sub>,λ</sub>
         */
        double getCjLambda(const int j) {
            if (j < 1 || j >= jMax) {
                return 0.;
            }
            return cjlambda[j];
        }

        /** Get the coefficient S<sup>j</sup><sub>,λ</sub>.
         * \param j j index
         * \return S<sup>j</sup><sub>,λ</sub>
         */
        double getSjLambda(const int j) {
            if (j < 1 || j >= jMax) {
                return 0.;
            }
            return sjlambda[j];
        }
    };




    ////////////////////////////////////////////////////////////////////
    /////////////////////////////// Slot ///////////////////////////////
    ////////////////////////////////////////////////////////////////////

    /** Coefficients valid for one time slot. */
    class Slot : public Equatable {

    public:
        /** The coefficients C<sub>i</sub><sup>j</sup>.
         * <p>
         * The index order is cij[j][i] <br/>
         * i corresponds to the equinoctial element, as follows: <br/>
         * - i=0 for a <br/>
         * - i=1 for k <br/>
         * - i=2 for h <br/>
         * - i=3 for q <br/>
         * - i=4 for p <br/>
         * - i=5 for λ <br/>
         * </p>
         */
        std::vector< ShortPeriodicsInterpolatedCoefficient > cij;

        /** The coefficients S<sub>i</sub><sup>j</sup>.
         * <p>
         * The index order is sij[j][i] <br/>
         * i corresponds to the equinoctial element, as follows: <br/>
         * - i=0 for a <br/>
         * - i=1 for k <br/>
         * - i=2 for h <br/>
         * - i=3 for q <br/>
         * - i=4 for p <br/>
         * - i=5 for λ <br/>
         * </p>
         */
        std::vector< ShortPeriodicsInterpolatedCoefficient > sij;

    public:
        Slot() {}

        /** Simple constructor.
         *  @param jMax maximum value for j index
         *  @param interpolationPoints number of points used in the interpolation process
         */
        Slot(const int jMax, const int interpolationPoints);

        bool operator == ( const Slot &other ) const
        {
            return cij == other.cij && sij == other.sij;
        }
    };



    //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// GeneratingFunctionCoefficients ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////


    /** This class computes the coefficients for the generating function S and its derivatives.
     * <p>
     * The form of the generating functions is: <br>
     *  S = C⁰ + &Sigma;<sub>j=1</sub><sup>N+1</sup>(C<sup>j</sup> * cos(jF) + S<sup>j</sup> * sin(jF)) <br>
     *  The coefficients C⁰, C<sup>j</sup>, S<sup>j</sup> are the Fourrier coefficients
     *  presented in Danielson 4.2-14,15 except for the case j=1 where
     *  C¹ = C¹<sub>Fourier</sub> - hU and
     *  S¹ = S¹<sub>Fourier</sub> + kU <br>
     *  Also the coefficients of the derivatives of S by a, k, h, α, β, γ and λ
     *  are computed end expressed in a similar manner. The formulas used are 4.2-19, 20, 23, 24
     * </p>
     */
    class GeneratingFunctionCoefficients {

    private:
        DSSTThirdBody * const owner;

        /** Maximum value of j index. */
        int jMax;

        /** The Fourier coefficients as presented in Danielson 4.2-14,15. */
        FourierCjSjCoefficients cjsjFourier;

        /** The coefficients C<sup>j</sup> of the function S and its derivatives.
         * <p>
         * The index j belongs to the interval [0,jMax]. The coefficient C⁰ is the free coefficient.<br>
         * Each column of the matrix contains the coefficient corresponding to the following functions: <br/>
         * - S <br/>
         * - dS / da <br/>
         * - dS / dk <br/>
         * - dS / dh <br/>
         * - dS / dα <br/>
         * - dS / dβ <br/>
         * - dS / dγ <br/>
         * - dS / dλ
         * </p>
         */
        DoubleVectord cjCoefs;

        /** The coefficients S<sup>j</sup> of the function S and its derivatives.
         * <p>
         * The index j belongs to the interval [0,jMax].<br>
         * Each column of the matrix contains the coefficient corresponding to the following functions: <br/>
         * - S <br/>
         * - dS / da <br/>
         * - dS / dk <br/>
         * - dS / dh <br/>
         * - dS / dα <br/>
         * - dS / dβ <br/>
         * - dS / dγ <br/>
         * - dS / dλ
         * </p>
         */
        DoubleVectord sjCoefs;

        /**
         * Compute the coefficients for the generating function S and its derivatives.
         */
        void computeGeneratingFunctionCoefficients();

    public:
        /**
         * Standard constructor.
         *
         * @param nMax maximum value of n index
         * @param sMax maximum value of s index
         * @param jMax maximum value of j index
         */
        GeneratingFunctionCoefficients(DSSTThirdBody &owner, const int nMax, const int sMax, const int jMax);

        /** Get the coefficient C<sup>j</sup> for the function S.
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function S
         */
        double getSCj(const int j) const {
            return cjCoefs[0][j];
        }

        /** Get the coefficient S<sup>j</sup> for the function S.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the function S
         */
        double getSSj(const int j) const {
            return sjCoefs[0][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/da.
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/da
         */
        double getdSdaCj(const int j) const {
            return cjCoefs[1][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/da.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/da
         */
        double getdSdaSj(const int j) const {
            return sjCoefs[1][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/dk
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/dk
         */
        double getdSdkCj(const int j) const {
            return cjCoefs[2][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/dk.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/dk
         */
        double getdSdkSj(const int j) const {
            return sjCoefs[2][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/dh
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/dh
         */
        double getdSdhCj(const int j) const {
            return cjCoefs[3][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/dh.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/dh
         */
        double getdSdhSj(const int j) const {
            return sjCoefs[3][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/dα
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/dα
         */
        double getdSdalphaCj(const int j) const {
            return cjCoefs[4][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/dα.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/dα
         */
        double getdSdalphaSj(const int j) const {
            return sjCoefs[4][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/dβ
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/dβ
         */
        double getdSdbetaCj(const int j) const {
            return cjCoefs[5][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/dβ.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/dβ
         */
        double getdSdbetaSj(const int j) const {
            return sjCoefs[5][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/dγ
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/dγ
         */
        double getdSdgammaCj(const int j) const {
            return cjCoefs[6][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/dγ.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/dγ
         */
        double getdSdgammaSj(const int j) const {
            return sjCoefs[6][j];
        }

        /** Get the coefficient C<sup>j</sup> for the derivative dS/dλ
         * <br>
         * Possible values for j are within the interval [0,jMax].
         * The value 0 is used to obtain the free coefficient C⁰
         * @param j j index
         * @return C<sup>j</sup> for the function dS/dλ
         */
        double getdSdlambdaCj(const int j) const {
            return cjCoefs[7][j];
        }

        /** Get the coefficient S<sup>j</sup> for the derivative dS/dλ.
         * <br>
         * Possible values for j are within the interval [1,jMax].
         * @param j j index
         * @return S<sup>j</sup> for the derivative dS/dλ
         */
        double getdSdlambdaSj(const int j) const {
            return sjCoefs[7][j];
        }
    };



    //////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// ThirdBodyShortPeriodicCoefficients ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * The coefficients used to compute the short periodic contribution for the Third body perturbation.
     * <p>
     * The short periodic contribution for the Third Body is expressed in Danielson 4.2-25.<br>
     * The coefficients C<sub>i</sub>⁰, C<sub>i</sub><sup>j</sup>, S<sub>i</sub><sup>j</sup>
     * are computed by replacing the corresponding values in formula 2.5.5-10.
     * </p>
     */
    class ThirdBodyShortPeriodicCoefficients : public ShortPeriodTerms {

    private:

        /** Maximal value for j. */
        int jMax;

        /** Number of points used in the interpolation process. */
        int interpolationPoints;

        /** Max frequency of F. */
        int    maxFreqF;

        /** Coefficients prefix. */
        std::string prefix;

        /** All coefficients slotz. */
        TimeSpanMap<Slot> slotz;

    public:
        ThirdBodyShortPeriodicCoefficients() { }

        /**
         * Standard constructor.
         *  @param interpolationPoints number of points used in the interpolation process
         * @param jMax maximal value for j
         * @param maxFreqF Max frequency of F
         * @param bodyName third body name
         * @param slotz all coefficients slotz
         */
        ThirdBodyShortPeriodicCoefficients(const int jMax, const int interpolationPoints,
                                           const int maxFreqF, const std::string bodyName,
                                           const TimeSpanMap<Slot> slotz);

        /** Get the slot valid for some date.
         * @param meanStates mean states defining the slot
         * @return slot valid at the specified date
         */
        Slot createSlot(const std::vector< SpacecraftState > meanStates);

        /** {@inheritDoc} */
        SingleVectord value(const SpacecraftState spacecraftState);

        /** {@inheritDoc} */
        std::string getCoefficientsKeyPrefix() {
            return prefix;
        }

        /** {@inheritDoc}
         * <p>
         * For third body attraction forces,there are maxFreqF + 1 cj coefficients,
         * maxFreqF sj coefficients where maxFreqF depends on the orbit.
         * The j index is the integer multiplier for the eccentric longitude argument
         * in the cj and sj coefficients.
         * </p>
         */
        std::map<std::string, SingleVectord> getCoefficients(
                const AbsoluteDate date, const std::set<std::string> selected);

    private:
        /** Put a coefficient in a map if selected.
         * @param map map to populate
         * @param selected set of coefficients that should be put in the map
         * (empty set means all coefficients are selected)
         * @param value coefficient value
         * @param id coefficient identifier
         * @param indices list of coefficient indices
         */
        void storeIfSelected(std::map<std::string, SingleVectord> &map, const std::set<std::string> selected,
                             const SingleVectord value, const std::string id, const std::vector<int> indices);


        /** Internal class used only for serialization. */
        class DataTransferObject {

        private:
            /** Maximum value for j index. */
            int jMax;

            /** Number of points used in the interpolation process. */
            int interpolationPoints;

            /** Max frequency of F. */
            int    maxFreqF;

            /** Coefficients prefix. */
            std::string prefix;

            /** Transitions dates. */
            std::vector< AbsoluteDate > transitionDates;

            /** All slotz. */
            std::vector< Slot > allSlotz;

        public:
            /** Simple constructor.
              * @param jMax maximum value for j index
              * @param interpolationPoints number of points used in the interpolation process
              * @param maxFreqF max frequency of F
              * @param prefix prefix for coefficients keys
              * @param transitionDates transitions dates
              * @param allSlotz all slotz
              */
            DataTransferObject(const int jMax, const int interpolationPoints,
                               const int maxFreqF, const std::string prefix,
                               const std::vector< AbsoluteDate > transitionDates,
                               const std::vector< Slot > allSlotz) :
                jMax(jMax),
                interpolationPoints(interpolationPoints),
                maxFreqF(maxFreqF),
                prefix(prefix),
                transitionDates(transitionDates),
                allSlotz(allSlotz) {}

        private:
            /** Replace the deserialized data transfer object with a {@link ThirdBodyShortPeriodicCoefficients}.
              * @return replacement {@link ThirdBodyShortPeriodicCoefficients}
              */
            ThirdBodyShortPeriodicCoefficients readResolve();

        };


        /** Replace the instance with a data transfer object for serialization.
         * @return data transfer object that will be serialized
         * @exception NotSerializableException if an additional state provider is not serializable
         */
        DataTransferObject writeReplace();

    };


    void updateShortPeriodTerms( std::vector< SpacecraftState > meanStates );


private:
    //! Short period terms
    ThirdBodyShortPeriodicCoefficients shortPeriods;

};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_THIRDBODY_H
