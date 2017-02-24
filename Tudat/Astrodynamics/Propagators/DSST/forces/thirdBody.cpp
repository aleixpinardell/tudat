/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "thirdBody.h"
#include "Tudat/Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "tudat/Tudat/Astrodynamics/Propagators/DSST/utilities/derivativeStructure.h"
#include "tudat/Tudat/Astrodynamics/Propagators/DSST/utilities/jacobiPolynomials.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

//! Constructor.
DSSTThirdBody::DSSTThirdBody( const std::string bodyName, const Vector3 bodyPosition,
                              const double gravitationalParameter /*, const double distance*/ ):
    bodyName(bodyName), bodyPosition(bodyPosition), gm(gravitationalParameter), R3(bodyPosition.norm()),
    maxAR3Pow(INT_MIN), maxEccPow(INT_MIN)
{
    Vns = CoefficientsFactory().computeVns( MAX_POWER );
}


std::vector<ShortPeriodTerms> DSSTThirdBody::initialize( AuxiliaryElements auxiliaryElements, const bool meanOnly )
{
    // Initializes specific parameters.
    initializeStep( auxiliaryElements );

    // Truncation tolerance.
    const double aor = a / R3;
    const double tol = ( aor > 0.3 || ( aor > 0.15  && ecc > 0.25 ) ) ?
                BIG_TRUNCATION_TOLERANCE : SMALL_TRUNCATION_TOLERANCE;

    // Utilities for truncation
    // Set a lower bound for eccentricity
    const double eo2  = std::max( 0.0025, ecc / 2.0 );
    const double x2o2 = XX / 2.;
    std::vector< double > eccPwr( MAX_POWER );
    std::vector< double > chiPwr( MAX_POWER );
    eccPwr[0] = 1.0;
    chiPwr[0] = X;
    for ( int i = 1; i < MAX_POWER; i++ ) {
        eccPwr[i] = eccPwr[i - 1] * eo2;
        chiPwr[i] = chiPwr[i - 1] * x2o2;
    }

    // Auxiliary quantities.
    const double ao2rxx = aor / (2.0 * XX);
    double xmuarn       = ao2rxx * ao2rxx * gm / (X * R3);
    double term         = 0.0;

    // Compute max power for a/R3 and e.
    maxAR3Pow = 2;
    maxEccPow = 0;
    int n     = 2;
    int m     = 2;
    int nsmd2 = 0;

    do {
        // Upper bound for Tnm.
        term =  xmuarn *
                ( factorial( n + m ) / ( factorial( nsmd2 ) * factorial( nsmd2 + m ) ) ) *
                ( factorial( n + m + 1 ) / ( factorial( m ) * factorial( n + 1 ) ) ) *
                ( factorial( n - m + 1 ) / factorial( n + 1 ) ) *
                eccPwr[ m ] * upper_bounds::getDnl( XX, chiPwr[ m ], n + 2, m );

        if ( term < tol ) {
            if (m == 0) {
                break;
            } else if ( m < 2 ) {
                xmuarn *= ao2rxx;
                m = 0;
                n++;
                nsmd2++;
            } else {
                m -= 2;
                nsmd2++;
            }
        } else {
            maxAR3Pow = n;
            maxEccPow = std::max( m, maxEccPow );
            xmuarn *= ao2rxx;
            m++;
            n++;
        }
    } while ( n < MAX_POWER );

    maxEccPow = std::min( maxAR3Pow, maxEccPow );

    // allocate the array aoR3Pow
    aoR3Pow = std::vector< double >( maxAR3Pow + 1 );

    maxFreqF = maxAR3Pow + 1;
    maxEccPowShort = MAX_ECCPOWER_SP;

    Qns = CoefficientsFactory().computeQns( gamma, maxAR3Pow, std::max( maxEccPow, maxEccPowShort ) );

    const int jMax = maxAR3Pow + 1;

    shortPeriods = ThirdBodyShortPeriodicCoefficients(
                jMax, INTERPOLATION_POINTS, maxFreqF, bodyName, TimeSpanMap<Slot>(
                    Slot(jMax,INTERPOLATION_POINTS ) ) );

    std::vector<ShortPeriodTerms> list;
    list.push_back(shortPeriods);
    return list;
}


void DSSTThirdBody::initializeStep( const AuxiliaryElements auxiliaryElements )
{
    // Equinoctial elements
    a = auxiliaryElements.sma;
    k = auxiliaryElements.k;
    h = auxiliaryElements.h;
    q = auxiliaryElements.q;
    p = auxiliaryElements.p;

    // Eccentricity
    ecc = auxiliaryElements.ecc;

    // Distance from center of mass of the central body to the 3rd body
    R3 = bodyPosition.norm();

    // Direction cosines
    const Vector3 bodyDir = bodyPosition / R3;
    alpha = bodyDir.dot( auxiliaryElements.f );
    beta  = bodyDir.dot( auxiliaryElements.g );
    gamma = bodyDir.dot( auxiliaryElements.w );

    // Equinoctial coefficients
    A = auxiliaryElements.A;
    B = auxiliaryElements.B;
    C = auxiliaryElements.C;
    meanMotion = auxiliaryElements.n;

    //\Chi^{-2}
    BB = B * B;
    //\Chi^{-3}
    BBB = BB * B;

    //b = 1 / (1 + B)
    b = 1.0 / (1.0 + B);

    // \Chi
    X = 1.0 / B;
    XX = X * X;
    XXX = X * XX;
    // -2 * a / A
    m2aoA = -2.0 * a / A;
    // B / A
    BoA = B / A;
    // 1 / AB
    ooAB = 1.0 / (A * B);
    // -C / 2AB
    mCo2AB = -C * ooAB / 2.0;
    // B / A(1 + B)
    BoABpo = BoA / (1.0 + B);

    // mu3 / R3
    muoR3 = gm / R3;

    //h * \Chi³
    hXXX = h * XXX;
    //k * \Chi³
    kXXX = k * XXX;
}



v_double DSSTThirdBody::getMeanElementRates( )
{
    // Qns coefficients
    Qns = CoefficientsFactory().computeQns( gamma, maxAR3Pow, maxEccPow );

    // a / R3 up to power maxAR3Pow
    const double aoR3 = a / R3;
    aoR3Pow[ 0 ] = 1.0;
    for ( int i = 1; i <= maxAR3Pow; i++ ) {
        aoR3Pow[ i ] = aoR3 * aoR3Pow[ i - 1 ];
    }

    // Compute potential U derivatives
    const v_double dU  = computeUDerivatives();
    const double dUda  = dU[0];
    const double dUdk  = dU[1];
    const double dUdh  = dU[2];
    const double dUdAl = dU[3];
    const double dUdBe = dU[4];
    const double dUdGa = dU[5];

    // Compute cross derivatives [Eq. 2.2-(8)]
    // U(alpha,gamma) = alpha * dU/dgamma - gamma * dU/dalpha
    const double UAlphaGamma   = alpha * dUdGa - gamma * dUdAl;
    // U(beta,gamma) = beta * dU/dgamma - gamma * dU/dbeta
    const double UBetaGamma    =  beta * dUdGa - gamma * dUdBe;
    // Common factor
    const double pUAGmIqUBGoAB = ( p * UAlphaGamma - I * q * UBetaGamma ) * ooAB;

    // Compute mean elements rates [Eq. 3.1-(1)]
    v_double meanElementRates;
    meanElementRates = { 0.0,                                                               // da
                         BoA * dUdk + k * pUAGmIqUBGoAB,                                    // dk
                         -BoA * dUdh - h * pUAGmIqUBGoAB,                                    // dh
                         mCo2AB * UBetaGamma,                                               // dq
                         mCo2AB * UAlphaGamma * I,                                          // dp
                         m2aoA * dUda + BoABpo * ( h * dUdh + k * dUdk ) + pUAGmIqUBGoAB }; // dlambda
    return meanElementRates;
}


v_double DSSTThirdBody::computeUDerivatives()
{
    // Gs and Hs coefficients
    Eigen::Matrix2Xd GsHs = CoefficientsFactory().computeGsHs( k, h, alpha, beta, maxEccPow );

    // Gs partial derivatives
    std::map< std::string, Eigen::VectorXd > dGs =
            CoefficientsFactory().computeGsDerivatives( k, h, alpha, beta, maxEccPow );

    // Initialise U.
    U = 0.0;

    // Potential derivatives
    double dUda  = 0.0;
    double dUdk  = 0.0;
    double dUdh  = 0.0;
    double dUdAl = 0.0;
    double dUdBe = 0.0;
    double dUdGa = 0.0;

    // Get kernels of Hansen coefficients and their derivatives
    const NSSetPair K0nsAndDerivatives = CoefficientsFactory().computeK0nsAndDerivatives( X, maxAR3Pow, maxEccPow );
    const vv_double K0ns  = K0nsAndDerivatives.first;
    const vv_double dK0ns = K0nsAndDerivatives.second;

    for ( int s = 0; s <= maxEccPow; s++ ) {
        // Get the current Gs coefficient
        const double gs = GsHs( 0, s );

        // Compute Gs partial derivatives
        double dGsdh  = dGs[ "h" ]( s );
        double dGsdk  = dGs[ "k" ]( s );
        double dGsdAl = dGs[ "alpha" ]( s );
        double dGsdBe = dGs[ "beta" ]( s );

        // Kronecker symbol (2 - delta(0,s))
        const double delta0s = ( s == 0 ) ? 1.0 : 2.0;

        for ( int n = std::max( 2, s ); n <= maxAR3Pow; n++ ) {
            // (n - s) must be even
            if ( (n - s) % 2 == 0 ) {
                // Extract data from previous computation :
                const double kns   = K0ns[n][s];
                const double dkns  = dK0ns[n][s];

                const double vns   = Vns[ NSKey(n, s) ];
                const double coef0 = delta0s * aoR3Pow[n] * vns;
                const double coef1 = coef0 * Qns[n][s];
                const double coef2 = coef1 * kns;
                // dQns/dGamma = Q(n, s + 1) from Equation 3.1-(8)
                // for n = s, Q(n, n + 1) = 0. (Cefola & Broucke, 1975)
                const double dqns = ( n == s ) ? 0. : Qns[n][s + 1];

                //Compute U:
                U += coef2 * gs;

                // Compute dU / da :
                dUda += coef2 * n * gs;
                // Compute dU / dh
                dUdh += coef1 * (kns * dGsdh + hXXX * gs * dkns);
                // Compute dU / dk
                dUdk += coef1 * (kns * dGsdk + kXXX * gs * dkns);
                // Compute dU / dAlpha
                dUdAl += coef2 * dGsdAl;
                // Compute dU / dBeta
                dUdBe += coef2 * dGsdBe;
                // Compute dU / dGamma
                dUdGa += coef0 * kns * dqns * gs;
            }
        }
    }

    // multiply by mu3 / R3
    U *= muoR3;

    // Return U derivatives
    v_double UDerivatives = { dUda  * muoR3 / a,
                             dUdk  * muoR3,
                             dUdh  * muoR3,
                             dUdAl * muoR3,
                             dUdBe * muoR3,
                             dUdGa * muoR3 };
    return UDerivatives;
}


void DSSTThirdBody::updateShortPeriodTerms( std::vector< SpacecraftState > meanStates )
{
    Slot slot = shortPeriods.createSlot(meanStates);

    for (const SpacecraftState meanState : meanStates) {

        initializeStep( AuxiliaryElements(meanState, I) );

        // a / R3 up to power maxAR3Pow
        const double aoR3 = a / R3;
        aoR3Pow[0] = 1.;
        for (int i = 1; i <= maxAR3Pow; i++) {
            aoR3Pow[i] = aoR3 * aoR3Pow[i - 1];
        }

        // Qns coefficients
        Qns = CoefficientsFactory().computeQns(gamma, maxAR3Pow, std::max(maxEccPow, maxEccPowShort));
        const GeneratingFunctionCoefficients gfCoefs =
                GeneratingFunctionCoefficients(*this, maxAR3Pow, MAX_ECCPOWER_SP, maxAR3Pow + 1);

        //Compute additional quantities
        // 2 * a / An
        const double ax2oAn  = -m2aoA / meanMotion;
        // B / An
        const double BoAn  = BoA / meanMotion;
        // 1 / ABn
        const double ooABn = ooAB / meanMotion;
        // C / 2ABn
        const double Co2ABn = -mCo2AB / meanMotion;
        // B / (A * (1 + B) * n)
        const double BoABpon = BoABpo / meanMotion;
        // -3 / n²a² = -3 / nA
        const double m3onA = -3 / (A * meanMotion);

        //Compute the C<sub>i</sub><sup>j</sup> and S<sub>i</sub><sup>j</sup> coefficients.
        for (unsigned int j = 1; j < slot.cij.size(); j++) {
            // First compute the C<sub>i</sub><sup>j</sup> coefficients
            v_double currentCij(6);

            // Compute the cross derivatives operator :
            const double SAlphaGammaCj   = alpha * gfCoefs.getdSdgammaCj(j) - gamma * gfCoefs.getdSdalphaCj(j);
            const double SAlphaBetaCj    = alpha * gfCoefs.getdSdbetaCj(j)  - beta  * gfCoefs.getdSdalphaCj(j);
            const double SBetaGammaCj    =  beta * gfCoefs.getdSdgammaCj(j) - gamma * gfCoefs.getdSdbetaCj(j);
            const double ShkCj           =     h * gfCoefs.getdSdkCj(j)     -  k    * gfCoefs.getdSdhCj(j);
            const double pSagmIqSbgoABnCj = (p * SAlphaGammaCj - I * q * SBetaGammaCj) * ooABn;
            const double ShkmSabmdSdlCj  =  ShkCj - SAlphaBetaCj - gfCoefs.getdSdlambdaCj(j);

            currentCij[0] =  ax2oAn * gfCoefs.getdSdlambdaCj(j);
            currentCij[1] =  -(BoAn * gfCoefs.getdSdhCj(j) + h * pSagmIqSbgoABnCj + k * BoABpon * gfCoefs.getdSdlambdaCj(j));
            currentCij[2] =    BoAn * gfCoefs.getdSdkCj(j) + k * pSagmIqSbgoABnCj - h * BoABpon * gfCoefs.getdSdlambdaCj(j);
            currentCij[3] =  Co2ABn * (q * ShkmSabmdSdlCj - I * SAlphaGammaCj);
            currentCij[4] =  Co2ABn * (p * ShkmSabmdSdlCj - SBetaGammaCj);
            currentCij[5] = -ax2oAn * gfCoefs.getdSdaCj(j) + BoABpon * (h * gfCoefs.getdSdhCj(j) + k * gfCoefs.getdSdkCj(j)) + pSagmIqSbgoABnCj + m3onA * gfCoefs.getSCj(j);

            // add the computed coefficients to the interpolators
            slot.cij[j].addGridPoint(meanState.getDate(), currentCij);

            // Compute the S<sub>i</sub><sup>j</sup> coefficients
            v_double currentSij(6);

            // Compute the cross derivatives operator :
            const double SAlphaGammaSj   = alpha * gfCoefs.getdSdgammaSj(j) - gamma * gfCoefs.getdSdalphaSj(j);
            const double SAlphaBetaSj    = alpha * gfCoefs.getdSdbetaSj(j)  - beta  * gfCoefs.getdSdalphaSj(j);
            const double SBetaGammaSj    =  beta * gfCoefs.getdSdgammaSj(j) - gamma * gfCoefs.getdSdbetaSj(j);
            const double ShkSj           =     h * gfCoefs.getdSdkSj(j)     -  k    * gfCoefs.getdSdhSj(j);
            const double pSagmIqSbgoABnSj = (p * SAlphaGammaSj - I * q * SBetaGammaSj) * ooABn;
            const double ShkmSabmdSdlSj  =  ShkSj - SAlphaBetaSj - gfCoefs.getdSdlambdaSj(j);

            currentSij[0] =  ax2oAn * gfCoefs.getdSdlambdaSj(j);
            currentSij[1] =  -(BoAn * gfCoefs.getdSdhSj(j) + h * pSagmIqSbgoABnSj + k * BoABpon * gfCoefs.getdSdlambdaSj(j));
            currentSij[2] =    BoAn * gfCoefs.getdSdkSj(j) + k * pSagmIqSbgoABnSj - h * BoABpon * gfCoefs.getdSdlambdaSj(j);
            currentSij[3] =  Co2ABn * (q * ShkmSabmdSdlSj - I * SAlphaGammaSj);
            currentSij[4] =  Co2ABn * (p * ShkmSabmdSdlSj - SBetaGammaSj);
            currentSij[5] = -ax2oAn * gfCoefs.getdSdaSj(j) + BoABpon * (h * gfCoefs.getdSdhSj(j) + k * gfCoefs.getdSdkSj(j)) + pSagmIqSbgoABnSj + m3onA * gfCoefs.getSSj(j);

            // add the computed coefficients to the interpolators
            slot.sij[j].addGridPoint(meanState.getDate(), currentSij);

            if (j == 1) {
                //Compute the C⁰ coefficients using Danielson 2.5.2-15a.
                v_double value(6);
                for (int i = 0; i < 6; ++i) {
                    value[i] = currentCij[i] * k / 2. + currentSij[i] * h / 2.;
                }
                slot.cij[0].addGridPoint(meanState.getDate(), value);
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// FourierCjSjCoefficients ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

DSSTThirdBody::FourierCjSjCoefficients::FourierCjSjCoefficients(
        DSSTThirdBody &owner, const int nMax, const int sMax, const int jMax ) :
    owner(&owner), gns(GnsCoefficients(owner, nMax, sMax)), wnsjEtomjmsCoefficient(WnsjEtomjmsCoefficient(owner)),
    ABDECoefficients(CjSjAlphaBetaKH(owner)), cj(new_vv_double(7, jMax + 1)), sj(new_vv_double(7, jMax + 1)),
    cjlambda(new_v_double(jMax)), sjlambda(new_v_double(jMax)), nMax(nMax), sMax(sMax), jMax(jMax)
{
    computeCoefficients();
}

void DSSTThirdBody::FourierCjSjCoefficients::computeCoefficients()
{
    for (int j = 1; j <= jMax; j++) {
        // initialise the coefficients
        for (int i = 0; i <= 6; i++) {
            cj[i][j] = 0.;
            sj[i][j] = 0.;
        }
        if (j < jMax) {
            // initialise the C<sup>j</sup><sub>,λ</sub> and S<sup>j</sup><sub>,λ</sub> coefficients
            cjlambda[j] = 0.;
            sjlambda[j] = 0.;
        }
        for (int s = 0; s <= sMax; s++) {

            // Compute the coefficients A<sub>j, s</sub>, B<sub>j, s</sub>, D<sub>j, s</sub> and E<sub>j, s</sub>
            ABDECoefficients.computeCoefficients(j, s);

            // compute starting value for n
            const int minN = std::max(2, std::max(j - 1, s));

            for (int n = minN; n <= nMax; n++) {
                // check if n-s is even
                if ((n - s) % 2 == 0) {
                    // compute the coefficient e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n+1, s</sup>
                    // and its derivatives
                    const v_double wjnp1semjms =
                            wnsjEtomjmsCoefficient.computeWjnsEmjmsAndDeriv(j, s, n + 1);

                    // compute the coefficient e<sup>-|j-s|</sup>*w<sub>-j</sub><sup>n+1, s</sup>
                    // and its derivatives
                    const v_double wmjnp1semjms =
                            wnsjEtomjmsCoefficient.computeWjnsEmjmsAndDeriv(-j, s, n + 1);

                    // compute common factors
                    const double coef1 = -(wjnp1semjms[0] * ABDECoefficients.getCoefA()
                            + wmjnp1semjms[0] * ABDECoefficients.getCoefB());
                    const double coef2 =   wjnp1semjms[0] * ABDECoefficients.getCoefD()
                            + wmjnp1semjms[0] * ABDECoefficients.getCoefE();

                    //Compute C<sup>j</sup>
                    cj[0][j] += gns.getGns(n, s) * coef1;
                    //Compute dC<sup>j</sup> / da
                    cj[1][j] += gns.getdGnsda(n, s) * coef1;
                    //Compute dC<sup>j</sup> / dk
                    cj[2][j] += -gns.getGns(n, s) *
                            ( wjnp1semjms[1] * ABDECoefficients.getCoefA() +
                            wjnp1semjms[0] * ABDECoefficients.getdCoefAdk() +
                            wmjnp1semjms[1] * ABDECoefficients.getCoefB() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefBdk()
                            );
                    //Compute dC<sup>j</sup> / dh
                    cj[3][j] += -gns.getGns(n, s) *
                            ( wjnp1semjms[2] * ABDECoefficients.getCoefA() +
                            wjnp1semjms[0] * ABDECoefficients.getdCoefAdh() +
                            wmjnp1semjms[2] * ABDECoefficients.getCoefB() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefBdh()
                            );
                    //Compute dC<sup>j</sup> / dα
                    cj[4][j] += -gns.getGns(n, s) *
                            ( wjnp1semjms[0] * ABDECoefficients.getdCoefAdalpha() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefBdalpha()
                            );
                    //Compute dC<sup>j</sup> / dβ
                    cj[5][j] += -gns.getGns(n, s) *
                            ( wjnp1semjms[0] * ABDECoefficients.getdCoefAdbeta() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefBdbeta()
                            );
                    //Compute dC<sup>j</sup> / dγ
                    cj[6][j] += gns.getdGnsdgamma(n, s) * coef1;

                    //Compute S<sup>j</sup>
                    sj[0][j] += gns.getGns(n, s) * coef2;
                    //Compute dS<sup>j</sup> / da
                    sj[1][j] += gns.getdGnsda(n, s) * coef2;
                    //Compute dS<sup>j</sup> / dk
                    sj[2][j] += gns.getGns(n, s) *
                            ( wjnp1semjms[1] * ABDECoefficients.getCoefD() +
                            wjnp1semjms[0] * ABDECoefficients.getdCoefDdk() +
                            wmjnp1semjms[1] * ABDECoefficients.getCoefE() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefEdk()
                            );
                    //Compute dS<sup>j</sup> / dh
                    sj[3][j] += gns.getGns(n, s) *
                            ( wjnp1semjms[2] * ABDECoefficients.getCoefD() +
                            wjnp1semjms[0] * ABDECoefficients.getdCoefDdh() +
                            wmjnp1semjms[2] * ABDECoefficients.getCoefE() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefEdh()
                            );
                    //Compute dS<sup>j</sup> / dα
                    sj[4][j] += gns.getGns(n, s) *
                            ( wjnp1semjms[0] * ABDECoefficients.getdCoefDdalpha() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefEdalpha()
                            );
                    //Compute dS<sup>j</sup> / dβ
                    sj[5][j] += gns.getGns(n, s) *
                            ( wjnp1semjms[0] * ABDECoefficients.getdCoefDdbeta() +
                            wmjnp1semjms[0] * ABDECoefficients.getdCoefEdbeta()
                            );
                    //Compute dS<sup>j</sup> / dγ
                    sj[6][j] += gns.getdGnsdgamma(n, s) * coef2;

                    //Check if n is greater or equal to j and j is at most jMax-1
                    if (n >= j && j < jMax) {
                        // compute the coefficient e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup>
                        // and its derivatives
                        const v_double wjnsemjms = wnsjEtomjmsCoefficient.computeWjnsEmjmsAndDeriv(j, s, n);

                        // compute the coefficient e<sup>-|j-s|</sup>*w<sub>-j</sub><sup>n, s</sup>
                        // and its derivatives
                        const v_double wmjnsemjms = wnsjEtomjmsCoefficient.computeWjnsEmjmsAndDeriv(-j, s, n);

                        //Compute C<sup>j</sup><sub>,λ</sub>
                        cjlambda[j] += gns.getGns(n, s) * (wjnsemjms[0] * ABDECoefficients.getCoefD() +
                                wmjnsemjms[0] * ABDECoefficients.getCoefE());
                        //Compute S<sup>j</sup><sub>,λ</sub>
                        sjlambda[j] += gns.getGns(n, s) * (wjnsemjms[0] * ABDECoefficients.getCoefA() +
                                wmjnsemjms[0] * ABDECoefficients.getCoefB());
                    }
                }
            }
        }
        // Divide by j
        for (int i = 0; i <= 6; i++) {
            cj[i][j] /= j;
            sj[i][j] /= j;
        }
    }
    //The C⁰ coefficients are not computed here.
    //They are evaluated at the final point.

    //C⁰<sub>,λ</sub>
    cjlambda[0] = owner->k * cjlambda[1] / 2. + owner->h * sjlambda[1] / 2.;
}



//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// WnsjEtomjmsCoefficient ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

DSSTThirdBody::WnsjEtomjmsCoefficient::WnsjEtomjmsCoefficient( DSSTThirdBody &owner ) : owner(&owner)
{
    //initialise fields
    c = owner.ecc * owner.b;
    c2 = c * c;

    //b² * χ
    const double b2Chi = owner.b * owner.b * owner.X;
    //Compute derivatives of b
    dbdh = owner.h * b2Chi;
    dbdk = owner.k * b2Chi;

    //Compute derivatives of e
    dedh = owner.h / owner.ecc;
    dedk = owner.k / owner.ecc;

    //Compute derivatives of c
    dcdh = owner.ecc * dbdh + owner.b * dedh;
    dcdk = owner.ecc * dbdk + owner.b * dedk;

    //Compute the powers (1 - c²)<sup>n</sup> and (1 + c²)<sup>n</sup>
    omc2tn = new_v_double(owner.maxAR3Pow + owner.maxFreqF + 2);
    opc2tn = new_v_double(owner.maxAR3Pow + owner.maxFreqF + 2);
    const double omc2 = 1. - c2;
    const double opc2 = 1. + c2;
    omc2tn[0] = 1.;
    opc2tn[0] = 1.;
    for (int i = 1; i <= owner.maxAR3Pow + owner.maxFreqF + 1; i++) {
        omc2tn[i] = omc2tn[i - 1] * omc2;
        opc2tn[i] = opc2tn[i - 1] * opc2;
    }

    //Compute the powers of b
    btjms = new_v_double(owner.maxAR3Pow + owner.maxFreqF + 1);
    btjms[0] = 1.;
    for (int i = 1; i <= owner.maxAR3Pow + owner.maxFreqF; i++) {
        btjms[i] = btjms[i - 1] * owner.b;
    }
}

v_double DSSTThirdBody::WnsjEtomjmsCoefficient::computeWjnsEmjmsAndDeriv( const int j, const int s, const int n )
{
    v_double wjnsemjms = {0., 0., 0.};

    // |j|
    const int absJ = std::abs(j);
    // |s|
    const int absS = std::abs(s);
    // |j - s|
    const int absJmS = std::abs(j - s);
    // |j + s|
    const int absJpS = std::abs(j + s);

    //The lower index of P. Also the power of (1 - c²)
    int l;
    // the factorial ratio coefficient or 1. if |s| <= |j|
    double factCoef;
    if (absS > absJ) {
        factCoef = (factorial(n + s) / factorial(n + j)) * (factorial(n - s) / factorial(n - j));
        l = n - absS;
    } else {
        factCoef = 1.;
        l = n - absJ;
    }

    // (-1)<sup>|j-s|</sup>
    const double sign = absJmS % 2 != 0 ? -1. : 1.;
    //(1 - c²)<sup>n-|s|</sup> / (1 + c²)<sup>n</sup>
    const double coef1 = omc2tn[l] / opc2tn[n];
    //-b<sup>|j-s|</sup>
    const double coef2 = sign * btjms[absJmS];
    // P<sub>l</sub><sup>|j-s|, |j+s|</sup>(χ)
    const DerivativeStructure jac =
            JacobiPolynomials::getValue(l, absJmS, absJpS, DerivativeStructure(1, 1, 0, owner->X));

    // the derivative of coef1 by c
    const double dcoef1dc = -coef1 * 2. * c * (((double) n) / opc2tn[1] + ((double) l) / omc2tn[1]);
    // the derivative of coef1 by h
    const double dcoef1dh = dcoef1dc * dcdh;
    // the derivative of coef1 by k
    const double dcoef1dk = dcoef1dc * dcdk;

    // the derivative of coef2 by b
    const double dcoef2db = absJmS == 0 ? 0 : sign * (double) absJmS * btjms[absJmS - 1];
    // the derivative of coef2 by h
    const double dcoef2dh = dcoef2db * dbdh;
    // the derivative of coef2 by k
    const double dcoef2dk = dcoef2db * dbdk;

    // the jacobi polinomial value
    const double jacobi = jac.getValue();
    // the derivative of the Jacobi polynomial by h
    const double djacobidh = jac.getPartialDerivative({1}) * owner->hXXX;
    // the derivative of the Jacobi polynomial by k
    const double djacobidk = jac.getPartialDerivative({1}) * owner->kXXX;

    //group the above coefficients to limit the mathematical operations
    const double term1 = factCoef * coef1 * coef2;
    const double term2 = factCoef * coef1 * jacobi;
    const double term3 = factCoef * coef2 * jacobi;

    //compute e<sup>-|j-s|</sup>*w<sub>j</sub><sup>n, s</sup> and its derivatives by k and h
    wjnsemjms[0] = term1 * jacobi;
    wjnsemjms[1] = dcoef1dk * term3 + dcoef2dk * term2 + djacobidk * term1;
    wjnsemjms[2] = dcoef1dh * term3 + dcoef2dh * term2 + djacobidh * term1;

    return wjnsemjms;
}




///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// GnsCoefficients ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////


DSSTThirdBody::GnsCoefficients::GnsCoefficients(DSSTThirdBody &owner, const int nMax, const int sMax) : owner(&owner)
{
    this->nMax = nMax;
    this->sMax = sMax;

    const int rows    = nMax + 1;
    const int columns = sMax + 1;
    this->gns          = new_vv_double(rows,columns);
    this->dgnsda       = new_vv_double(rows,columns);
    this->dgnsdgamma   = new_vv_double(rows,columns);

    // Generate the coefficients
    generateCoefficients();
}


void DSSTThirdBody::GnsCoefficients::generateCoefficients()
{
    for (int s = 0; s <= sMax; s++) {
        // The n index is always at least the maximum between 2 and s
        const int minN = std::max(2, s);
        for (int n = minN; n <= nMax; n++) {
            // compute the coefficients only if (n - s) % 2 == 0
            if ( (n - s) % 2 == 0 ) {
                // Kronecker symbol (2 - delta(0,s))
                const double delta0s = (s == 0) ? 1. : 2.;
                const double vns   = owner->Vns[ NSKey(n, s) ];
                const double coef0 = delta0s * owner->aoR3Pow[n] * vns * owner->muoR3;
                const double coef1 = coef0 * owner->Qns[n][s];
                // dQns/dGamma = Q(n, s + 1) from Equation 3.1-(8)
                // for n = s, Q(n, n + 1) = 0. (Cefola & Broucke, 1975)
                const double dqns = (n == s) ? 0. : owner->Qns[n][s + 1];

                //Compute the coefficient and its derivatives.
                this->gns[n][s] = coef1;
                this->dgnsda[n][s] = coef1 * n / owner->a;
                this->dgnsdgamma[n][s] = coef0 * dqns;
            } else {
                // the coefficient and its derivatives is 0
                this->gns[n][s] = 0.;
                this->dgnsda[n][s] = 0.;
                this->dgnsdgamma[n][s] = 0.;
            }
        }
    }
}




///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// CjSjAlphaBetaKH ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void DSSTThirdBody::CjSjAlphaBetaKH::computeCoefficients(const int j, const int s)
{
    // sign of j-s
    const int sign = j < s ? -1 : 1;

    //|j-s|
    const int absJmS = std::abs(j - s);

    //j+s
    const int jps = j + s;

    //Compute the coefficient A and its derivatives
    coefAandDeriv[0] = sign * cjsjalbe.getCj(s) * cjsjkh.getSj(absJmS) + cjsjalbe.getSj(s) * cjsjkh.getCj(absJmS);
    coefAandDeriv[1] = sign * cjsjalbe.getCj(s) * cjsjkh.getDsjDk(absJmS) + cjsjalbe.getSj(s) * cjsjkh.getDcjDk(absJmS);
    coefAandDeriv[2] = sign * cjsjalbe.getCj(s) * cjsjkh.getDsjDh(absJmS) + cjsjalbe.getSj(s) * cjsjkh.getDcjDh(absJmS);
    coefAandDeriv[3] = sign * cjsjalbe.getDcjDk(s) * cjsjkh.getSj(absJmS) + cjsjalbe.getDsjDk(s) * cjsjkh.getCj(absJmS);
    coefAandDeriv[4] = sign * cjsjalbe.getDcjDh(s) * cjsjkh.getSj(absJmS) + cjsjalbe.getDsjDh(s) * cjsjkh.getCj(absJmS);

    //Compute the coefficient B and its derivatives
    coefBandDeriv[0] = cjsjalbe.getCj(s) * cjsjkh.getSj(jps) - cjsjalbe.getSj(s) * cjsjkh.getCj(jps);
    coefBandDeriv[1] = cjsjalbe.getCj(s) * cjsjkh.getDsjDk(jps) - cjsjalbe.getSj(s) * cjsjkh.getDcjDk(jps);
    coefBandDeriv[2] = cjsjalbe.getCj(s) * cjsjkh.getDsjDh(jps) - cjsjalbe.getSj(s) * cjsjkh.getDcjDh(jps);
    coefBandDeriv[3] = cjsjalbe.getDcjDk(s) * cjsjkh.getSj(jps) - cjsjalbe.getDsjDk(s) * cjsjkh.getCj(jps);
    coefBandDeriv[4] = cjsjalbe.getDcjDh(s) * cjsjkh.getSj(jps) - cjsjalbe.getDsjDh(s) * cjsjkh.getCj(jps);

    //Compute the coefficient D and its derivatives
    coefDandDeriv[0] = cjsjalbe.getCj(s) * cjsjkh.getCj(absJmS) - sign * cjsjalbe.getSj(s) * cjsjkh.getSj(absJmS);
    coefDandDeriv[1] = cjsjalbe.getCj(s) * cjsjkh.getDcjDk(absJmS) - sign * cjsjalbe.getSj(s) * cjsjkh.getDsjDk(absJmS);
    coefDandDeriv[2] = cjsjalbe.getCj(s) * cjsjkh.getDcjDh(absJmS) - sign * cjsjalbe.getSj(s) * cjsjkh.getDsjDh(absJmS);
    coefDandDeriv[3] = cjsjalbe.getDcjDk(s) * cjsjkh.getCj(absJmS) - sign * cjsjalbe.getDsjDk(s) * cjsjkh.getSj(absJmS);
    coefDandDeriv[4] = cjsjalbe.getDcjDh(s) * cjsjkh.getCj(absJmS) - sign * cjsjalbe.getDsjDh(s) * cjsjkh.getSj(absJmS);

    //Compute the coefficient E and its derivatives
    coefEandDeriv[0] = cjsjalbe.getCj(s) * cjsjkh.getCj(jps) + cjsjalbe.getSj(s) * cjsjkh.getSj(jps);
    coefEandDeriv[1] = cjsjalbe.getCj(s) * cjsjkh.getDcjDk(jps) + cjsjalbe.getSj(s) * cjsjkh.getDsjDk(jps);
    coefEandDeriv[2] = cjsjalbe.getCj(s) * cjsjkh.getDcjDh(jps) + cjsjalbe.getSj(s) * cjsjkh.getDsjDh(jps);
    coefEandDeriv[3] = cjsjalbe.getDcjDk(s) * cjsjkh.getCj(jps) + cjsjalbe.getDsjDk(s) * cjsjkh.getSj(jps);
    coefEandDeriv[4] = cjsjalbe.getDcjDh(s) * cjsjkh.getCj(jps) + cjsjalbe.getDsjDh(s) * cjsjkh.getSj(jps);
}


//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// GeneratingFunctionCoefficients ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

DSSTThirdBody::GeneratingFunctionCoefficients::GeneratingFunctionCoefficients(
        DSSTThirdBody &owner, const int nMax, const int sMax, const int jMax) :
    owner(&owner), jMax(jMax), cjsjFourier(FourierCjSjCoefficients(owner, nMax, sMax, jMax)),
    cjCoefs(new_vv_double(8,jMax + 1)), sjCoefs(new_vv_double(8,jMax + 1))
{
     computeGeneratingFunctionCoefficients();
}

void DSSTThirdBody::GeneratingFunctionCoefficients::computeGeneratingFunctionCoefficients()
{
    // Compute potential U and derivatives
    const v_double dU  = owner->computeUDerivatives();

    //Compute the C<sup>j</sup> coefficients
    for (int j = 1; j <= jMax; j++) {
        //Compute the C<sup>j</sup> coefficients
        cjCoefs[0][j] = cjsjFourier.getCj(j);
        cjCoefs[1][j] = cjsjFourier.getdCjda(j);
        cjCoefs[2][j] = cjsjFourier.getdCjdk(j) - (cjsjFourier.getSjLambda(j - 1) - cjsjFourier.getSjLambda(j + 1)) / 2;
        cjCoefs[3][j] = cjsjFourier.getdCjdh(j) - (cjsjFourier.getCjLambda(j - 1) + cjsjFourier.getCjLambda(j + 1)) / 2;
        cjCoefs[4][j] = cjsjFourier.getdCjdalpha(j);
        cjCoefs[5][j] = cjsjFourier.getdCjdbeta(j);
        cjCoefs[6][j] = cjsjFourier.getdCjdgamma(j);
        cjCoefs[7][j] = cjsjFourier.getCjLambda(j);

        //Compute the S<sup>j</sup> coefficients
        sjCoefs[0][j] = cjsjFourier.getSj(j);
        sjCoefs[1][j] = cjsjFourier.getdSjda(j);
        sjCoefs[2][j] = cjsjFourier.getdSjdk(j) + (cjsjFourier.getCjLambda(j - 1) - cjsjFourier.getCjLambda(j + 1)) / 2;
        sjCoefs[3][j] = cjsjFourier.getdSjdh(j) - (cjsjFourier.getSjLambda(j - 1) + cjsjFourier.getSjLambda(j + 1)) / 2;
        sjCoefs[4][j] = cjsjFourier.getdSjdalpha(j);
        sjCoefs[5][j] = cjsjFourier.getdSjdbeta(j);
        sjCoefs[6][j] = cjsjFourier.getdSjdgamma(j);
        sjCoefs[7][j] = cjsjFourier.getSjLambda(j);

        //In the special case j == 1 there are some additional terms to be added
        if (j == 1) {
            //Additional terms for C<sup>j</sup> coefficients
            cjCoefs[0][j] += - owner->h * owner->U;
            cjCoefs[1][j] += - owner->h * dU[0];
            cjCoefs[2][j] += - owner->h * dU[1];
            cjCoefs[3][j] += -(owner->h * dU[2] + owner->U + cjsjFourier.getC0Lambda());
            cjCoefs[4][j] += - owner->h * dU[3];
            cjCoefs[5][j] += - owner->h * dU[4];
            cjCoefs[6][j] += - owner->h * dU[5];

            //Additional terms for S<sup>j</sup> coefficients
            sjCoefs[0][j] += owner->k * owner->U;
            sjCoefs[1][j] += owner->k * dU[0];
            sjCoefs[2][j] += owner->k * dU[1] + owner->U + cjsjFourier.getC0Lambda();
            sjCoefs[3][j] += owner->k * dU[2];
            sjCoefs[4][j] += owner->k * dU[3];
            sjCoefs[5][j] += owner->k * dU[4];
            sjCoefs[6][j] += owner->k * dU[5];
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// ThirdBodyShortPeriodicCoefficients ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::ThirdBodyShortPeriodicCoefficients(
        const int jMax, const int interpolationPoints, const int maxFreqF,
        const std::string bodyName, const TimeSpanMap<Slot> slotz) :
    jMax(jMax), interpolationPoints(interpolationPoints), maxFreqF(maxFreqF),
    prefix("DSST-3rd-body-" + bodyName + "-"), slotz(slotz) { }


DSSTThirdBody::Slot DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::createSlot(
        const std::vector<SpacecraftState> meanStates)
{
    const Slot         slot  = Slot(jMax, interpolationPoints);
    const AbsoluteDate first = (*meanStates.begin()).getDate();
    const AbsoluteDate last  = (*meanStates.end()).getDate();
    if (first.compareTo(last) <= 0) {
        slotz.addValidAfter(slot, first);
    } else {
        slotz.addValidBefore(slot, first);
    }
    return slot;
}


v_double DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::value(const SpacecraftState spacecraftState)
{
    // select the coefficients slot
    Slot slot = slotz.get(spacecraftState.getDate());

    // the current eccentric longitude
    const double F = orbital_element_conversions::getEccentricLongitude( spacecraftState.equinoctialComponents );

    //initialize the short periodic contribution with the corresponding C⁰ coeficient
    v_double shortPeriodic = slot.cij[0].value(spacecraftState.getDate());

    // Add the cos and sin dependent terms
    for (int j = 1; j <= maxFreqF; j++) {
        //compute cos and sin
        const double cosjF = std::cos(j * F);
        const double sinjF = std::sin(j * F);

        const v_double c = slot.cij[j].value(spacecraftState.getDate());
        const v_double s = slot.sij[j].value(spacecraftState.getDate());
        for (int i = 0; i < 6; i++) {
            shortPeriodic[i] += c[i] * cosjF + s[i] * sinjF;
        }
    }

    return shortPeriodic;
}

std::map<std::string, v_double> DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::getCoefficients(
        const AbsoluteDate date, const std::set<std::string> selected )
{
    // select the coefficients slot
    Slot slot = slotz.get(date);

    std::map<std::string, v_double> coefficients; // FIXME (2 * maxFreqF + 1)
    storeIfSelected(coefficients, selected, slot.cij[0].value(date), "c", {0});
    for (int j = 1; j <= maxFreqF; j++) {
        storeIfSelected(coefficients, selected, slot.cij[j].value(date), "c", {j});
        storeIfSelected(coefficients, selected, slot.sij[j].value(date), "s", {j});
    }
    return coefficients;
}


void DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::storeIfSelected(
        std::map<std::string, v_double> &map, const std::set<std::string> selected,
        const v_double value, const std::string id, const std::vector<int> indices )
{
    std::stringstream keyBuilder;
    keyBuilder << getCoefficientsKeyPrefix() << id;
    for (int index : indices) {
        keyBuilder << "[" << index << "]";
    }
    const std::string key = keyBuilder.str();
    if ( selected.size() == 0 || selected.find(key) != selected.end() ) {
        map[key] = value;
    }
}


DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::DataTransferObject
DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::writeReplace()
{
    // slotz transitions
    std::set< const Transition<Slot>, ChronologicalComparator > transitions = slotz.getTransitions();
    std::vector< AbsoluteDate >              transitionDates( transitions.size() );
    std::vector< Slot >                      allSlots( transitions.size() + 1 );
    unsigned int i = 0;
    for (const Transition<Slot> transition : transitions) {
        if (i == 0) {
            // slot before the first transition
            allSlots[i] = transition.getBefore();
        }
        if (i < transitionDates.size()) {
            transitionDates[i] = transition.getDate();
            allSlots[++i]      = transition.getAfter();
        }
    }

    return DataTransferObject(jMax, interpolationPoints, maxFreqF, prefix, transitionDates, allSlots);
}

DSSTThirdBody::ThirdBodyShortPeriodicCoefficients
DSSTThirdBody::ThirdBodyShortPeriodicCoefficients::DataTransferObject::readResolve()
{
    TimeSpanMap<Slot> slotz = TimeSpanMap<Slot>(allSlotz[0]);
    for (unsigned int i = 0; i < transitionDates.size(); ++i) {
        slotz.addValidAfter(allSlotz[i + 1], transitionDates[i]);
    }

    return ThirdBodyShortPeriodicCoefficients(jMax, interpolationPoints, maxFreqF, prefix, slotz);
}



////////////////////////////////////////////////////////////////////
/////////////////////////////// Slot ///////////////////////////////
////////////////////////////////////////////////////////////////////

DSSTThirdBody::Slot::Slot(const int jMax, const int interpolationPoints)
{
    // allocate the coefficients arrays
    cij = std::vector< ShortPeriodicsInterpolatedCoefficient >( jMax + 1 );
    sij = std::vector< ShortPeriodicsInterpolatedCoefficient >( jMax + 1 );
    for (int j = 0; j <= jMax; j++) {
        cij[j] = ShortPeriodicsInterpolatedCoefficient(interpolationPoints);
        sij[j] = ShortPeriodicsInterpolatedCoefficient(interpolationPoints);
    }
}


} // namespace dsst

} // namespace propagators

} // namespace tudat
