#ifndef TUDAT_DSST_CONSERVATIVETHIRDBODYPERTURBED_H
#define TUDAT_DSST_CONSERVATIVETHIRDBODYPERTURBED_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "thirdBodyPertubed.h"
#include "conservative.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

//! Abstract class for perturbations caused by third bodies that can be modelled as a potential
class ConservativeThirdBodyPerturbedForceModel :
    public ThirdBodyPerturbedForceModel, public ConservativeForceModel
{
protected:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     * \param thirdBody The third body exerting the acceleration.
     */
    ConservativeThirdBodyPerturbedForceModel( AuxiliaryElements &auxiliaryElements, Body &thirdBody ) :
        ForceModel(                   auxiliaryElements ),
        ThirdBodyPerturbedForceModel( auxiliaryElements, thirdBody ),
        ConservativeForceModel(       auxiliaryElements ) { }


private:

    void setUp() {
        ConservativeForceModel::setUp();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    void updateMembers( ) {
    ThirdBodyPerturbedForceModel::updateMembers();
    ConservativeForceModel::updateMembers();
	Ufactor = computeUfactor();
    }

    //! Set the values of N and S for the series expansion of the disturbing potential.
    void determineTruncationValues( );

    //! Function to compute the U derivatives for the current state [ Eq. 3.2-(2) ]
    Eigen::Vector6d computeMeanDisturbingFunctionPartialDerivatives( );

    //! Get the short period terms for the current auxiliary elements [ Sections 2.5.2 and 4.2 ]
    Eigen::Vector6d computeShortPeriodTerms( );


    // Terms used in several parts of the class

    //! a / R3 up to power maxAR3Pow
    Vectord aR3pow;


    // Constants used to determine the suitable N and S for the serie expansion

    //! Max power for summation
    virtual unsigned int MAX_N() { return 22; }

    //! Truncation tolerance for big, eccentric  orbits
    virtual double BIG_TRUNCATION_TOLERANCE() { return 1.e-1; }

    //! Truncation tolerance for small orbits
    virtual double SMALL_TRUNCATION_TOLERANCE() { return 1.9e-6; }

    //! Number of points for interpolation
    // virtual unsigned int INTERPOLATION_POINTS() { return 3; }

    //! Maximum power for eccentricity used in short periodic computation
    // virtual unsigned int MAX_ECCPOWER_SP() { return 4; }


    // Shared by third body central gravity and solar radiation pressure without eclipses

    //! Factor for U outside the serie expansion [ see Eq. 3.2-(1) ]
    double Ufactor;

    //! Function that updates `Ufactor`
    virtual double computeUfactor( ) = 0;

    //! Minimum value for n in the mean disturbing potential serie expansion
    virtual unsigned int nMin( unsigned int s ) = 0;

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_CONSERVATIVETHIRDBODYPERTURBED_H
