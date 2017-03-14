#ifndef TUDAT_DSST_PROPAGATORS_FORCEMODELS_CONSERVATIVESOLARRADIATIONPRESSURE_H
#define TUDAT_DSST_PROPAGATORS_FORCEMODELS_CONSERVATIVESOLARRADIATIONPRESSURE_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

// #include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "abstract/conservativeThirdBodyPertubed.h"

// #include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

namespace force_models
{


//! Final class for the contribution of a third body's central gravity
class ConservativeRadiationPressure final : public ConservativeThirdBodyPerturbed
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     * \param thirdBody The third body exerting the acceleration.
     */
    ConservativeRadiationPressure(
            AuxiliaryElements &auxiliaryElements, boost::shared_ptr< CelestialBody > thirdBody,
            boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface
            ) :
        ForceModel( auxiliaryElements ),
        ConservativeThirdBodyPerturbed( auxiliaryElements, thirdBody ),
        radiationPressureInterface( radiationPressureInterface ), mass( auxiliaryElements.mass ) { }


private:

    //! Pointer to the radiation pressure interface
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface;

    //! Propagated body's mass [ kg ]
    double &mass;

    double computeUfactor() {
        using namespace mathematical_constants;
        using namespace physical_constants;

        // radiationPressureInterface->updateInterface();
        /// FIXME: Who calls this? Probably not needed for this case...
        /// Why? Because Ufactor has a constant value, does not depend on the distance (r - r3) !!
        /// And eclipses ignored, so no need to check for ocultations...

        boost::function< double( ) > power = radiationPressureInterface->getSourcePowerFunction();
        const double area = radiationPressureInterface->getArea();
        const double CR = radiationPressureInterface->getRadiationPressureCoefficient();

        // Eq. 3.5-(7)
        const double T = 0.5 * CR * area / mass * power() / ( 4 * PI * SPEED_OF_LIGHT );
        // Eq. 3.5-(10)
        return - T / R3;
    }

    //! Minimum value for n in the expansion
    unsigned int nMin( unsigned int s ) {
        return std::max( 1, (int) s );
    }

};


} // namespace force_models

} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_PROPAGATORS_FORCEMODELS_CONSERVATIVERADIATIONPRESSURE_H
