#ifndef TUDAT_DSST_SHORTPERIODTERMS_H
#define TUDAT_DSST_SHORTPERIODTERMS_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/absoluteDate.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/spacecraftState.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

/** Additive short period terms contributing to the mean to osculating orbit mapping.
 * <p>
 * Each instance contains a set of several terms that are computed together.
 * </p>
 * @see DSSTForceModel
 */
class ShortPeriodTerms {

    /** Evaluate the contributions of the short period terms.
     * @param meanOrbit mean orbit to which the short period contribution applies
     * @return short period terms contributions
     */
    virtual v_double value(SpacecraftState spacecraftState);

    /** Get the prefix for short period coefficients keys.
     * <p>
     * This prefix is used to identify the coefficients of the
     * current force model from the coefficients pertaining to
     * other force models. All the keys in the map returned by
     * {@link #getCoefficients(AbsoluteDate, Set)}
     * start with this prefix, which must be unique among all
     * providers.
     * </p>
     * @return the prefix for short periodic coefficients keys
     * @see #getCoefficients(AbsoluteDate, Set)
     */
    virtual std::string getCoefficientsKeyPrefix();

    /** Computes the coefficients involved in the contributions.
     * <p>
     * This method is intended mainly for validation purposes. Its output
     * is highly dependent on the implementation details in each force model
     * and may change from version to version. It is <em>not</em> recommended
     * to use it for any operational purposes.
     * </p>
     * @param date current date
     * @param selected set of coefficients that should be put in the map
     * (empty set means all coefficients are selected)
     * @return the selected coefficients of the short periodic variations,
     * in a map where all keys start with {@link #getCoefficientsKeyPrefix()}
     */
    virtual std::map<std::string, v_double> getCoefficients(AbsoluteDate date, std::set<std::string> selected);

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_SHORTPERIODTERMS_H
