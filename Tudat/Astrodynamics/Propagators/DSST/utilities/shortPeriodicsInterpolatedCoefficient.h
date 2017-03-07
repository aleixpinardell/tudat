#ifndef TUDAT_DSST_SHORTPERIODICSINTERPOLATEDCOEFFICIENT_H
#define TUDAT_DSST_SHORTPERIODICSINTERPOLATEDCOEFFICIENT_H

// #include "Tudat/Mathematics/Interpolators/hermiteCubicSplineInterpolator.h"

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "absoluteDate.h"
#include "hermiteInterpolator.h"
#include "orbit.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

/** Interpolated short periodics coefficients.
 * <p>
 * Representation of a coefficient that need to be interpolated over time.
 * </p><p>
 * The short periodics coefficients can be interpolated for faster computation.
 * This class stores computed values of the coefficients through the method
 * {@link #addGridPoint} and gives an interpolated result through the method
 * {@link #value}.
 * </p>
 *
 */
class ShortPeriodicsInterpolatedCoefficient : public Equatable {

private:
    /** Values of the already computed coefficients.*/
    NestedVectord values;

    /** Grid points.*/
    std::vector<AbsoluteDate> abscissae;

    /** Number of points used in the interpolation.*/
    int interpolationPoints;

    /** Index of the latest closest neighbor.*/
    int latestClosestNeighbor = 0;


    /** Find the closest point from a specific date amongst the available points.
     * @param date date of interest
     * @return index of the closest abscissa from the date of interest
     */
    int getClosestNeighbor(const AbsoluteDate date) {
        //the starting point is the latest result of a call to this method.
        //Indeed, as this class is meant to be called during an integration process
        //with an input date evolving often continuously in time, there is a high
        //probability that the result will be the same as for last call of
        //this method.
        int closestNeighbor = latestClosestNeighbor;

        //case where the date is before the available points
        if (date.compareTo(abscissae[0]) <= 0) {
            closestNeighbor = 0;
        }
        //case where the date is after the available points
        else if (date.compareTo(abscissae[abscissae.size() - 1]) >= 0) {
            closestNeighbor = abscissae.size() - 1;
        }
        //general case: one is looking for the two consecutives entries that surround the input date
        //then one choose the closest one
        else {
            int lowerBorder = latestClosestNeighbor;
            int upperBorder = latestClosestNeighbor;

            const int searchDirection = date.compareTo(abscissae[latestClosestNeighbor]);
            if (searchDirection > 0) {
                upperBorder++;
                while (date.compareTo(abscissae[upperBorder]) > 0) {
                    upperBorder++;
                    lowerBorder++;
                }
            }
            else {
                lowerBorder--;
                while (date.compareTo(abscissae[lowerBorder]) < 0) {
                    upperBorder--;
                    lowerBorder--;
                }
            }

            const double lowerDistance = std::fabs(date.durationFrom(abscissae[lowerBorder]));
            const double upperDistance = std::fabs(date.durationFrom(abscissae[upperBorder]));

            closestNeighbor = (lowerDistance < upperDistance) ? lowerBorder : upperBorder;
        }

        //The result is stored in order to speed up the next call to the function
        //Indeed, it is highly likely that the requested result will be the same
        latestClosestNeighbor = closestNeighbor;
        return closestNeighbor;
    }


    /** Find the closest available points from the specified date.
     * @param date date of interest
     * @return indices corresponding to the closest points on the time scale
     */
    std::vector< int > getNeighborsIndices(const AbsoluteDate date) {
        const int sizeofNeighborhood = std::min(interpolationPoints, (int) abscissae.size());
        std::vector< int > neighborsIndices(sizeofNeighborhood);

        //If the size of the complete sample is less than
        //the desired number of interpolation points,
        //then the entire sample is considered as the neighborhood
        if (interpolationPoints >= (int) abscissae.size()) {
            for (int i = 0; i < sizeofNeighborhood; i++) {
                neighborsIndices[i] = i;
            }
        } else {
            // get indices around closest neighbor
            int inf = getClosestNeighbor(date);
            int sup = inf + 1;

            while (sup - inf < interpolationPoints) {
                if (inf == 0) { //This means that we have reached the earliest date
                    sup++;
                } else if (sup >= (int) abscissae.size()) { //This means that we have reached the latest date
                    inf--;
                } else { //the choice is made between the two next neighbors
                    const double lowerNeighborDistance = std::fabs(abscissae[inf - 1].durationFrom(date));
                    const double upperNeighborDistance = std::fabs(abscissae[sup].durationFrom(date));

                    if (lowerNeighborDistance <= upperNeighborDistance) {
                        inf--;
                    } else {
                        sup++;
                    }
                }
            }

            for (int i = 0; i < interpolationPoints; ++i) {
                neighborsIndices[i] = inf + i;
            }

        }

        return neighborsIndices;
    }


public:
    ShortPeriodicsInterpolatedCoefficient() { }

    /** Simple constructor.
     * @param interpolationPoints number of points used in the interpolation
     */
    ShortPeriodicsInterpolatedCoefficient(const int interpolationPoints) : interpolationPoints(interpolationPoints) {}


    /** Compute the value of the coefficient.
     * @param date date at which the coefficient should be computed
     * @return value of the coefficient
     */
    Vectord value(const AbsoluteDate date) {
        //Get the closest points from the input date
        const std::vector<int> neighbors = getNeighborsIndices(date);

        /*
        //Creation and set up of the interpolator
        v_double independentVariables;
        v_double dependentVariables;
        v_double derivatives;

        for (int i : neighbors) {
            independentVariables.push_back( abscissae[i].durationFrom(date) );
            dependentVariables.push_back( values[0][i] );
            derivatives.push_back( values[1][i] );
            // dependentVariables.push_back( values[i][0] );
            // derivatives.push_back( values[i][1] );
        }

        interpolators::HermiteCubicSplineInterpolatorDouble interpolator(
                    independentVariables, dependentVariables, derivatives );

        //interpolation
        return { interpolator.interpolate( 0.0 ) };
        */

        HermiteInterpolator interpolator;
        for (int i : neighbors) {
            interpolator.addSamplePoint( abscissae[i].durationFrom(date), {values[i]} );
        }

        //interpolation
        return interpolator.value(0.0);

    }


    /** Clear the recorded values from the interpolation grid.
     */
    void clearHistory() {
        abscissae.clear();
        values.clear();
    }

    /** Add a point to the interpolation grid.
     * @param date abscissa of the point
     * @param value value of the element
     */
    void addGridPoint(const AbsoluteDate date, const Vectord value) {
        ptrdiff_t pos = find(abscissae.begin(), abscissae.end(), date) - abscissae.begin();
        //If the grid is empty, the value is directly added to both arrays
        if (abscissae.size() == 0) {
            abscissae.push_back(date);
            values.push_back(value);
        }
        //If the grid already contains this point, only its value is changed
        else if (pos < (int) abscissae.size()) {
            values[pos] = value;
        }
        //If the grid does not contain this point, the position of the point
        //in the grid is computed first
        else {
            const int closestNeighbor = getClosestNeighbor(date);
            const int index = (date.compareTo(abscissae[closestNeighbor]) < 0) ?
                        closestNeighbor : closestNeighbor + 1;
            abscissae.insert( abscissae.begin() + index, date );
            values.insert( values.begin() + index, value );
        }
    }


    bool operator == ( const ShortPeriodicsInterpolatedCoefficient &other ) const {
        return values == other.values;
    }

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_SHORTPERIODICSINTERPOLATEDCOEFFICIENT_H
