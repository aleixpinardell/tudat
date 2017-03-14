/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GAUSSIAN_INTEGRATOR_H
#define TUDAT_GAUSSIAN_INTEGRATOR_H

#include <vector>
#include <map>

#include <boost/function.hpp>
// #include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{

namespace numerical_quadrature
{


//! Trapezoid numerical quadrature wrapper class.
/*!
 *  Numerical method that uses the trapezoid method to compute definite integrals of a dataset.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class GaussianQuadrature : public NumericalQuadrature< IndependentVariableType , DependentVariableType >
{
public:

    //! Constructor.
    /*!
     * Constructor
     *
     */
    GaussianQuadrature( boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                        IndependentVariableType lowerLimit, IndependentVariableType upperLimit,
                        const unsigned int numberOfAbscissae )
    {
        /*        resetData( integrand, lowerLimit, upperLimit, numberOfAbscissae );
    }


    void resetData( boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                    IndependentVariableType lowerLimit, IndependentVariableType upperLimit,
                    const unsigned int numberOfAbscissae )
    {
*/
        this->lowerLimit = lowerLimit;
        this->upperLimit = upperLimit;

        // Determine the values of the auxiliary independent variable (abscissae)
        const IndependentVariableArray abscissae = getAbscissae( numberOfAbscissae );

        // Determine the values of the weight factors
        const IndependentVariableArray weightFactors = getWeightFactors( numberOfAbscissae );

        // Change of variable -> from range [-1, 1] to range [lowerLimit, upperLimit]
        const IndependentVariableArray independentVariables =
                0.5 * ( ( upperLimit - lowerLimit ) * abscissae + upperLimit + lowerLimit );

        // Determine the value of the dependent variable
        weighedDependentVariables.resizeLike( independentVariables );
        for ( unsigned int i = 0; i < numberOfAbscissae; i++ ) {
            weighedDependentVariables( i ) = weightFactors( i ) * integrand( independentVariables( i ) );
        }

        performQuadrature( );
    }

    //! Function to return computed value of the quadrature.
    /*!
     *  Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     *  \return Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     */
    DependentVariableType getQuadrature( )
    {
        return quadratureResult;
    }


    typedef Eigen::Array<   DependentVariableType, 1, Eigen::Dynamic >   DependentVariableArray;
    typedef Eigen::Array< IndependentVariableType, 1, Eigen::Dynamic > IndependentVariableArray;

    IndependentVariableArray getAbscissae( const unsigned int n )
    {
        IndependentVariableArray abscissae( n );

        // Include abscissa 0.0 if n is odd
        unsigned int i = 0;
        if ( n % 2 == 1 ) {
            abscissae.col( i++ ) = 0.0;
        }

        // Include ± abscissae
        IndependentVariableArray uniqueAbscissae = getUniqueAbscissae( n );
        for ( unsigned int j = 0; j < uniqueAbscissae.size(); j++ ) {
            abscissae.col( i++ ) = -uniqueAbscissae[ j ];
            abscissae.col( i++ ) =  uniqueAbscissae[ j ];
        }

        return abscissae;
    }

    IndependentVariableArray getWeightFactors( const unsigned int n )
    {
        IndependentVariableArray weightFactors( n );

        IndependentVariableArray uniqueWeightFactors = getUniqueWeightFactors( n );

        // Include non-repeated weight factor if n is odd
        unsigned int i = 0;
        unsigned int j = 0;
        if ( n % 2 == 1 ) {
            weightFactors.col( i++ ) = uniqueWeightFactors[ j++ ];
        }

        // Include repeated weight factors
        for ( ; j < uniqueWeightFactors.size(); j++ ) {
            weightFactors.col( i++ ) = uniqueWeightFactors[ j ];
            weightFactors.col( i++ ) = uniqueWeightFactors[ j ];
        }

        return weightFactors;
    }


protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     * Function that is called to perform the numerical quadrature. Sets the result in the quadratureResult_ local
     * variable.
     */
    void performQuadrature( )
    {
        quadratureResult = 0.5 * ( upperLimit - lowerLimit ) * weighedDependentVariables.sum();
    }

private:

    //!
    std::map< unsigned int, IndependentVariableArray > uniqueAbscissae;

    //!
    std::map< unsigned int, IndependentVariableArray > uniqueWeightFactors;

    //! Lower limit for the integral.
    IndependentVariableType lowerLimit;

    //! Upper limit for the integral.
    IndependentVariableType upperLimit;

    //! Dependent variables times the weight factors.
    DependentVariableArray weighedDependentVariables;

    //! Computed value of the quadrature, as computed by last call to performQuadrature.
    DependentVariableType quadratureResult;


    IndependentVariableArray getUniqueAbscissae( const unsigned int n )
    {
        if ( uniqueAbscissae.count( n ) == 0 )
        {
            readAbscissae();
        }
        return uniqueAbscissae.at( n );
    }

    IndependentVariableArray getUniqueWeightFactors( const unsigned int n )
    {
        if ( uniqueWeightFactors.count( n ) == 0 )
        {
            readWeightFactors();
        }
        return uniqueWeightFactors.at( n );
    }


    // Transform from map of std::vector to map of Eigen::Array
    std::map< unsigned int, IndependentVariableArray > vector2eigen(
            std::map< unsigned int, std::vector< IndependentVariableType > > map )
    {
        std::map< unsigned int, IndependentVariableArray > eigenMap;
        for ( auto ent: map )
        {
            IndependentVariableArray array( ent.second.size() );
            for ( unsigned int i = 0; i < array.cols(); i++ )
            {
                array.col( i ) = ent.second.at( i );
            }
            eigenMap[ ent.first ] = array;
        }
        return eigenMap;
    }

    void readAbscissae()
    {
        auto abscissae = input_output::readMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianAbscissae.txt" );
        uniqueAbscissae = vector2eigen( abscissae );
    }

    void readWeightFactors()
    {
        auto weightFactors = input_output::readMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianWeightFactors.txt" );
        uniqueWeightFactors = vector2eigen( weightFactors );
    }


};

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_GAUSSIAN_INTEGRATOR_H
