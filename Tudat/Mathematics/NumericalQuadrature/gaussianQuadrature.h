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


//! Gaussian numerical quadrature wrapper class.
/*!
 * Numerical method that uses the Gaussian abscissae and weight factors to compute definite integrals of a function.
 * The Gaussian abscissae and weight factors are not calculated, but read from text files. The number of abscissae (or
 * weight factors) has to be at least n = 2. The current text files contain tabulated values up to n = 64.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class GaussianQuadrature : public NumericalQuadrature< IndependentVariableType , DependentVariableType >
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param integrand Function to be integrated numerically.
     * \param lowerLimit Lower limit for the integral.
     * \param upperLimit Upper limit for the integral.
     * \param numberOfAbscissae Number of abscissae at which the integrand will be evaluated. Must be an integer value
     * between 2 and 64.
     */
    GaussianQuadrature( boost::function< DependentVariableType( IndependentVariableType ) > integrand,
                        IndependentVariableType lowerLimit, IndependentVariableType upperLimit,
                        const unsigned int numberOfAbscissae )
    {
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
        weighedIntegrands.resizeLike( independentVariables );
        for ( unsigned int i = 0; i < numberOfAbscissae; i++ ) {
            weighedIntegrands( i ) = weightFactors( i ) * integrand( independentVariables( i ) );
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

    //! Get all the abscissae (i.e. n abscissae for nth order) from uniqueAbscissae
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

    //! Get all the weight factors (i.e. n weight factors for nth order) from uniqueWeightFactors
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
        quadratureResult = 0.5 * ( upperLimit - lowerLimit ) * weighedIntegrands.sum();
    }

private:

    //! Map containing the abscissae read from the text file (currently up to `n = 64`).
    //! The following relation holds: `size( uniqueAbscissae[n] ) = floor( n / 2 )`
    //! For the actual abscissae, the following must hold: `size( abscissae[n] ) = n`
    //! The actual abscissae are generated from `uniqueAbscissae` by `getAbscissae()`
    std::map< unsigned int, IndependentVariableArray > uniqueAbscissae;

    //! Map containing the weight factors read from the text file (currently up to `n = 64`).
    //! The following relation holds: `size( uniqueWeightFactors[n] ) = ceil( n / 2 )`
    //! For the actual weight factors, the following must hold: `size( uniqueWeightFactors[n] ) = n`
    //! The actual weight factors are generated from `uniqueWeightFactors` by `getWeightFactors()`
    std::map< unsigned int, IndependentVariableArray > uniqueWeightFactors;

    //! Lower limit for the integral.
    IndependentVariableType lowerLimit;

    //! Upper limit for the integral.
    IndependentVariableType upperLimit;

    //! Integrands at the abscissae, times the respective weight factors.
    DependentVariableArray weighedIntegrands;

    //! Computed value of the quadrature, as computed by last call to performQuadrature.
    DependentVariableType quadratureResult;

    //! Get the unique abscissae for a specified order `n`.
    /*!
     * \param n The number of abscissae or weight factors.
     * \return `uniqueAbscissae[n]`, after reading the text file with the tabulated abscissae if necessary.
     */
    IndependentVariableArray getUniqueAbscissae( const unsigned int n )
    {
        if ( uniqueAbscissae.count( n ) == 0 )
        {
            readAbscissae();
        }
        return uniqueAbscissae.at( n );
    }

    //! Get the unique weight factors for a specified order `n`.
    /*!
     * \param n The number of abscissae or weight factors.
     * \return `uniqueWeightFactors[n]`, after reading the text file with the tabulated abscissae if necessary.
     */
    IndependentVariableArray getUniqueWeightFactors( const unsigned int n )
    {
        if ( uniqueWeightFactors.count( n ) == 0 )
        {
            readWeightFactors();
        }
        return uniqueWeightFactors.at( n );
    }


    //! Transform from map of std::vector (output of text file reader) to map of Eigen::Array
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

    //! Read Gaussian abscissae from text file
    void readAbscissae()
    {
        auto abscissae = input_output::readMapFromFile< unsigned int, IndependentVariableType >(
                    input_output::getTudatRootPath( ) + "/Mathematics/NumericalQuadrature/gaussianAbscissae.txt" );
        uniqueAbscissae = vector2eigen( abscissae );
    }

    //! Read Gaussian weight factors from text file
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
