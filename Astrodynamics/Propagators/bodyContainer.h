/*! \file bodyContainer.h
 *    Header file that defines the class for all bodies that will be propagated
 *    by the Tudat propagators.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 20 September, 2010
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    author              comment
 *      100920    K. Kumar            File created.
 *      100926    K. Kumar            Filename changed, Doxygen comments added.
 *      100927    K. Kumar            Set functions removed, bodyContainer.cpp
 *                                    merged.
 *      100929    J. Melman           Deleted destructor and constructor
 *                                    implementations.
 *                                    stateVectorRangeStart_ ->
 *                                    stateVectorStartIndex_.
 *      100929    K. Kumar            EigenRoutines.h replaced in include
 *                                    statements by linearAlgebra.h.
 */

#ifndef BODYCONTAINER_H
#define BODYCONTAINER_H

// Include statements.
#include <vector>
#include "forceModel.h"
#include "propagator.h"
#include "linearAlgebra.h"

// Forward declarations.
class Propagator;
class NumericalPropagator;

//! Body container class.
/*!
 * Class containing properties of bodies to be propagated by Tudat propagators.
 */
class BodyContainer
{
public:

    // Definition of friendships.
    // Set Propagator class as friend.
    friend class Propagator;

    // Set NumericalPropagator class as friend.
    friend class NumericalPropagator;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    BodyContainer( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~BodyContainer( );

protected:

private:

    //! Starting index of stateVector.
    /*!
     * Starting index of stateVector in assembled state vector.
     */
    unsigned int stateVectorStartIndex_;

    //! Size of state vector.
    /*!
     * Size of state vector in assembled state vector.
     */
    unsigned int sizeOfStateVector_;

    //! Initial state vector.
    /*!
     * Initial state vector.
     */
    VectorXd initialStateVector_;

    //! State vector.
    /*!
     * State vector.
     */
    VectorXd stateVector_;

    //! Final state vector.
    /*!
     * Final state vector.
     */
    VectorXd finalStateVector_;

    //! A map of propagation history.
    /*!
     * A map of propagation history with propagation time taken as key.
     */
    std::map < double, VectorXd > propagationHistory_;

    //! Vector container of pointers to force models.
    /*!
     * Vector container of pointers to force models.
     */
    std::vector < ForceModel* > vectorContainerOfPointersToForceModels_;

    //! Pointer to Body class.
    /*!
     * Pointer to Body class.
     */
    Body* pointerToBody_;

    //! Pointer to Propagator class.
    /*!
     * Pointer to Propagator class.
     */
    Propagator* pointerToPropagator_;
};

#endif // BODYCONTAINER_H

// End of file.