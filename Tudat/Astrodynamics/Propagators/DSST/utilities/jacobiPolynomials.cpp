/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "jacobiPolynomials.h"
#include "polynomialUtil.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{

DerivativeStructure JacobiPolynomials::getValue(const int l, const int v, const int w, const DerivativeStructure gamma)
{

    std::vector<PolynomialFunction> polyList;
    const JacobiKey key = JacobiKey(v, w);

    // Check the existence of the corresponding key in the map.
    if (MAP.count(key) == 0) {
        MAP[key] = std::vector<PolynomialFunction>();
    }

    polyList = MAP[key];

    // If the l-th degree polynomial has not been computed yet, the polynomials
    // up to this degree are computed.
    for (int degree = polyList.size(); degree <= l; degree++) {
        polyList.insert( polyList.begin() + degree, PolynomialsUtils::createJacobiPolynomial(degree, v, w));
    }

    // compute value and derivative
    return polyList[l].value(gamma);
}

} // namespace dsst

} // namespace propagators

} // namespace tudat
