#ifndef TUDAT_DSST_DERIVATIVESTRUCTURE_H
#define TUDAT_DSST_DERIVATIVESTRUCTURE_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"
#include "DSCompiler.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

class DerivativeStructure {
private:
    DSCompiler compiler;

    Vectord data;


public:

    DerivativeStructure(const DerivativeStructure &ds) : compiler(ds.compiler), data(ds.data) { }

    DerivativeStructure(const DSCompiler compiler) : compiler(compiler), data(compiler.getSize()) { }

    DerivativeStructure(const int parameters, const int order) :
        DerivativeStructure( DSCompiler::getCompiler( parameters, order ) ) { }

    DerivativeStructure(const int parameters, const int order, const double value) :
        DerivativeStructure(parameters, order) {
        this->data[0] = value;
    }

    DerivativeStructure(const int parameters, const int order, const int index, const double value) :
        DerivativeStructure(parameters, order, value)
    {
        if (index >= parameters) {
            throw std::runtime_error( "Number too large when constructing DerivativeStructure" );
        }

        if (order > 0) {
            // the derivative of the variable with respect to itself is 1.
            data[DSCompiler::getCompiler(index, order).getSize()] = 1.0;
        }

    }


    /** Get the number of free parameters.
     * @return number of free parameters
     */
    int getFreeParameters() const {
        return compiler.getFreeParameters();
    }

    /** Get the derivation order.
     * @return derivation order
     */
    int getOrder() const {
        return compiler.getOrder();
    }


    /** {@inheritDoc}
     */
    /*
    DerivativeStructure multiply(const double a) {
        DerivativeStructure ds(*this);
        for (unsigned int i = 0; i < ds.data.size(); ++i) {
            ds.data[i] *= a;
        }
        return ds;
    }
    */

    /** {@inheritDoc} */
    /*
    DerivativeStructure multiply(const int n) {
        return multiply((double) n);
    }
    */

    DerivativeStructure add(const double a) {
        DerivativeStructure ds = DerivativeStructure(*this);
        ds.data[0] += a;
        return ds;
    }

    DerivativeStructure multiply(const DerivativeStructure a)
    {
        // FIXME compiler.checkCompatibility(a.compiler);
        DerivativeStructure result = DerivativeStructure(compiler);
        compiler.multiply(data, 0, a.data, 0, result.data, 0);
        return result;
    }



    /** Get the value part of the derivative structure.
      * @return value part of the derivative structure
      * @see #getPartialDerivative(int...)
      */
     double getValue() const {
         return data[0];
     }

     /** Get a partial derivative.
      * @param orders derivation orders with respect to each variable (if all orders are 0,
      * the value is returned)
      * @return partial derivative
      * @see #getValue()
      * @exception MathIllegalArgumentException if the numbers of variables does not
      * match the instance
      * @exception MathIllegalArgumentException if sum of derivation orders is larger
      * than the instance limits
      */
     double getPartialDerivative(const Vectori orders) const
     {
         return data[compiler.getPartialDerivativeIndex(orders)];
     }



    void operator = (DerivativeStructure ds)
    {
        compiler = ds.compiler;
        data = ds.data;
    }

};




} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_DERIVATIVESTRUCTURE_H
