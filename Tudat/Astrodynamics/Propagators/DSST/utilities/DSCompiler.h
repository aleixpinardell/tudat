#ifndef TUDAT_DSST_DSCOMPILER_H
#define TUDAT_DSST_DSCOMPILER_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{

class DSCompiler {
public:
    typedef std::vector< std::vector< std::reference_wrapper < DSCompiler > > >
    DSCompilerReferences;

private:
    static DSCompilerReferences compilers;

    /** Number of free parameters. */
    int parameters;

    /** Derivation order. */
    int order;

    /** Number of partial derivatives (including the single 0 order derivative element). */
    NestedVectori sizes;

    /** Indirection array of the lower derivative elements. */
    Vectori lowerIndirection;

    /** Indirection arrays for multiplication. */
    TripleVectori multIndirection;



    NestedVectori compileSizes(const int parameters, const int order, const DSCompiler valueCompiler) {
        NestedVectori sizes = NestedVectori(parameters + 1,order + 1);
        if (parameters == 0) {
            sizes[0] = Vectori(order + 1, 1);
        } else {
            sizes = valueCompiler.sizes;
            sizes[parameters][0] = 1;
            for (int i = 0; i < order; ++i) {
                sizes[parameters][i + 1] = sizes[parameters][i] + sizes[parameters - 1][i + 1];
            }
        }
        return sizes;
    }


    /** Compile the lower derivatives indirection array.
      * <p>
      * This indirection array contains the indices of all elements
      * except derivatives for last derivation order.
      * </p>
      * @param parameters number of free parameters
      * @param order derivation order
      * @param valueCompiler compiler for the value part
      * @param derivativeCompiler compiler for the derivative part
      * @return lower derivatives indirection array
      */
     Vectori compileLowerIndirection(const int parameters, const int order,
                                   const DSCompiler valueCompiler, const DSCompiler derivativeCompiler) {

         if (parameters == 0 || order <= 1) {
             return Vectori( { 0 } );
         }

         // this is an implementation of definition 6 in Dan Kalman's paper.
         const int vSize = valueCompiler.lowerIndirection.size();
         const int dSize = derivativeCompiler.lowerIndirection.size();
         Vectori lowerIndirection(vSize + dSize);

         for ( int i = 0; i < vSize; i++ )
         {
             lowerIndirection[i] = valueCompiler.lowerIndirection[i];
         }

         for ( int i = 0; i < dSize; ++i ) {
             lowerIndirection[vSize + i] = valueCompiler.getSize() + derivativeCompiler.lowerIndirection[i];
         }

         return lowerIndirection;

     }

    /** Compile the multiplication indirection array.
     * <p>
     * This indirection array contains the indices of all pairs of elements
     * involved when computing a multiplication. This allows a straightforward
     * loop-based multiplication (see {@link #multiply(double[], int, double[], int, double[], int)}).
     * </p>
     * @param parameters number of free parameters
     * @param order derivation order
     * @param valueCompiler compiler for the value part
     * @param derivativeCompiler compiler for the derivative part
     * @param lowerIndirection lower derivatives indirection array
     * @return multiplication indirection array
     */
    TripleVectori compileMultiplicationIndirection(const int parameters, const int order, const DSCompiler valueCompiler,
                                             const DSCompiler derivativeCompiler, const Vectori lowerIndirection) {

        if ((parameters == 0) || (order == 0)) {
            return { { { 1, 0, 0 } } };
        }

        // this is an implementation of definition 3 in Dan Kalman's paper.
        const int vSize = valueCompiler.multIndirection.size();
        const int dSize = derivativeCompiler.multIndirection.size();
        TripleVectori multIndirection = TripleVectori(vSize + dSize, 0, 0);

        for ( int i = 0; i < vSize; i++ )
        {
            multIndirection[i] = valueCompiler.multIndirection[i];
        }

        for ( int i = 0; i < dSize; ++i ) {
            const NestedVectori dRow = derivativeCompiler.multIndirection[i];
            NestedVectori row(dRow.size() * 2);
            for (unsigned int j = 0; j < dRow.size(); ++j) {
                row.push_back( Vectori( { dRow[j][0], lowerIndirection[dRow[j][1]], vSize + dRow[j][2] } ) );
                row.push_back( Vectori( { dRow[j][0], vSize + dRow[j][1], lowerIndirection[dRow[j][2]] } ) );
            }

            // combine terms with similar derivation orders
            NestedVectori combined(row.size());
            for ( int j = 0; j < (int) row.size(); ++j ) {
                Vectori termJ = row[j];
                if (termJ[0] > 0) {
                    for (unsigned int k = j + 1; k < row.size(); ++k) {
                        Vectori termK = row[k];
                        if (termJ[1] == termK[1] && termJ[2] == termK[2]) {
                            // combine termJ and termK
                            termJ[0] += termK[0];
                            // make sure we will skip termK later on in the outer loop
                            termK[0] = 0;
                        }
                    }
                    combined.push_back(termJ);
                }
            }

            multIndirection[vSize + i] = combined;

        }

        return multIndirection;

    }


    /** Get the index of a partial derivative in an array.
     * @param parameters number of free parameters
     * @param order derivation order
     * @param sizes sizes array
     * @param orders derivation orders with respect to each parameter
     * (the length of this array must match the number of parameters)
     * @return index of the partial derivative
     * @exception MathIllegalArgumentException if sum of derivation orders is larger
     * than the instance limits
     */
    static int getPartialDerivativeIndex(
            const int parameters, const int order, const NestedVectori sizes, const Vectori orders)
    {

        // the value is obtained by diving into the recursive Dan Kalman's structure
        // this is theorem 2 of his paper, with recursion replaced by iteration
        int index     = 0;
        int m         = order;
        int ordersSum = 0;
        for (int i = parameters - 1; i >= 0; --i) {

            // derivative order for current free parameter
            int derivativeOrder = orders[i];

            // safety check
            ordersSum += derivativeOrder;
            if (ordersSum > order) {
                throw std::runtime_error( "Number too large when getting partial derivative index" );
            }

            while (derivativeOrder-- > 0) {
                // as long as we differentiate according to current free parameter,
                // we have to skip the value part and dive into the derivative part
                // so we add the size of the value part to the base index
                index += sizes[i][m--];
            }

        }

        return index;

    }


public:
    DSCompiler() { }

    DSCompiler(const int parameters, const int order,
               const DSCompiler valueCompiler, const DSCompiler derivativeCompiler) :
        parameters(parameters), order(order) {
        sizes = compileSizes( parameters, order, valueCompiler );
        lowerIndirection = compileLowerIndirection(parameters, order, valueCompiler, derivativeCompiler);
        multIndirection = compileMultiplicationIndirection(
                    parameters, order, valueCompiler, derivativeCompiler, lowerIndirection);
    }

    static DSCompiler getCompiler(int parameters, int order) {

        // get the cached compilers
        DSCompilerReferences cache = compilers;
        if ( (int)cache.size() > parameters )
        {
            if ( (int)cache[parameters].size() > order )
            {
                return cache[parameters][order];
            }
        }

        // we need to create more compilers
        const int maxParameters = std::max(parameters, (int) cache.size());
        const int maxOrder      = std::max(order,      cache.size() == 0 ? 0 : (int) cache[0].size());
        DSCompiler blankCompiler;
        std::reference_wrapper < DSCompiler > null = std::ref( blankCompiler );
        DSCompilerReferences newCache( maxParameters + 1,
                                       std::vector< std::reference_wrapper < DSCompiler > >(
                                           maxOrder + 1, null ) );

        // preserve the already created compilers
        for (unsigned int i = 0; i < cache.size(); ++i) {
            newCache[i] = cache[i];
        }

        // create the array in increasing diagonal order
        for (int diag = 0; diag <= parameters + order; ++diag) {
            for (int o = std::max(0, diag - parameters); o <= std::min(order, diag); ++o) {
                const int p = diag - o;
                if ( (int)newCache.size() > p ) { continue; }
                if ( (int)newCache[p].size() > o ) { continue; }
                const DSCompiler valueCompiler      = (p == 0) ? null : newCache[p - 1][o];
                const DSCompiler derivativeCompiler = (o == 0) ? null : newCache[p][o - 1];
                DSCompiler newCompiler = DSCompiler(p, o, valueCompiler, derivativeCompiler);
                newCache[p][o] = std::ref(newCompiler);
            }
        }

        // atomically reset the cached compilers array
        compilers = newCache;

        return newCache[parameters][order];

    }


    /** Get the number of free parameters.
     * @return number of free parameters
     */
    int getFreeParameters() const {
        return parameters;
    }

    /** Get the derivation order.
     * @return derivation order
     */
    int getOrder() const {
        return order;
    }

    /** Get the array size required for holding partial derivatives data.
     * <p>
     * This number includes the single 0 order derivative element, which is
     * guaranteed to be stored in the first element of the array.
     * </p>
     * @return array size required for holding partial derivatives data
     */
    int getSize() const {
        return sizes[parameters][order];
    }


    void multiply(const Vectord lhs, const int lhsOffset, const Vectord rhs,
                  const int rhsOffset, Vectord &result, const int resultOffset)
    {
        for (unsigned int i = 0; i < multIndirection.size(); ++i) {
            const NestedVectori mappingI = multIndirection[i];
            double r = 0;
            for (unsigned int j = 0; j < mappingI.size(); ++j) {
                r += mappingI[j][0] *
                     lhs[lhsOffset + mappingI[j][1]] *
                     rhs[rhsOffset + mappingI[j][2]];
            }
            result[resultOffset + i] = r;
        }
    }


    int getPartialDerivativeIndex(const Vectori orders) const
    {
        if ( (int) orders.size() != getFreeParameters() )
        {
            throw std::runtime_error( "orders.size() and getFreeParameters() do not match" );
        }
        return getPartialDerivativeIndex(parameters, order, sizes, orders);
    }

    /*
    void set(DSCompiler compiler)
    {
        parameters = compiler.parameters;
        order = compiler.order;
        sizes = compiler.sizes;
        lowerIndirection = compiler.lowerIndirection;
        multIndirection = compiler.multIndirection;
    }
    */

};


} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_DSCOMPILER_H
