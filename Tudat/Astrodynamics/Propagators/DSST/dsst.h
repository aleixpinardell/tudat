#ifndef TUDAT_DSST_H
#define TUDAT_DSST_H


#include <cmath>
#include <map>
#include <vector>
#include <set>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


typedef Eigen::Matrix< double, 3, 1 > Vector3;
typedef Eigen::Matrix< double, 6, 1 > Vector6;


typedef std::vector< double > v_double;
typedef std::vector< std::vector< double > > vv_double;

static v_double new_v_double( int m )
{
    return v_double( m );
}

static vv_double new_vv_double( int m, int n )
{
    return vv_double( m, v_double( n ) );
}


typedef std::vector< int > v_int;
typedef std::vector< std::vector< int > > vv_int;
typedef std::vector< std::vector< std::vector< int > > > vvv_int;

static v_int new_v_int( int m )
{
    return v_int( m );
}

static vv_int new_vv_int( int m, int n )
{
    return vv_int( m, v_int( n ) );
}

static vvv_int new_vvv_int( int m, int n, int o )
{
    return vvv_int( m, new_vv_int( n, o ) );
}


static double factorial( double x ) { return tgamma( x + 1 ); }


class Equatable
{
public:
    virtual bool operator == ( const Equatable &other ) const;

    virtual bool operator != ( const Equatable &other ) const {
        return ! (*this == other);
    }
};


} // namespace dsst

} // namespace propagators

} // namespace tudat


#endif // TUDAT_DSST_H
