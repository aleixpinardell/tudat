#ifndef TUDAT_PROPAGATORS_DSST_BODY_H
#define TUDAT_PROPAGATORS_DSST_BODY_H


#include <cmath>
#include <map>
#include <vector>
#include <set>

#include "Eigen/Eigen"

#include "boost/boost/function.hpp"

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


//! Class for a body with a fixed name, and variable or constant position, mass, area, CD and CR
class Body
{
public:

    //! Simple constructor
    Body( const std::string &name ) : name( name ) { }

    //! Body's name
    const std::string name;


    //! Set body's position function
    void definePosition( boost::function< Eigen::Vector3d() > function ) {
        position = function;
    }

    //! Get body's position [ m, m, m ]
    Eigen::Vector3d getPosition() const {
        return position();
    }

    //! Get position relative to another body [ m, m, m ]
    Eigen::Vector3d getPositionFrom( boost::shared_ptr< Body > body ) const {
        return getPosition() - body->getPosition();
    }

    //! Get distance between to another body [ m ]
    double getDistanceTo( boost::shared_ptr< Body > body ) const {
        return getPositionFrom( body ).norm();
    }


    //! Set body's mass function [ kg ]
    void defineMass( boost::function< double() > function ) {
        mass = function;
    }

    //! Get body's mass [ kg ]
    virtual double getMass() const {
        return mass();
    }


    //! Set body's surface area function [ m^2 ]
    void defineSurfaceArea( boost::function< double() > function ) {
        surfaceArea = function;
    }

    //! Get body's surface area [ m^2 ]
    virtual double getSurfaceArea() const {
        return surfaceArea();
    }


    //! Set body's cross-sectional area function [ m^2 ]
    void defineCrossSectionalArea( boost::function< double() > function ) {
        crossSectionalArea = function;
    }

    //! Get body's cross-sectional area [ m^2 ]
    virtual double getCrossSectionalArea() const {
        return crossSectionalArea();
    }


    //! Set body's drag coefficient function
    void defineDragCoefficient( boost::function< double() > function ) {
        dragCoefficient = function;
    }

    //! Get body's drag coefficient
    virtual double getDragCoefficient() const {
        return dragCoefficient();
    }


    //! Set body's radiation pressure coefficient function
    void defineRadiationPressureCoefficient( boost::function< double() > function ) {
        radiationPressureCoefficient = function;
    }

    //! Get body's radiation pressure coefficient
    virtual double getRadiationPressureCoefficient() const {
        return radiationPressureCoefficient();
    }


private:

    //! Body's position in the reference system used in the propagation [ m, m, m ]
    boost::function< Eigen::Vector3d() > position;

    //! Body's mass [ kg ]
    boost::function< double() > mass;

    //! Body's surface area [ m^2 ]
    boost::function< double() > surfaceArea;

    //! Body's cross-sectional area [ m^2 ]
    boost::function< double() > crossSectionalArea;

    //! Body's drag coefficient
    boost::function< double() > dragCoefficient;

    //! Body's radiation pressure coefficient
    boost::function< double() > radiationPressureCoefficient;

};


//! Class for a body whose mass is defined through its gravitational parameter ( µ known more accurately than mass )
class CelestialBody : public virtual Body
{
public:

    //! Inherit constructor
    using Body::Body;


    //! Set body's gravitational parameter function [ m^3/s^2 ]
    void defineGravitationalParameter( boost::function< double() > function ) {
        gravitationalParameter = function;
    }

    //! Get body's gravitational parameter [ m^3/s^2 ]
    double getGravitationalParameter() const {
        return gravitationalParameter();
    }


    //! Set body's rotational velocity function [ rad/s, rad/s, rad/s ]
    void defineRotationalVelocity( boost::function< Eigen::Vector3d() > function ) {
        rotationalVelocity = function;
    }

    //! Get body's rotational velocity [ rad/s, rad/s, rad/s ]
    Eigen::Vector3d getRotationalVelocity() const {
        return rotationalVelocity();
    }


    //! Get body's mass [ kg ]
    virtual double getMass() const {
        return getGravitationalParameter() / physical_constants::GRAVITATIONAL_CONSTANT;
    }


private:

    //! Body's gravitational parameter [ m^3/s^2 ]
    boost::function< double() > gravitationalParameter;

    //! Body's rotational velocity [ rad/s, rad/s, rad/s ]
    boost::function< Eigen::Vector3d() > rotationalVelocity;


};


//! Class for a spherical body
class SphericalBody : public virtual Body
{
public:

    //! Inherit constructor
    using Body::Body;


    //! Set body's radius function [ m ]
    void defineRadius( boost::function< double() > function ) {
        radius = function;
    }

    //! Get body's radius [ m ]
    double getRadius() const {
        return radius();
    }


    //! Get body's surface area [ m^2 ]
    virtual double getSurfaceArea() const {
        return 4 * getCrossSectionalArea();
    }


    //! Get body's cross-sectional area [ m^2 ]
    virtual double getCrossSectionalArea() const {
        return mathematical_constants::PI * std::pow( getRadius(), 2 );
    }


private:

    //! Body's radius [ m ]
    boost::function< double() > radius;


};


//! Class for a radiative, spherical, celestial body
class SphericalCelestialBody : public virtual CelestialBody, public virtual SphericalBody
{
public:

    //! Constructor
    SphericalCelestialBody( const std::string &name ) : Body( name ), CelestialBody( name ), SphericalBody( name ) { }

};





/*
//! Class for a radiative, spherical, celestial body
class RadiativeSphericalCelestialBody : public virtual SphericalCelestialBody, public virtual RadiativeBody
{
public:

    //! Inherit constructor
    using Body::Body;

};

typedef RadiativeSphericalCelestialBody Star;
*/


} // namespace dsst

} // namespace propagators

} // namespace tudat


#endif // TUDAT_PROPAGATORS_DSST_BODY_H
