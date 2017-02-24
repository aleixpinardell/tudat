#ifndef TUDAT_DSST_ABSOLUTEDATE_H
#define TUDAT_DSST_ABSOLUTEDATE_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

namespace tudat
{

namespace propagators
{

namespace dsst
{


class AbsoluteDate /*: public TimeStamped*/ {

private:
    /** Reference epoch in seconds from 2000-01-01T12:00:00 TAI.
     * <p>Beware, it is not {@link #J2000_EPOCH} since it is in TAI and not in TT.</p> */
    double epoch = 0.0;

    /** Offset from the reference epoch in seconds. */
    double offset = 0.0;

public:

    // static const AbsoluteDate J2000_EPOCH = AbsoluteDate();

    /** Create an instance with a default value ({@link #J2000_EPOCH}).
     */
    AbsoluteDate() {
        /* FIME
        epoch  = J2000_EPOCH.epoch;
        offset = J2000_EPOCH.offset;
        */
    }

    AbsoluteDate( const double epoch ) : epoch(epoch) { }

    /** Compute the physically elapsed duration between two instants.
     * <p>The returned duration is the number of seconds physically
     * elapsed between the two instants, measured in a regular time
     * scale with respect to surface of the Earth (i.e either the {@link
     * TAIScale TAI scale}, the {@link TTScale TT scale} or the {@link
     * GPSScale GPS scale}). It is the only method that gives a
     * duration with a physical meaning.</p>
     * <p>This method gives the same result (with less computation)
     * as calling {@link #offsetFrom(AbsoluteDate, TimeScale)}
     * with a second argument set to one of the regular scales cited
     * above.</p>
     * <p>This method is the reverse of the {@link #AbsoluteDate(AbsoluteDate,
     * double)} constructor.</p>
     * @param instant instant to subtract from the instance
     * @return offset in seconds between the two instants (positive
     * if the instance is posterior to the argument)
     * @see #offsetFrom(AbsoluteDate, TimeScale)
     * @see #AbsoluteDate(AbsoluteDate, double)
     */
    double durationFrom(const AbsoluteDate instant) const {
        return (epoch - instant.epoch) + (offset - instant.offset);
    }

    /** Compare the instance with another date.
     * @param date other date to compare the instance to
     * @return a negative integer, zero, or a positive integer as this date
     * is before, simultaneous, or after the specified date.
     */
    int compareTo(const AbsoluteDate date) const {
        const double t = durationFrom(date);
        return (t < 0) ? -1 : (t > 0);
    }

    /*
    AbsoluteDate getDate() {
        return *this;
    }
    */

    bool operator == ( const AbsoluteDate &other ) const {
        return compareTo( other ) == 0;
    }

};


class TimeStamped : public Equatable {
public:
    /** Get the date.
     * @return date attached to the object
     */
    virtual AbsoluteDate getDate() const;


    bool operator == (const TimeStamped& rhs) const {
        return getDate().compareTo(rhs.getDate()) == 0;
    }

    bool operator < (const TimeStamped& rhs) const {
        return getDate().compareTo(rhs.getDate()) < 0;
    }

    bool operator > (const TimeStamped& rhs) const {
        return getDate().compareTo(rhs.getDate()) > 0;
    }

    bool operator <= (const TimeStamped& rhs) const {
        return ( *this < rhs || *this == rhs );
    }

    bool operator >= (const TimeStamped& rhs) const {
        return ( *this > rhs || *this == rhs );
    }


    /*
    bool operator() (const TimeStamped& lhs, const TimeStamped& rhs) const
    {
        return lhs.getDate().compareTo(rhs.getDate()) < 0;
    }
    */
};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_ABSOLUTEDATE_H
