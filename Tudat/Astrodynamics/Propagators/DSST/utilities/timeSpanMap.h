#ifndef TUDAT_DSST_TIMESPANMAP_H
#define TUDAT_DSST_TIMESPANMAP_H

#include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "absoluteDate.h"


namespace tudat
{

namespace propagators
{

namespace dsst
{


/** Local class holding transition times. */
template <class S>
class Transition : public TimeStamped {

private:
    /** Transition date. */
    const AbsoluteDate date;

    /** Entry valid before the transition. */
    const S before;

    /** Entry valid after the transition. */
    const S after;

public:
    Transition() { }

    // Transition( const Transition &t ) : date(t.date), before(t.before), after(t.after) { }

    /** Simple constructor.
     * @param date transition date
     * @param before entry valid before the transition
     * @param after entry valid after the transition
     */
    Transition(const AbsoluteDate date, const S before, const S after) : date(date), before(before), after(after) { }

    /** Get the transition date.
     * @return transition date
     */
    AbsoluteDate getDate() const {
        return date;
    }

    /** Get the entry valid before transition.
     * @return entry valid before transition
     */
    S getBefore() const {
        return before;
    }

    /** Get the entry valid after transition.
     * @return entry valid after transition
     */
    S getAfter() const {
        return after;
    }

    /*
    Transition& operator = ( Transition rhs ) const {
        return *this;
    }
    */

};



struct ChronologicalComparator : public std::binary_function<TimeStamped, TimeStamped, bool> {
    bool operator()(const TimeStamped& lhs, const TimeStamped& rhs) const
    {
        return lhs.getDate().compareTo(rhs.getDate()) < 0;
    }
};


template<class Q, class R>
Q floor( std::set<Q,R> elements, Q target ) {
    for ( auto iter = elements.rbegin(); iter != elements.rend(); ++iter ) {
        auto element = *iter;
        if ( element <= target ) {
            return element;
        }
    }
    return Q();
}

template<class Q, class R>
Q ceil( std::set<Q,R> elements, Q target ) {
    for ( auto iter = elements.begin(); iter != elements.end(); ++iter ) {
        auto element = *iter;
        if ( element >= target ) {
            return element;
        }
    }
    return Q();
}


/** Container for objects that apply to spans of time.

 * @param <T> Type of the data.

 */
template<class T>
class TimeSpanMap {

private:
    /** Container for the data. */
    std::set< Transition<T>, ChronologicalComparator > data;

public:
    TimeSpanMap() { }

    /** Create a map containing a single object, initially valid throughout the timeline.
     * <p>
     * The real validity of this first entry will be truncated as other
     * entries are either {@link #addValidBefore(Object, AbsoluteDate)
     * added before} it or {@link #addValidAfter(Object, AbsoluteDate)
     * added after} it.
     * </p>
     * @param entry entry (initially valid throughout the timeline)
     */
    TimeSpanMap(const T entry)
    {
        data.insert( Transition<T>(AbsoluteDate(), entry, entry) );
    }

    /** Add an entry valid before a limit date.
     * <p>
     * As an entry is valid, it truncates the validity of the neighboring
     * entries already present in the map.
     * </p>
     * <p>
     * The transition dates should be entered only once, either
     * by a call to this method or by a call to {@link #addValidAfter(Object,
     * AbsoluteDate)}. Repeating a transition date will lead to unexpected
     * result and is not supported.
     * </p>
     * @param entry entry to add
     * @param latestValidityDate date before which the entry is valid
     * (sould be different from <em>all</em> dates already used for transitions)
     */
    void addValidBefore(const T entry, const AbsoluteDate latestValidityDate) {

        if (data.size() == 1) {
            const Transition<T> single = *(data.begin());
            if (single.getBefore() == single.getAfter()) {
                // the single entry was a dummy one, without a real transition
                // we replace it entirely
                data.clear();
                data.insert(Transition<T>(latestValidityDate, entry, single.getAfter()));
                return;
            }
        }

        const Transition<T> previous = floor( data, Transition<T>(latestValidityDate, entry, T()) );
        if ( previous == Transition<T>() ) {
            // the new transition will be the first one
            data.insert(Transition<T>(latestValidityDate, entry, (*data.begin()).getBefore()));
        } else {
            // the new transition will be after the previous one
            data.erase(previous);
            data.insert(Transition<T>(previous.getDate(), previous.getBefore(), entry));
            data.insert(Transition<T>(latestValidityDate, entry,                previous.getAfter()));
        }

    }

    /** Add an entry valid after a limit date.
     * <p>
     * As an entry is valid, it truncates the validity of the neighboring
     * entries already present in the map.
     * </p>
     * <p>
     * The transition dates should be entered only once, either
     * by a call to this method or by a call to {@link #addValidBefore(Object,
     * AbsoluteDate)}. Repeating a transition date will lead to unexpected
     * result and is not supported.
     * </p>
     * @param entry entry to add
     * @param earliestValidityDate date after which the entry is valid
     * (sould be different from <em>all</em> dates already used for transitions)
     */
    void addValidAfter(const T entry, const AbsoluteDate earliestValidityDate) {

        if (data.size() == 1) {
            const Transition<T> single = *(data.begin());
            if (single.getBefore() == single.getAfter()) {
                // the single entry was a dummy one, without a real transition
                // we replace it entirely
                data.clear();
                data.insert(Transition<T>(earliestValidityDate, single.getBefore(), entry));
                return;
            }
        }

        const Transition<T> next = ceil( data, Transition<T>(earliestValidityDate, entry, T()) );
        if ( next == Transition<T>() ) {
            // the new transition will be the last one
            data.insert(Transition<T>(earliestValidityDate, (*data.end()).getAfter(), entry));
        } else {
            // the new transition will be before the next one
            data.erase(next);
            data.insert(Transition<T>(earliestValidityDate, next.getBefore(), entry));
            data.insert(Transition<T>(next.getDate(),       entry,            next.getAfter()));
        }

    }

    /** Get the entry valid at a specified date.
     * @param date date at which the entry must be valid
     * @return valid entry at specified date
     */
    T get(const AbsoluteDate date) {
        const Transition<T> previous = floor( data, Transition<T>(date, T(), T()) );
        if ( previous == Transition<T>() ) {
            // there are no transition before the specified date
            // return the first valid entry
            return (*data.begin()).getBefore();
        } else {
            return previous.getAfter();
        }
    }

    /** Get an unmodifiable view of the sorted transitions.
     * @return unmodifiable view of the sorted transitions
     */
    std::set< const Transition<T>, ChronologicalComparator > getTransitions() {
        std::set< const Transition<T>, ChronologicalComparator > transitions;
        for ( auto transition : data ) {
            transitions.insert( transition );
        }
        return transitions;
    }


};



} // namespace dsst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_TIMESPANMAP_H
