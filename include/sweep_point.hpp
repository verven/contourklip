#ifndef CONTOURKLIP_SWEEPPOINT_HPP
#define CONTOURKLIP_SWEEPPOINT_HPP

#include <set>
#include "geometry_base.hpp"
#include "bezier_utils.hpp"

namespace contourklip {
    enum BooleanOpType {
        UNION = 0,
        INTERSECTION = 1,
        DIFFERENCE = 2,
        XOR = 3,
        DIVIDE = 4
    };

    std::ostream &operator<<(std::ostream &o, const BooleanOpType &p) {
        switch (p) {
            case INTERSECTION:
                return o << "intersection";
            case UNION:
                return o << "union";
            case DIFFERENCE:
                return o << "difference";
            case XOR:
                return o << "xor";
            case DIVIDE:
                return o << "divide";
        }
        return o;
    }

    namespace detail {
        enum EdgeType {
            NORMAL, SAME_TRANSITION, DIFFERENT_TRANSITION
        };
        enum PolygonType {
            SUBJECT = 0,
            CLIPPING = 1
        };

        struct SweepPoint {
            SweepPoint() = default;

            bool left = false;
            Point2d point;
            SweepPoint *other_point = nullptr;
            PolygonType ptype{};
            std::size_t contourid = 0;

            //index at which SL iterator to other sp is stored
            std::size_t other_iter_pos = 0;
            bool in_out = false; // if for a ray passing upwards into the edge, it is an inside-outside transition
            bool other_in_out = false; // inout for the closest edge downward in SL that is from the other polygon

            // indicates if the edge associated to this point is in the result contour of the clipping
            // operation given by the index. 0 = default, 1 = INTERSECTION, 2 = DIFFERENCE
            std::array<bool, 3> in_result{};
            // follows the same principle, but instead maps to the previous SweepPoint* (downwards)
            // which is in the result.
            std::array<SweepPoint *, 3> prev_in_result{};

            //the following fields are used when connecting the edges
            std::size_t pos = 0; // position in the result array
            bool result_in_out = false; //if the associated edge is an in out transition into its result contour
            std::size_t result_contour_id = 0;
            EdgeType edgetype = NORMAL;

            //used if segment is curve
            bool curve = false;
            bool islast = false;
            Point2d controlp;
            Point2d initial_controlp;
            SweepPoint *start = nullptr;

            SweepPoint(bool left, const Point2d &point,
                       SweepPoint *otherPoint) :
                    left(left), point(point), other_point(otherPoint) {}

            explicit SweepPoint(const Point2d &point) :
                    point(point) {}

            bool vertical() const { return point.x() == other_point->point.x(); }

            void set_if_left() {
                left = increasing(point, other_point->point);
                this->other_point->left = !left;
            }
        };

        std::ostream &operator<<(std::ostream &o, const SweepPoint &p) {
            if (p.other_point) {
                o << "[" << p.point << "->" << p.other_point->point;
                if (p.left && p.curve) o << "\n--->c" << p.controlp << "->" << p.other_point->controlp;
                o << ", l " << p.left
                  << ", res " << p.in_result[0]
                  << ", ptype " << p.ptype
                  << ", cid " << p.result_contour_id
                  << ", c " << p.curve
                  << "]";
                return o;
            } else {
                return o << "[" << p.point << "->" << " [nullptr] " << p.left << "]";
            }
        }

        bool overlapping(const SweepPoint *e1, const SweepPoint *e2) {
            if (e1->curve != e2->curve) {
                return false;
            }
            bool overlapping_ends = e1->point == e2->point
                                    && e1->other_point->point == e2->other_point->point;
            if (!e1->curve) {
                return overlapping_ends;
            }

            return overlapping_ends
                   && e1->controlp == e2->controlp
                   && e1->other_point->controlp == e2->other_point->controlp;
        }

        //Returns true iff the segment associated with e1 is below a point p.
        auto curve_below_point = [](const SweepPoint *e1, const Point2d &p) -> bool {

            CubicBezier c{e1->point, e1->controlp,
                          e1->other_point->controlp, e1->other_point->point};

            double a = 0., b = 1.;
            Point2d left = c.p0;
            Point2d right = c.p3;
            Point2d sample;
            while ((b - a) > 1e-10) {
                if (left.y() < p.y() && right.y() < p.y()) {
                    return true;
                }
                if (left.y() >= p.y() && right.y() >= p.y()) {
                    return false;
                }
                double mid = 0.5 * (a + b);
                sample = beziermap(c, mid);
                if (sample.x() < p.x()) {
                    a = mid;
                    left = sample;
                } else {
                    b = mid;
                    right = sample;
                }
            }
            return sample.y() < p.y();
        };


        // check if the associated bezier curve of the first is below the associated bezier curve of the second.
        // preconditions: both are curves, share the endpoint (may be start or end), and otherwise do not intersect.
        template<typename Orient2DF = LeftOfLine, typename CollinearF = IsCollinear>
        bool curve_below(const SweepPoint *e1, const SweepPoint *e2, CollinearF f = {}) {
            if (!f(e1->point, e1->controlp, e2->controlp)
                && e1->point != e1->controlp) { // special case with vanishing derivative
                if (vertical(e1->point, e1->controlp)) {
                    return e1->controlp.y() < e1->point.y();
                }
                return above_line<Orient2DF>(e1->point, e1->controlp, e2->controlp);
            }
            const SweepPoint *sp = e1;
            const SweepPoint *other = e2;
            bool reversed;
            // we want to sample the point on the curve which is shortest in the x direction
            if ((reversed = e1->left == (e1->other_point->point.x()
                                         > e2->other_point->point.x()))) {
                std::swap(sp, other);
            }
            Point2d a = beziermap(sp->point, sp->controlp,
                                  sp->other_point->controlp, sp->other_point->point, 0.5);

            CubicBezier tmp{other->point, other->controlp,
                            other->other_point->controlp, other->other_point->point};
            double t = t_from_x(tmp, a.x());

            double y_other = beziermap(tmp, t).y();
            return reversed ? a.y() > y_other : a.y() < y_other;
        }

        // This is the comparator used for the sweeppoint queue. returns true iff e1 < e2.
        template<typename Orient2DF = LeftOfLine, typename CollinearF = IsCollinear>
        bool queue_comp(const SweepPoint *e1, const SweepPoint *e2) {
            CollinearF collinear{};
            if (e1->point.x() != e2->point.x())
                return e1->point.x() < e2->point.x();
            if (e1->point.y() !=
                e2->point.y())
                return e1->point.y() < e2->point.y();

            if (e1->left != e2->left) {
                //right endpoint is processed first.
                return e2->left;
            }

            if (overlapping(e1, e2)) {
                return e1->ptype < e2->ptype;
            }
            // Both events represent lines
            if (!e1->curve && !e2->curve) {
                if (!collinear(e1->point, e1->other_point->point, e2->other_point->point)) {
                    // the event associate to the bottom segment is processed first
                    return above_line<Orient2DF>(e1->point, e1->other_point->point, e2->other_point->point);
                }
                return e1->ptype < e2->ptype;
            }

            // very special case where curve_below fails to differentiate due to round-off.
            // In that case we use some consistent criterion.
            // At this point the segments do not exactly overlap.
            bool a = curve_below<Orient2DF, CollinearF>(e1, e2);
            bool b = curve_below<Orient2DF, CollinearF>(e2, e1);
            if (a == b) {
                if (e1->ptype != e2->ptype) return e1->ptype < e2->ptype;
                if (e1->other_point->point != e2->other_point->point)
                    return increasing(e1->other_point->point, e2->other_point->point);
                if (e1->controlp != e2->controlp) return increasing(e1->controlp, e2->controlp);
                return increasing(e1->other_point->controlp, e2->other_point->controlp);
            }
            //at least one point is from a curve
            return a;
        }

        // Comparator used for the sweep line. returns true iff le1 < le2.
        // Note that only left events can be in the sweep line.
        template<typename Orient2DF, typename CollinearF>
        struct SComp {
            bool operator()(const SweepPoint *le1, const SweepPoint *le2) const {
                if (le1 == le2)
                    return false;

                if (overlapping(le1, le2)) {
                    return le1->ptype < le2->ptype;
                }

                if (le1->point == le2->point) {
                    return queue_comp<Orient2DF, CollinearF>(le1, le2);
                }

                if (le1->point.x() == le2->point.x()) {
                    return le1->point.y() < le2->point.y();
                }

                if (!le1->curve && !le2->curve) {
                    if (queue_comp<Orient2DF, CollinearF>(le1, le2)) {
                        return above_line<Orient2DF>(le1->point,
                                                     le1->other_point->point, le2->point);
                    }
                    return above_line<Orient2DF>(le2->point,
                                                 le1->point, le2->other_point->point);
                }

                //one of the segments is a curve.
                if (queue_comp<Orient2DF, CollinearF>(le1, le2)) {
                    // le1 has been inserted first.
                    return curve_below_point(le1, le2->point);
                }
                return !curve_below_point(le2, le1->point);
            }
        };
    }
}
#endif //CONTOURKLIP_SWEEPPOINT_HPP