#ifndef CONTOURKLIP_BEZIER_UTILS_HPP
#define CONTOURKLIP_BEZIER_UTILS_HPP

#ifdef DEBUG
#include <iostream>
#endif
#include <array>
#include "geometry_base.hpp"
#include "polynomial_solver.hpp"
#include "direct_solvers.hpp"

namespace contourklip::detail {

        template<typename T>
        T beziermap(const T &p0, const T &p1, const T &p2, const T &p3, const T &t) {
            return (1 - t) * (1 - t) * (1 - t) * p0 + 3 * (1 - t) * (1 - t) * t * p1 + 3 * (1 - t) * t * t * p2 +
                   t * t * t * p3;
        }

        Point2d beziermap(const CubicBezier &c, const double t) {
            return {
                    beziermap(c.p0.x(), c.p1.x(), c.p2.x(), c.p3.x(), t),
                    beziermap(c.p0.y(), c.p1.y(), c.p2.y(), c.p3.y(), t)
            };
        }

        Point2d beziermap(const Point2d &p0, const Point2d &p1, const Point2d &p2, const Point2d &p3, const double t) {
            return {
                    beziermap(p0.x(), p1.x(), p2.x(), p3.x(), t),
                    beziermap(p0.y(), p1.y(), p2.y(), p3.y(), t)
            };
        }

        std::pair<CubicBezier, CubicBezier> casteljau_split(const CubicBezier &c, double t) {

            double x0 = c.p0.x(), y0 = c.p0.y();
            double x1 = c.p1.x(), y1 = c.p1.y();
            double x2 = c.p2.x(), y2 = c.p2.y();
            double x3 = c.p3.x(), y3 = c.p3.y();

            double x12 = (x1 - x0) * t + x0;
            double y12 = (y1 - y0) * t + y0;

            double x23 = (x2 - x1) * t + x1;
            double y23 = (y2 - y1) * t + y1;

            double x34 = (x3 - x2) * t + x2;
            double y34 = (y3 - y2) * t + y2;

            double x123 = (x23 - x12) * t + x12;
            double y123 = (y23 - y12) * t + y12;

            double x234 = (x34 - x23) * t + x23;
            double y234 = (y34 - y23) * t + y23;

            //instead of computing it from the above values, this way ensures consistency wrt. roundoff.
            Point2d split = beziermap(c, t);

            return {
                    CubicBezier{{x0, y0}, {x12, y12}, {x123, y123}, split},
                    CubicBezier{split, {x234, y234}, {x34, y34}, {x3, y3}}
            };
        }

        Point2d scaleto(const Point2d &a, const Point2d &b, double scale) {

            double d = dist(a, b);
            double vx = b.x() - a.x(), vy = b.y() - a.y();
            return {scale * 1. / d * vx + a.x(), scale * 1. / d * vy + a.y()};
        }

        CubicBezier bezier_merge(const CubicBezier &a, const CubicBezier &b) {

            double l_b = dist(b.p0, b.p1);
            double r_a = dist(a.p2, a.p3);
            double ratio = l_b / r_a;

            double q1_x = (ratio + 1) * a.p1.x() - ratio * a.p0.x();
            double q1_y = (ratio + 1) * a.p1.y() - ratio * a.p0.y();

            double q2_x = (1. / ratio) * ((ratio + 1) * b.p2.x() - b.p3.x());
            double q2_y = (1. / ratio) * ((ratio + 1) * b.p2.y() - b.p3.y());

            return {a.p0, {q1_x, q1_y}, {q2_x, q2_y}, b.p3};
        }

        std::optional<CubicBezier> sub_bezier(const CubicBezier &c, double t, double u) {
#ifdef DEBUG
            assert( ! (c.p0 == c.p1 && c.p3 == c.p2) );
            assert(t != u);
#endif

            auto[l1, r1] = casteljau_split(c, t);
            if (u == 1) {
                return r1;
            }
            auto[l2, r2] = casteljau_split(c, u);
            if (t == 0) {
                return l2;
            }

            // we have a quadratic bezier
            if (c.p3 == c.p2 || c.p0 == c.p1) {
                Point2d inter = basic_intersection(r1.p0, r1.p1, l2.p2, l2.p3);
                Point2d p1 = linear_map(r1.p0, inter, 2. / 3.);
                Point2d p2 = linear_map(inter, l2.p3, 1. / 3.);
                return {{r1.p0, p1, p2, l2.p3}};
            }

            double d = dist(r2.p0, r2.p1);
            double rx = r1.p2.x() - r2.p2.x();
            double ry = r1.p2.y() - r2.p2.y();
            double ratio;
            if (std::abs(rx) > std::abs(ry)) {
                ratio = (r2.p2.x() - r2.p3.x()) / (rx);
            } else {
                ratio = (r2.p2.y() - r2.p3.y()) / (ry);
            }

            if (std::isnan(ratio) || ratio == 0) {
#ifdef DEBUG
                std::cout << "ratio is nan, input is: " << c << ", [" << t << ", " << u << "] " << '\n';
                std::cout << "d, ratio:,  " << d << " " << ratio << '\n';
                std::cout << "r1, r2 " << r1 <<",  " << r2 << '\n';
                std::abort();
#endif
                return {};
            }

            double right_dist = d / ratio;
            double p1_x = (r1.p1.x() + ratio * r1.p0.x()) / (ratio + 1.);
            double p1_y = (r1.p1.y() + ratio * r1.p0.y()) / (ratio + 1.);
            Point2d p1{p1_x, p1_y};
            Point2d p2 = scaleto(l2.p3, l2.p2, right_dist);

            return {CubicBezier{l1.p3, p1, p2, r2.p0}};
        }

        inline double control_poly_area(const CubicBezier &c) {
            return quadri_area(c.p0, c.p1, c.p2, c.p3);
        }

        void reverse(CubicBezier &c) {
            Point2d tmp = c.p0;
            c.p0 = c.p3;
            c.p3 = tmp;
            tmp = c.p1;
            c.p1 = c.p2;
            c.p2 = tmp;
        }

//        bool degenerate(const CubicBezier &c) {
//            return is_collinear(c.p0, c.p1, c.p2) && is_collinear(c.p1, c.p2, c.p3);
//        }

        std::pair<std::array<double, 4>, std::array<double, 4>> bezier_transpose(
                const CubicBezier &c) {
            return {{c.p0.x(), c.p1.x(), c.p2.x(), c.p3.x()},
                    {c.p0.y(), c.p1.y(), c.p2.y(), c.p3.y()}};
        }

        template<typename PointMap>
        CubicBezier transform(const CubicBezier &c, PointMap &f) {
            CubicBezier out{};
            out.p0 = f(c.p0);
            out.p1 = f(c.p1);
            out.p2 = f(c.p2);
            out.p3 = f(c.p3);
            return out;
        }

        inline bool unit_inter(double t) {
            return 0 <= t && t <= 1;
        }

        inline Point2d rotate_pt(const Point2d &in, const Point2d &rpoint, double cos_val, double sin_val) {
            double offx = in.x() - rpoint.x(), offy = in.y() - rpoint.y();
            return {offx * cos_val - offy * sin_val + rpoint.x(),
                    offx * sin_val + offy * cos_val + rpoint.y()};
        }

        inline Point2d
        rotate_translate_pt(const Point2d &in, const Point2d &rpoint, const Point2d &subtract, double cos_val,
                            double sin_val) {
            auto t = rotate_pt(in, rpoint, cos_val, sin_val);
            return {t.x() - subtract.x(), t.y() - subtract.y()};
        }

        std::pair<double, double>
        axis_align_rotate_vals(const Point2d p_start, const Point2d p_end, double abstol = 1e-8) {
            if (std::abs(p_start.y() - p_end.y()) < abstol) {
                return {1., 0.};
            } else {
                double vx = (p_end.x() - p_start.x()), vy = (p_end.y() - p_start.y());
                double cos_val, sin_val;
                if (std::abs(vx) < abstol) {
                    cos_val = 0;
                    sin_val = vy > 0 ? -1 : 1;
                } else {

                    double d = -vy / vx;
                    cos_val = 1. / std::sqrt(1 + d * d);
                    sin_val = d * cos_val;
                }
                return {cos_val, sin_val};
            }
        }

        template<typename T>
        std::tuple<T, T, T, T> cubic_line_coeffs(const T &a, const T &b, const T &c, const T &d) {
            T a_0 = a;
            T a_1 = -3 * a + 3 * b;
            T a_2 = 3 * a - 6 * b + 3 * c;
            T a_3 = -1 * a + 3 * b - 3 * c + d;
            return {a_0, a_1, a_2, a_3};
        }

        template<int D>
        std::tuple<double, double, double, double> cubic_line_coeffs(const CubicBezier &c) {
            static_assert(D == 0 || D == 1);
            double a_0 = get<D>(c.p0);
            double a_1 = -3 * get<D>(c.p0) + 3 * get<D>(c.p1);
            double a_2 = 3 * get<D>(c.p0) - 6 * get<D>(c.p1) + 3 * get<D>(c.p2);
            double a_3 = -1 * get<D>(c.p0) + 3 * get<D>(c.p1) - 3 * get<D>(c.p2) + get<D>(c.p3);
            return {a_0, a_1, a_2, a_3};
        }

        // returns the smallest parametric t value in the bezier interval such that B_x(t) = x.
        double t_from_x(const CubicBezier &c, double x, double cubic_tol = 1e-9) {
            auto[a_0, a_1, a_2, a_3] = cubic_line_coeffs<0>(c);
            a_0 -= x;
            double val = 1;
            auto rootprocess = [&val](double t) {
                if (in_range(t, 0, 1)) {
                    val = std::min(val, t);
                }
            };
            directsolvers::solve_cubic_real(a_0, a_1, a_2, a_3, rootprocess, cubic_tol);
            return val;
        }

        template<typename ConsumerF>
        void t_vals_from_x(const CubicBezier &c, double x, ConsumerF &f, double cubic_tol = 1e-9) {
            auto[a_0, a_1, a_2, a_3] = cubic_line_coeffs<0>(c);
            a_0 -= x;
            directsolvers::solve_cubic_real(a_0, a_1, a_2, a_3, f, cubic_tol);
        }

        bool curve_below(const CubicBezier &c, const Point2d &p) {
            double t = t_from_x(c, p.x());
            return beziermap(c, t).y() < p.y();
        }

        bool bezier_direction(const CubicBezier &c) {
            double area = quadri_area(c.p0, c.p1, c.p2, c.p3);
            if (area == 0) {
                return increasing(c.p0, c.p3);
            }
            return area > 0;
        }

        enum Extremity_Direction {
            X_Extremity = 0,
            Y_Extremity = 1
        };

        template<Extremity_Direction D, typename OutF>
        void bezier_extremities(const CubicBezier &bt, OutF f) {
            static_assert(D == 0 || D == 1);

            //quadratic coefficients
            double a = (-get<D>(bt.p0) + 3 * get<D>(bt.p1) + -3 * get<D>(bt.p2) + get<D>(bt.p3));
            double b = 2 * (get<D>(bt.p0) - 2 * get<D>(bt.p1) + get<D>(bt.p2));
            double c = get<D>(bt.p1) - get<D>(bt.p0);

            std::pair<double, double> res{};
            switch (directsolvers::solve_quadratic(c, b, a, res)) {
                case 1:
                    if (unit_inter(res.first)) {
                        f(res.first);
                    }
                    break;
                case 2:
                    if (unit_inter(res.first)) {
                        f(res.first);
                    }
                    if (unit_inter(res.second)) {
                        f(res.second);
                    }
                    break;
                default:
                    return;
            }
        }

        //returns sorted t values ti of a bezier such that each curve subinterval
        // [t_i, t_i+1] is weakly monotonic (increasing/decreasing) in the x and in the y
        // direction. Only inner values (different from 0 and 1) are returned.
        // Note that there can be max. 5 subcurves, which follows from the fact
        // that there are max 4 roots.
        template<typename OutF>
        void bezier_monotonic_split(const CubicBezier &bt, OutF f, double abstol = 1e-10) {

            if ((in_box(bt.p0, bt.p3, bt.p1) || bt.p0 == bt.p1)
                && (in_box(bt.p0, bt.p3, bt.p2) || bt.p3 == bt.p2)) {
                return;
            }

            struct ExtremityValue {
                double t = 1.5;
                Extremity_Direction direction;
            };

            std::array<ExtremityValue, 4> t_vals{};
            std::size_t num = 0;
            Extremity_Direction d = X_Extremity;

            auto t_report = [&t_vals, &num, &d](double t) {
                t_vals[num++] = {t, d};
            };

            bezier_extremities<X_Extremity>(bt, t_report);
            d = Y_Extremity;
            bezier_extremities<Y_Extremity>(bt, t_report);

            std::sort(t_vals.begin(), t_vals.end(),
                      [](const ExtremityValue &a, const ExtremityValue &b) { return a.t < b.t; });

            for (std::size_t i = 0; i < num; ++i) {
                // we want values strictly inside the curve interval. What is more, in the case where
                // a control point_ overlaps with an endpoint, the derivative is also 0, even though it is not
                // an extrema.
                if ((t_vals[i].t < abstol) || (t_vals[i].t + abstol > 1)) {
                    continue;
                }
                f(t_vals[i].t, t_vals[i].direction);
            }
        }

        bool in_box_special(const Segment &seg, const Point2d &p) {
            bool x_ok = in_range_closed(std::min(seg.first.x(), seg.second.x()),
                                        std::max(seg.first.x(), seg.second.x()), p.x());
            bool y_ok = in_range_closed(std::min(seg.first.y(), seg.second.y()),
                                        std::max(seg.first.y(), seg.second.y()), p.y());
            return ((x_ok && y_ok)
                    || (seg.first.x() == seg.second.x() && y_ok)
                    || (seg.first.y() == seg.second.y() && x_ok)
            );
        }

        template<typename IntersectionFunctor>
        void line_bezier_inter_impl(const Segment &seg, const CubicBezier &c,
                                    IntersectionFunctor &inter, double abstol = 1e-10, double cubic_tol = 1e-9) {

            //robustness: we first check if we have some overlapping points/points that lie on the segment.
            //Then we can use exact t=0, 1 and ignore these cases later.
            bool overlapping_start = c.p0 == seg.first || c.p0 == seg.second,
                    overlapping_end = c.p3 == seg.first || c.p3 == seg.second;

            std::size_t start_offset = overlapping_start, end_offset = overlapping_end;

            auto rotate = axis_align_rotate_vals(seg.first, seg.second);
            auto normalize_f = [rotate, seg](const Point2d &p) -> Point2d {
                return rotate_translate_pt(p, seg.first, seg.first, rotate.first, rotate.second);
            };
            CubicBezier normalized = transform(c, normalize_f);

            auto[a_0, a_1, a_2, a_3] = cubic_line_coeffs<1>(normalized);

            std::array<double, 3> roots{1.5, 1.5, 1.5};
            std::size_t num = 0;
            auto addroot = [&roots, &num, &abstol](double r) {
                if (std::abs(r) < abstol) {
                    if (num == 0 || (roots[num - 1] != 0)) {
                        roots[num++] = 0;
                    }
                    return;
                }
                if ((r > 1 && r < 1 + abstol) || (r < 1 && r + abstol > 1)) {
                    if (num == 0 || (roots[num - 1] != 1)) {
                        roots[num++] = 1;
                    }
                    return;
                }
                if (r >= 0 && r <= 1) {
                    if (num == 0 || (roots[num - 1] != r)) {
                        roots[num++] = r;
                    }
                    return;
                }
            };
            //note that roots are sorted in increasing order
            directsolvers::solve_cubic_real(a_0, a_1, a_2, a_3, addroot, cubic_tol);
            // we take into account overlapping start/end points and points on the segment,
            // to produce a consistent output.
            std::size_t bound = num >= end_offset ? num - end_offset : 0;
            for (std::size_t i = start_offset; i < bound; ++i) {
                double r = roots[i];
                if (in_range_closed(0, 1, r)) {
                    //checking if intersection in segment interval.
                    Point2d mapped = beziermap(c, r);
                    if (in_box_special(seg, mapped)) {
//                        if(make_bbox(seg).strict_contains(mapped)){
                        double segment_t = segment_tval(seg, mapped);
#ifdef DEBUG
                        //                    std::cout << "intersection at " << segment_t << ", "
                        //                    << r << ", " << mapped << '\n';
                        //                    std::cout << "with data " << seg << " " << c << '\n';
#endif
                        inter(segment_t, r, mapped);
                    }
                }
            }
        }

        template<typename IntersectionFunctor>
        void line_bezier_inter(const Segment &seg, const CubicBezier &c,
                               IntersectionFunctor &inter_f, double abstol = 1e-10, double cubic_tol = 1e-9) {

            BBox box{};
            make_hull_bbox(c, box);

            bool nointersect = !make_bbox(seg).weak_overlap(box);
            if (nointersect) {
                return;
            }

            Segment seg2 = seg;
            CubicBezier c2 = c;
            bool seg_reversed, bezier_reversed;
            if ((seg_reversed = !increasing(seg.first, seg.second))) {
                seg2 = {seg.second, seg.first};
            }
            if ((bezier_reversed = !bezier_direction(c))) {
                reverse(c2);
            }
            auto inter_f2 = [&](double t, double u, const Point2d &p) {
                t = seg_reversed ? 1 - t : t;
                u = bezier_reversed ? 1 - u : u;
                inter_f(t, u, p);
            };
            line_bezier_inter_impl(seg2, c2, inter_f2, abstol, cubic_tol);
        }

/*
 *
 *  the remainder of this file is concerned with (cubic) bezier-bezier intersections.
 *
*/

        //defines an ordering between 2 different bezier curves.
        bool bezier_ordering(const CubicBezier &b1, const CubicBezier &b2) {
            double area1 = control_poly_area(b1);
            double area2 = control_poly_area(b2);
            if (area1 != area2) {
                return area1 < area2;
            }
            if (b1.p0 != b2.p0) {
                return increasing(b1.p0, b2.p0);
            }
            if (b1.p1 != b2.p1) {
                return increasing(b1.p1, b2.p1);
            }
            if (b1.p2 != b2.p2) {
                return increasing(b1.p2, b2.p2);
            }
            return increasing(b1.p3, b2.p3);
        }

        double determinant3(double a, double b, double c,
                            double d, double e, double f,
                            double g, double h, double i) {
            using namespace directsolvers;
            return a * diff_of_products(e, i, f, h)
                   + b * diff_of_products(f, g, d, i)
                   + c * diff_of_products(d, h, e, g);
        }

        std::array<double, 3> determinant_coeffs(double a, double b, double c, double d) {
            return std::array<double, 3>{(b - d), (-a + c), directsolvers::diff_of_products(a, d, b, c)};
        }

        // a simple functor which maps a point to the parametric interval
        // using a rational function
        template<typename T>
        class maptoT {
            T inverta = 0;
            T invertb = 1;
            std::array<T, 3> a;
            std::array<T, 3> b;
        public:
            maptoT(std::array<T, 3> a, std::array<T, 3> b, bool isinverted) :
                    a(a), b(b) {
                if (isinverted) {
                    invertb = -1;
                    inverta = 1;
                }
            }

            T operator()(T x, T y) {
                return inverta + invertb * (a[0] * x + a[1] * y + a[2]) / (b[0] * x + b[1] * y + b[2]);
            }
        };

        template<typename collinearF = IsCollinear>
        maptoT<double> curveinverter(CubicBezier c, collinearF f = {}) {
            bool directionInverted = false;
            if ((directionInverted = f(c.p1, c.p2, c.p3))) {
                reverse(c);
            }

            auto mult_c = [](double c, const std::array<double, 3> &a) -> std::array<double, 3> {
                return {a[0] * c, a[1] * c, a[2] * c};
            };

            auto l31 = mult_c(3.0, determinant_coeffs(c.p3.x(), c.p3.y(), c.p1.x(), c.p1.y()));
            auto l30 = determinant_coeffs(c.p3.x(), c.p3.y(), c.p0.x(), c.p0.y());
            auto l21 = mult_c(9.0, determinant_coeffs(c.p2.x(), c.p2.y(), c.p1.x(), c.p1.y()));
            auto l20 = mult_c(3.0, determinant_coeffs(c.p2.x(), c.p2.y(), c.p0.x(), c.p0.y()));
            auto l10 = mult_c(3.0, determinant_coeffs(c.p1.x(), c.p1.y(), c.p0.x(), c.p0.y()));

            double d = 3.0 * determinant3(c.p1.x(), c.p1.y(), 1, c.p2.x(),
                                          c.p2.y(), 1, c.p3.x(), c.p3.y(), 1);

            double c1 = determinant3(c.p0.x(), c.p0.y(), 1, c.p1.x(), c.p1.y(), 1,
                                     c.p3.x(), c.p3.y(), 1) / d;
            double c2 = -1.0 * determinant3(c.p0.x(), c.p0.y(), 1, c.p2.x(),
                                            c.p2.y(), 1, c.p3.x(), c.p3.y(), 1) / d;

            std::array<double, 3> num_coeffs{}, den_coeffs{};
            auto compute_num = [&](std::size_t i) -> double {
                return (c1 * l30[i]) + (c2 * l20[i]) + l10[i];
            };
            auto compute_den = [&](std::size_t i) -> double {
                return num_coeffs[i] - (c2 * (l30[i] + l21[i]) + l20[i] + (c1 * l31[i]));
            };
            num_coeffs[0] = compute_num(0);
            num_coeffs[1] = compute_num(1);
            num_coeffs[2] = compute_num(2);
            den_coeffs[0] = compute_den(0);
            den_coeffs[1] = compute_den(1);
            den_coeffs[2] = compute_den(2);

            maptoT out{num_coeffs, den_coeffs, directionInverted};
            return out;
        }

        template<typename T>
        std::array<T, 10>
        generate_coefficients( const T &b0x, const T &b0y, const T &b1x, const T &b1y,
                               const T &b2x, const T &b2y, const T &b3x, const T &b3y,
                               const T &px1, const T &py1, const T &px2, const T &py2, const T &px3, const T &py3) noexcept {
            auto[ax, bx, cx, dx]
            = cubic_line_coeffs(b0x, b1x, b2x, b3x);
            auto[ay, by, cy, dy]
            = cubic_line_coeffs(b0y, b1y, b2y, b3y);

            // note: the following code has been automatically generated.

            auto Power = [](const T &t, const T &v) -> T {
                if (v == 2) return t * t;
                return t * t * t;
            };

            T tmp1 = 3*px1-3*px2+px3;
            T tmp2 = (3*py1-3*py2+py3);

            // clang-format off

            T c0=Power(ay,3)*Power(tmp1,3)-3*Power(ay,2)*(18*Power(px1,3)*py3+ax*Power(tmp1,2)*tmp2-3*px2*(Power(px3,2)*py1-3*px2*px3*py2+3*Power(px2,2)*py3)+3*px1*(2*Power(px3,2)*(3*py1-py2)+9*Power(px2,2)*(py1+py3)+3*px2*px3*(-3*py1-3*py2+py3))-9*Power(px1,2)*(2*px3*(py1-3*py2+py3)+3*px2*(py2+py3)))+3*ay*(Power(ax,2)*(tmp1)*Power(3*py1-3*py2+py3,2)+9*px1*(Power(px3,2)*Power(py1,2)+py3*(3*Power(px2,2)*py1-3*px1*px2*py2+Power(px1,2)*py3)+px3*(-3*px2*py1*py2+3*px1*Power(py2,2)-2*px1*py1*py3))+3*ax*(3*Power(px3,2)*py1*(2*py1-py2)+Power(px2,2)*(9*Power(py1,2)+9*py1*py3-6*py2*py3)+Power(px1,2)*(-9*Power(py2,2)+9*py2*py3+6*(2*py1-py3)*py3)-px1*px3*(12*Power(py1,2)-27*py1*py2+py2*(9*py2+py3))+px2*(3*px1*py3*(-9*py1+3*py2+py3)+px3*(-9*Power(py1,2)-9*py1*py2+6*Power(py2,2)+py1*py3))))-ax*(Power(ax,2)*Power(3*py1-3*py2+py3,3)+27*py1*(Power(px3,2)*Power(py1,2)+py3*(3*Power(px2,2)*py1-3*px1*px2*py2+Power(px1,2)*py3)+px3*(-3*px2*py1*py2+3*px1*Power(py2,2)-2*px1*py1*py3))-9*ax*(3*px3*(2*py1-py2)*(Power(py1,2)+Power(py2,2)-py1*(py2+py3))+px2*(-9*Power(py1,2)*(py2-2*py3)+3*Power(py2,2)*py3-py1*py3*(9*py2+2*py3))-px1*(6*Power(py1,2)*py3+py2*Power(py3,2)+py1*(-9*Power(py2,2)+9*py2*py3-6*Power(py3,2)))));T c1=-3*(Power(ay,2)*Power(tmp1,2)*(-(by*(tmp1))+bx*tmp2)+Power(ax,2)*Power(3*py1-3*py2+py3,2)*(-(by*(tmp1))+bx*tmp2)-9*(by*px1-bx*py1)*(Power(px3,2)*Power(py1,2)+py3*(3*Power(px2,2)*py1-3*px1*px2*py2+Power(px1,2)*py3)+px3*(-3*px2*py1*py2+3*px1*Power(py2,2)-2*px1*py1*py3))-3*ax*by*(3*Power(px3,2)*py1*(2*py1-py2)+Power(px2,2)*(9*Power(py1,2)+9*py1*py3-6*py2*py3)+Power(px1,2)*(-9*Power(py2,2)+9*py2*py3+6*(2*py1-py3)*py3)-px1*px3*(12*Power(py1,2)-27*py1*py2+py2*(9*py2+py3))+px2*(3*px1*py3*(-9*py1+3*py2+py3)+px3*(-9*Power(py1,2)-9*py1*py2+6*Power(py2,2)+py1*py3)))+6*ax*bx*(-3*px3*(2*py1-py2)*(Power(py1,2)+Power(py2,2)-py1*(py2+py3))+px2*(9*Power(py1,2)*(py2-2*py3)-3*Power(py2,2)*py3+py1*py3*(9*py2+2*py3))+px1*(6*Power(py1,2)*py3+py2*Power(py3,2)+py1*(-9*Power(py2,2)+9*py2*py3-6*Power(py3,2))))+ay*(2*ax*(tmp1)*tmp2*(by*(tmp1)-bx*tmp2)+6*by*(6*Power(px1,3)*py3-px2*(Power(px3,2)*py1-3*px2*px3*py2+3*Power(px2,2)*py3)+px1*(2*Power(px3,2)*(3*py1-py2)+9*Power(px2,2)*(py1+py3)+3*px2*px3*(-3*py1-3*py2+py3))-3*Power(px1,2)*(2*px3*(py1-3*py2+py3)+3*px2*(py2+py3)))-3*bx*(3*Power(px3,2)*py1*(2*py1-py2)+Power(px2,2)*(9*Power(py1,2)+9*py1*py3-6*py2*py3)+Power(px1,2)*(-9*Power(py2,2)+9*py2*py3+6*(2*py1-py3)*py3)-px1*px3*(12*Power(py1,2)-27*py1*py2+py2*(9*py2+py3))+px2*(3*px1*py3*(-9*py1+3*py2+py3)+px3*(-9*Power(py1,2)-9*py1*py2+6*Power(py2,2)+py1*py3)))));T c2=-3*(Power(ay,2)*Power(tmp1,2)*(-(cy*(tmp1))+cx*tmp2)+Power(ax,2)*Power(3*py1-3*py2+py3,2)*(-(cy*(tmp1))+cx*tmp2)-ay*(Power(by,2)*Power(tmp1,3)-54*cy*px1*Power(px2,2)*py1+36*cy*Power(px1,2)*px3*py1+54*cy*px1*px2*px3*py1-36*cy*px1*Power(px3,2)*py1+6*cy*px2*Power(px3,2)*py1+27*Power(bx,2)*px1*Power(py1,2)-27*Power(bx,2)*px2*Power(py1,2)+27*cx*Power(px2,2)*Power(py1,2)+9*Power(bx,2)*px3*Power(py1,2)-36*cx*px1*px3*Power(py1,2)-27*cx*px2*px3*Power(py1,2)+18*cx*Power(px3,2)*Power(py1,2)+54*cy*Power(px1,2)*px2*py2-108*cy*Power(px1,2)*px3*py2+54*cy*px1*px2*px3*py2-18*cy*Power(px2,2)*px3*py2+12*cy*px1*Power(px3,2)*py2-54*Power(bx,2)*px1*py1*py2+54*Power(bx,2)*px2*py1*py2-18*Power(bx,2)*px3*py1*py2+81*cx*px1*px3*py1*py2-27*cx*px2*px3*py1*py2-9*cx*Power(px3,2)*py1*py2+27*Power(bx,2)*px1*Power(py2,2)-27*cx*Power(px1,2)*Power(py2,2)-27*Power(bx,2)*px2*Power(py2,2)+9*Power(bx,2)*px3*Power(py2,2)-27*cx*px1*px3*Power(py2,2)+18*cx*px2*px3*Power(py2,2)-36*cy*Power(px1,3)*py3+54*cy*Power(px1,2)*px2*py3-54*cy*px1*Power(px2,2)*py3+18*cy*Power(px2,3)*py3+36*cy*Power(px1,2)*px3*py3-18*cy*px1*px2*px3*py3+18*Power(bx,2)*px1*py1*py3+36*cx*Power(px1,2)*py1*py3-18*Power(bx,2)*px2*py1*py3-81*cx*px1*px2*py1*py3+27*cx*Power(px2,2)*py1*py3+6*Power(bx,2)*px3*py1*py3+3*cx*px2*px3*py1*py3-18*Power(bx,2)*px1*py2*py3+27*cx*Power(px1,2)*py2*py3+18*Power(bx,2)*px2*py2*py3+27*cx*px1*px2*py2*py3-18*cx*Power(px2,2)*py2*py3-6*Power(bx,2)*px3*py2*py3-3*cx*px1*px3*py2*py3+3*Power(bx,2)*px1*Power(py3,2)-18*cx*Power(px1,2)*Power(py3,2)-3*Power(bx,2)*px2*Power(py3,2)+9*cx*px1*px2*Power(py3,2)+Power(bx,2)*px3*Power(py3,2)-2*bx*by*Power(tmp1,2)*tmp2-2*ax*(tmp1)*tmp2*(cy*(tmp1)-cx*tmp2))+3*(-3*cy*px1*Power(px3,2)*Power(py1,2)-6*Power(bx,2)*px3*Power(py1,3)+3*cx*Power(px3,2)*Power(py1,3)+9*cy*px1*px2*px3*py1*py2+9*Power(bx,2)*px2*Power(py1,2)*py2+9*Power(bx,2)*px3*Power(py1,2)*py2-9*cx*px2*px3*Power(py1,2)*py2-9*cy*Power(px1,2)*px3*Power(py2,2)-9*Power(bx,2)*px1*py1*Power(py2,2)-9*Power(bx,2)*px3*py1*Power(py2,2)+9*cx*px1*px3*py1*Power(py2,2)+3*Power(bx,2)*px3*Power(py2,3)-9*cy*px1*Power(px2,2)*py1*py3+6*cy*Power(px1,2)*px3*py1*py3+6*Power(bx,2)*px1*Power(py1,2)*py3-18*Power(bx,2)*px2*Power(py1,2)*py3+9*cx*Power(px2,2)*Power(py1,2)*py3+6*Power(bx,2)*px3*Power(py1,2)*py3-6*cx*px1*px3*Power(py1,2)*py3+9*cy*Power(px1,2)*px2*py2*py3+9*Power(bx,2)*px1*py1*py2*py3+9*Power(bx,2)*px2*py1*py2*py3-9*cx*px1*px2*py1*py2*py3-3*Power(bx,2)*px3*py1*py2*py3-3*Power(bx,2)*px2*Power(py2,2)*py3-3*cy*Power(px1,3)*Power(py3,2)-6*Power(bx,2)*px1*py1*Power(py3,2)+3*cx*Power(px1,2)*py1*Power(py3,2)+2*Power(bx,2)*px2*py1*Power(py3,2)+Power(bx,2)*px1*py2*Power(py3,2)+Power(by,2)*(6*Power(px1,3)*py3-px2*(Power(px3,2)*py1-3*px2*px3*py2+3*Power(px2,2)*py3)+px1*(2*Power(px3,2)*(3*py1-py2)+9*Power(px2,2)*(py1+py3)+3*px2*px3*(-3*py1-3*py2+py3))-3*Power(px1,2)*(2*px3*(py1-3*py2+py3)+3*px2*(py2+py3)))+bx*by*(3*Power(px3,2)*py1*(-2*py1+py2)+Power(px2,2)*(-9*Power(py1,2)-9*py1*py3+6*py2*py3)+3*Power(px1,2)*(3*Power(py2,2)-3*py2*py3+2*py3*(-2*py1+py3))+px1*px3*(12*Power(py1,2)-27*py1*py2+py2*(9*py2+py3))+px2*(3*px1*(9*py1-3*py2-py3)*py3+px3*(9*Power(py1,2)+9*py1*py2-6*Power(py2,2)-py1*py3))))+ax*(27*Power(bx,2)*Power(py1,3)-36*cx*px3*Power(py1,3)-81*Power(bx,2)*Power(py1,2)*py2+54*cx*px2*Power(py1,2)*py2+54*cx*px3*Power(py1,2)*py2+81*Power(bx,2)*py1*Power(py2,2)-54*cx*px1*py1*Power(py2,2)-54*cx*px3*py1*Power(py2,2)-27*Power(bx,2)*Power(py2,3)+18*cx*px3*Power(py2,3)+27*Power(bx,2)*Power(py1,2)*py3+36*cx*px1*Power(py1,2)*py3-108*cx*px2*Power(py1,2)*py3+36*cx*px3*Power(py1,2)*py3-54*Power(bx,2)*py1*py2*py3+54*cx*px1*py1*py2*py3+54*cx*px2*py1*py2*py3-18*cx*px3*py1*py2*py3+27*Power(bx,2)*Power(py2,2)*py3-18*cx*px2*Power(py2,2)*py3+9*Power(bx,2)*py1*Power(py3,2)-36*cx*px1*py1*Power(py3,2)+12*cx*px2*py1*Power(py3,2)-9*Power(bx,2)*py2*Power(py3,2)+6*cx*px1*py2*Power(py3,2)+Power(bx,2)*Power(py3,3)+Power(by,2)*Power(tmp1,2)*tmp2-2*bx*by*(tmp1)*Power(3*py1-3*py2+py3,2)-3*cy*(3*Power(px3,2)*py1*(2*py1-py2)+Power(px2,2)*(9*Power(py1,2)+9*py1*py3-6*py2*py3)+Power(px1,2)*(-9*Power(py2,2)+9*py2*py3+6*(2*py1-py3)*py3)-px1*px3*(12*Power(py1,2)-27*py1*py2+py2*(9*py2+py3))+px2*(3*px1*py3*(-9*py1+3*py2+py3)+px3*(-9*Power(py1,2)-9*py1*py2+6*Power(py2,2)+py1*py3)))));T c3=+(Power(by,3)*Power(tmp1,3)+162*ax*bx*cy*px1*Power(py1,2)+81*Power(ax,2)*dy*px1*Power(py1,2)-162*ax*bx*cy*px2*Power(py1,2)-81*Power(ax,2)*dy*px2*Power(py1,2)+81*bx*cy*Power(px2,2)*Power(py1,2)+81*ax*dy*Power(px2,2)*Power(py1,2)+54*ax*bx*cy*px3*Power(py1,2)+27*Power(ax,2)*dy*px3*Power(py1,2)-108*bx*cy*px1*px3*Power(py1,2)-108*ax*dy*px1*px3*Power(py1,2)-81*bx*cy*px2*px3*Power(py1,2)-81*ax*dy*px2*px3*Power(py1,2)+54*bx*cy*Power(px3,2)*Power(py1,2)+54*ax*dy*Power(px3,2)*Power(py1,2)+27*dy*px1*Power(px3,2)*Power(py1,2)-27*Power(bx,3)*Power(py1,3)-162*ax*bx*cx*Power(py1,3)-81*Power(ax,2)*dx*Power(py1,3)+108*bx*cx*px3*Power(py1,3)+108*ax*dx*px3*Power(py1,3)-27*dx*Power(px3,2)*Power(py1,3)-324*ax*bx*cy*px1*py1*py2-162*Power(ax,2)*dy*px1*py1*py2+324*ax*bx*cy*px2*py1*py2+162*Power(ax,2)*dy*px2*py1*py2-108*ax*bx*cy*px3*py1*py2-54*Power(ax,2)*dy*px3*py1*py2+243*bx*cy*px1*px3*py1*py2+243*ax*dy*px1*px3*py1*py2-81*bx*cy*px2*px3*py1*py2-81*ax*dy*px2*px3*py1*py2-81*dy*px1*px2*px3*py1*py2-27*bx*cy*Power(px3,2)*py1*py2-27*ax*dy*Power(px3,2)*py1*py2+81*Power(bx,3)*Power(py1,2)*py2+486*ax*bx*cx*Power(py1,2)*py2+243*Power(ax,2)*dx*Power(py1,2)*py2-162*bx*cx*px2*Power(py1,2)*py2-162*ax*dx*px2*Power(py1,2)*py2-162*bx*cx*px3*Power(py1,2)*py2-162*ax*dx*px3*Power(py1,2)*py2+81*dx*px2*px3*Power(py1,2)*py2+162*ax*bx*cy*px1*Power(py2,2)+81*Power(ax,2)*dy*px1*Power(py2,2)-81*bx*cy*Power(px1,2)*Power(py2,2)-81*ax*dy*Power(px1,2)*Power(py2,2)-162*ax*bx*cy*px2*Power(py2,2)-81*Power(ax,2)*dy*px2*Power(py2,2)+54*ax*bx*cy*px3*Power(py2,2)+27*Power(ax,2)*dy*px3*Power(py2,2)-81*bx*cy*px1*px3*Power(py2,2)-81*ax*dy*px1*px3*Power(py2,2)+81*dy*Power(px1,2)*px3*Power(py2,2)+54*bx*cy*px2*px3*Power(py2,2)+54*ax*dy*px2*px3*Power(py2,2)-81*Power(bx,3)*py1*Power(py2,2)-486*ax*bx*cx*py1*Power(py2,2)-243*Power(ax,2)*dx*py1*Power(py2,2)+162*bx*cx*px1*py1*Power(py2,2)+162*ax*dx*px1*py1*Power(py2,2)+162*bx*cx*px3*py1*Power(py2,2)+162*ax*dx*px3*py1*Power(py2,2)-81*dx*px1*px3*py1*Power(py2,2)+27*Power(bx,3)*Power(py2,3)+162*ax*bx*cx*Power(py2,3)+81*Power(ax,2)*dx*Power(py2,3)-54*bx*cx*px3*Power(py2,3)-54*ax*dx*px3*Power(py2,3)+108*ax*bx*cy*px1*py1*py3+54*Power(ax,2)*dy*px1*py1*py3+108*bx*cy*Power(px1,2)*py1*py3+108*ax*dy*Power(px1,2)*py1*py3-108*ax*bx*cy*px2*py1*py3-54*Power(ax,2)*dy*px2*py1*py3-243*bx*cy*px1*px2*py1*py3-243*ax*dy*px1*px2*py1*py3+81*bx*cy*Power(px2,2)*py1*py3+81*ax*dy*Power(px2,2)*py1*py3+81*dy*px1*Power(px2,2)*py1*py3+36*ax*bx*cy*px3*py1*py3+18*Power(ax,2)*dy*px3*py1*py3-54*dy*Power(px1,2)*px3*py1*py3+9*bx*cy*px2*px3*py1*py3+9*ax*dy*px2*px3*py1*py3-27*Power(bx,3)*Power(py1,2)*py3-162*ax*bx*cx*Power(py1,2)*py3-81*Power(ax,2)*dx*Power(py1,2)*py3-108*bx*cx*px1*Power(py1,2)*py3-108*ax*dx*px1*Power(py1,2)*py3+324*bx*cx*px2*Power(py1,2)*py3+324*ax*dx*px2*Power(py1,2)*py3-81*dx*Power(px2,2)*Power(py1,2)*py3-108*bx*cx*px3*Power(py1,2)*py3-108*ax*dx*px3*Power(py1,2)*py3+54*dx*px1*px3*Power(py1,2)*py3-108*ax*bx*cy*px1*py2*py3-54*Power(ax,2)*dy*px1*py2*py3+81*bx*cy*Power(px1,2)*py2*py3+81*ax*dy*Power(px1,2)*py2*py3+108*ax*bx*cy*px2*py2*py3+54*Power(ax,2)*dy*px2*py2*py3+81*bx*cy*px1*px2*py2*py3+81*ax*dy*px1*px2*py2*py3-81*dy*Power(px1,2)*px2*py2*py3-54*bx*cy*Power(px2,2)*py2*py3-54*ax*dy*Power(px2,2)*py2*py3-36*ax*bx*cy*px3*py2*py3-18*Power(ax,2)*dy*px3*py2*py3-9*bx*cy*px1*px3*py2*py3-9*ax*dy*px1*px3*py2*py3+54*Power(bx,3)*py1*py2*py3+324*ax*bx*cx*py1*py2*py3+162*Power(ax,2)*dx*py1*py2*py3-162*bx*cx*px1*py1*py2*py3-162*ax*dx*px1*py1*py2*py3-162*bx*cx*px2*py1*py2*py3-162*ax*dx*px2*py1*py2*py3+81*dx*px1*px2*py1*py2*py3+54*bx*cx*px3*py1*py2*py3+54*ax*dx*px3*py1*py2*py3-27*Power(bx,3)*Power(py2,2)*py3-162*ax*bx*cx*Power(py2,2)*py3-81*Power(ax,2)*dx*Power(py2,2)*py3+54*bx*cx*px2*Power(py2,2)*py3+54*ax*dx*px2*Power(py2,2)*py3+18*ax*bx*cy*px1*Power(py3,2)+9*Power(ax,2)*dy*px1*Power(py3,2)-54*bx*cy*Power(px1,2)*Power(py3,2)-54*ax*dy*Power(px1,2)*Power(py3,2)+27*dy*Power(px1,3)*Power(py3,2)-18*ax*bx*cy*px2*Power(py3,2)-9*Power(ax,2)*dy*px2*Power(py3,2)+27*bx*cy*px1*px2*Power(py3,2)+27*ax*dy*px1*px2*Power(py3,2)+6*ax*bx*cy*px3*Power(py3,2)+3*Power(ax,2)*dy*px3*Power(py3,2)-9*Power(bx,3)*py1*Power(py3,2)-54*ax*bx*cx*py1*Power(py3,2)-27*Power(ax,2)*dx*py1*Power(py3,2)+108*bx*cx*px1*py1*Power(py3,2)+108*ax*dx*px1*py1*Power(py3,2)-27*dx*Power(px1,2)*py1*Power(py3,2)-36*bx*cx*px2*py1*Power(py3,2)-36*ax*dx*px2*py1*Power(py3,2)+9*Power(bx,3)*py2*Power(py3,2)+54*ax*bx*cx*py2*Power(py3,2)+27*Power(ax,2)*dx*py2*Power(py3,2)-18*bx*cx*px1*py2*Power(py3,2)-18*ax*dx*px1*py2*Power(py3,2)-Power(bx,3)*Power(py3,3)-6*ax*bx*cx*Power(py3,3)-3*Power(ax,2)*dx*Power(py3,3)-3*bx*Power(by,2)*Power(tmp1,2)*tmp2+3*Power(ay,2)*Power(tmp1,2)*(dy*(tmp1)-dx*tmp2)+3*by*(-54*cy*px1*Power(px2,2)*py1+36*cy*Power(px1,2)*px3*py1+54*cy*px1*px2*px3*py1-36*cy*px1*Power(px3,2)*py1+6*cy*px2*Power(px3,2)*py1+27*Power(bx,2)*px1*Power(py1,2)-27*Power(bx,2)*px2*Power(py1,2)+27*cx*Power(px2,2)*Power(py1,2)+9*Power(bx,2)*px3*Power(py1,2)-36*cx*px1*px3*Power(py1,2)-27*cx*px2*px3*Power(py1,2)+18*cx*Power(px3,2)*Power(py1,2)+54*cy*Power(px1,2)*px2*py2-108*cy*Power(px1,2)*px3*py2+54*cy*px1*px2*px3*py2-18*cy*Power(px2,2)*px3*py2+12*cy*px1*Power(px3,2)*py2-54*Power(bx,2)*px1*py1*py2+54*Power(bx,2)*px2*py1*py2-18*Power(bx,2)*px3*py1*py2+81*cx*px1*px3*py1*py2-27*cx*px2*px3*py1*py2-9*cx*Power(px3,2)*py1*py2+27*Power(bx,2)*px1*Power(py2,2)-27*cx*Power(px1,2)*Power(py2,2)-27*Power(bx,2)*px2*Power(py2,2)+9*Power(bx,2)*px3*Power(py2,2)-27*cx*px1*px3*Power(py2,2)+18*cx*px2*px3*Power(py2,2)-36*cy*Power(px1,3)*py3+54*cy*Power(px1,2)*px2*py3-54*cy*px1*Power(px2,2)*py3+18*cy*Power(px2,3)*py3+36*cy*Power(px1,2)*px3*py3-18*cy*px1*px2*px3*py3+18*Power(bx,2)*px1*py1*py3+36*cx*Power(px1,2)*py1*py3-18*Power(bx,2)*px2*py1*py3-81*cx*px1*px2*py1*py3+27*cx*Power(px2,2)*py1*py3+6*Power(bx,2)*px3*py1*py3+3*cx*px2*px3*py1*py3-18*Power(bx,2)*px1*py2*py3+27*cx*Power(px1,2)*py2*py3+18*Power(bx,2)*px2*py2*py3+27*cx*px1*px2*py2*py3-18*cx*Power(px2,2)*py2*py3-6*Power(bx,2)*px3*py2*py3-3*cx*px1*px3*py2*py3+3*Power(bx,2)*px1*Power(py3,2)-18*cx*Power(px1,2)*Power(py3,2)-3*Power(bx,2)*px2*Power(py3,2)+9*cx*px1*px2*Power(py3,2)+Power(bx,2)*px3*Power(py3,2)+2*ay*Power(tmp1,2)*(cy*(tmp1)-cx*tmp2)-2*ax*(tmp1)*tmp2*(cy*(tmp1)-cx*tmp2))-3*ay*(2*bx*(tmp1)*tmp2*(cy*(tmp1)-cx*tmp2)+2*ax*(tmp1)*tmp2*(dy*(tmp1)-dx*tmp2)+6*dy*(6*Power(px1,3)*py3-px2*(Power(px3,2)*py1-3*px2*px3*py2+3*Power(px2,2)*py3)+px1*(2*Power(px3,2)*(3*py1-py2)+9*Power(px2,2)*(py1+py3)+3*px2*px3*(-3*py1-3*py2+py3))-3*Power(px1,2)*(2*px3*(py1-3*py2+py3)+3*px2*(py2+py3)))-3*dx*(3*Power(px3,2)*py1*(2*py1-py2)+Power(px2,2)*(9*Power(py1,2)+9*py1*py3-6*py2*py3)+Power(px1,2)*(-9*Power(py2,2)+9*py2*py3+6*(2*py1-py3)*py3)-px1*px3*(12*Power(py1,2)-27*py1*py2+py2*(9*py2+py3))+px2*(3*px1*py3*(-9*py1+3*py2+py3)+px3*(-9*Power(py1,2)-9*py1*py2+6*Power(py2,2)+py1*py3)))));T c4=-3*(27*ax*Power(cy,2)*Power(px1,2)*py1-54*ax*Power(cy,2)*px1*px2*py1+27*ax*Power(cy,2)*Power(px2,2)*py1+27*Power(cy,2)*px1*Power(px2,2)*py1+18*ax*Power(cy,2)*px1*px3*py1-18*Power(cy,2)*Power(px1,2)*px3*py1-18*ax*Power(cy,2)*px2*px3*py1-27*Power(cy,2)*px1*px2*px3*py1+3*ax*Power(cy,2)*Power(px3,2)*py1+18*Power(cy,2)*px1*Power(px3,2)*py1-3*Power(cy,2)*px2*Power(px3,2)*py1-27*Power(bx,2)*cy*px1*Power(py1,2)-54*ax*cx*cy*px1*Power(py1,2)-54*ax*bx*dy*px1*Power(py1,2)+27*Power(bx,2)*cy*px2*Power(py1,2)+54*ax*cx*cy*px2*Power(py1,2)+54*ax*bx*dy*px2*Power(py1,2)-27*cx*cy*Power(px2,2)*Power(py1,2)-27*bx*dy*Power(px2,2)*Power(py1,2)-9*Power(bx,2)*cy*px3*Power(py1,2)-18*ax*cx*cy*px3*Power(py1,2)-18*ax*bx*dy*px3*Power(py1,2)+36*cx*cy*px1*px3*Power(py1,2)+36*bx*dy*px1*px3*Power(py1,2)+27*cx*cy*px2*px3*Power(py1,2)+27*bx*dy*px2*px3*Power(py1,2)-18*cx*cy*Power(px3,2)*Power(py1,2)-18*bx*dy*Power(px3,2)*Power(py1,2)+27*Power(bx,2)*cx*Power(py1,3)+27*ax*Power(cx,2)*Power(py1,3)+54*ax*bx*dx*Power(py1,3)-18*Power(cx,2)*px3*Power(py1,3)-36*bx*dx*px3*Power(py1,3)-27*ax*Power(cy,2)*Power(px1,2)*py2+54*ax*Power(cy,2)*px1*px2*py2-27*Power(cy,2)*Power(px1,2)*px2*py2-27*ax*Power(cy,2)*Power(px2,2)*py2-18*ax*Power(cy,2)*px1*px3*py2+54*Power(cy,2)*Power(px1,2)*px3*py2+18*ax*Power(cy,2)*px2*px3*py2-27*Power(cy,2)*px1*px2*px3*py2+9*Power(cy,2)*Power(px2,2)*px3*py2-3*ax*Power(cy,2)*Power(px3,2)*py2-6*Power(cy,2)*px1*Power(px3,2)*py2+54*Power(bx,2)*cy*px1*py1*py2+108*ax*cx*cy*px1*py1*py2+108*ax*bx*dy*px1*py1*py2-54*Power(bx,2)*cy*px2*py1*py2-108*ax*cx*cy*px2*py1*py2-108*ax*bx*dy*px2*py1*py2+18*Power(bx,2)*cy*px3*py1*py2+36*ax*cx*cy*px3*py1*py2+36*ax*bx*dy*px3*py1*py2-81*cx*cy*px1*px3*py1*py2-81*bx*dy*px1*px3*py1*py2+27*cx*cy*px2*px3*py1*py2+27*bx*dy*px2*px3*py1*py2+9*cx*cy*Power(px3,2)*py1*py2+9*bx*dy*Power(px3,2)*py1*py2-81*Power(bx,2)*cx*Power(py1,2)*py2-81*ax*Power(cx,2)*Power(py1,2)*py2-162*ax*bx*dx*Power(py1,2)*py2+27*Power(cx,2)*px2*Power(py1,2)*py2+54*bx*dx*px2*Power(py1,2)*py2+27*Power(cx,2)*px3*Power(py1,2)*py2+54*bx*dx*px3*Power(py1,2)*py2-27*Power(bx,2)*cy*px1*Power(py2,2)-54*ax*cx*cy*px1*Power(py2,2)-54*ax*bx*dy*px1*Power(py2,2)+27*cx*cy*Power(px1,2)*Power(py2,2)+27*bx*dy*Power(px1,2)*Power(py2,2)+27*Power(bx,2)*cy*px2*Power(py2,2)+54*ax*cx*cy*px2*Power(py2,2)+54*ax*bx*dy*px2*Power(py2,2)-9*Power(bx,2)*cy*px3*Power(py2,2)-18*ax*cx*cy*px3*Power(py2,2)-18*ax*bx*dy*px3*Power(py2,2)+27*cx*cy*px1*px3*Power(py2,2)+27*bx*dy*px1*px3*Power(py2,2)-18*cx*cy*px2*px3*Power(py2,2)-18*bx*dy*px2*px3*Power(py2,2)+81*Power(bx,2)*cx*py1*Power(py2,2)+81*ax*Power(cx,2)*py1*Power(py2,2)+162*ax*bx*dx*py1*Power(py2,2)-27*Power(cx,2)*px1*py1*Power(py2,2)-54*bx*dx*px1*py1*Power(py2,2)-27*Power(cx,2)*px3*py1*Power(py2,2)-54*bx*dx*px3*py1*Power(py2,2)-27*Power(bx,2)*cx*Power(py2,3)-27*ax*Power(cx,2)*Power(py2,3)-54*ax*bx*dx*Power(py2,3)+9*Power(cx,2)*px3*Power(py2,3)+18*bx*dx*px3*Power(py2,3)+9*ax*Power(cy,2)*Power(px1,2)*py3+18*Power(cy,2)*Power(px1,3)*py3-18*ax*Power(cy,2)*px1*px2*py3-27*Power(cy,2)*Power(px1,2)*px2*py3+9*ax*Power(cy,2)*Power(px2,2)*py3+27*Power(cy,2)*px1*Power(px2,2)*py3-9*Power(cy,2)*Power(px2,3)*py3+6*ax*Power(cy,2)*px1*px3*py3-18*Power(cy,2)*Power(px1,2)*px3*py3-6*ax*Power(cy,2)*px2*px3*py3+9*Power(cy,2)*px1*px2*px3*py3+ax*Power(cy,2)*Power(px3,2)*py3-18*Power(bx,2)*cy*px1*py1*py3-36*ax*cx*cy*px1*py1*py3-36*ax*bx*dy*px1*py1*py3-36*cx*cy*Power(px1,2)*py1*py3-36*bx*dy*Power(px1,2)*py1*py3+18*Power(bx,2)*cy*px2*py1*py3+36*ax*cx*cy*px2*py1*py3+36*ax*bx*dy*px2*py1*py3+81*cx*cy*px1*px2*py1*py3+81*bx*dy*px1*px2*py1*py3-27*cx*cy*Power(px2,2)*py1*py3-27*bx*dy*Power(px2,2)*py1*py3-6*Power(bx,2)*cy*px3*py1*py3-12*ax*cx*cy*px3*py1*py3-12*ax*bx*dy*px3*py1*py3-3*cx*cy*px2*px3*py1*py3-3*bx*dy*px2*px3*py1*py3+27*Power(bx,2)*cx*Power(py1,2)*py3+27*ax*Power(cx,2)*Power(py1,2)*py3+54*ax*bx*dx*Power(py1,2)*py3+18*Power(cx,2)*px1*Power(py1,2)*py3+36*bx*dx*px1*Power(py1,2)*py3-54*Power(cx,2)*px2*Power(py1,2)*py3-108*bx*dx*px2*Power(py1,2)*py3+18*Power(cx,2)*px3*Power(py1,2)*py3+36*bx*dx*px3*Power(py1,2)*py3+18*Power(bx,2)*cy*px1*py2*py3+36*ax*cx*cy*px1*py2*py3+36*ax*bx*dy*px1*py2*py3-27*cx*cy*Power(px1,2)*py2*py3-27*bx*dy*Power(px1,2)*py2*py3-18*Power(bx,2)*cy*px2*py2*py3-36*ax*cx*cy*px2*py2*py3-36*ax*bx*dy*px2*py2*py3-27*cx*cy*px1*px2*py2*py3-27*bx*dy*px1*px2*py2*py3+18*cx*cy*Power(px2,2)*py2*py3+18*bx*dy*Power(px2,2)*py2*py3+6*Power(bx,2)*cy*px3*py2*py3+12*ax*cx*cy*px3*py2*py3+12*ax*bx*dy*px3*py2*py3+3*cx*cy*px1*px3*py2*py3+3*bx*dy*px1*px3*py2*py3-54*Power(bx,2)*cx*py1*py2*py3-54*ax*Power(cx,2)*py1*py2*py3-108*ax*bx*dx*py1*py2*py3+27*Power(cx,2)*px1*py1*py2*py3+54*bx*dx*px1*py1*py2*py3+27*Power(cx,2)*px2*py1*py2*py3+54*bx*dx*px2*py1*py2*py3-9*Power(cx,2)*px3*py1*py2*py3-18*bx*dx*px3*py1*py2*py3+27*Power(bx,2)*cx*Power(py2,2)*py3+27*ax*Power(cx,2)*Power(py2,2)*py3+54*ax*bx*dx*Power(py2,2)*py3-9*Power(cx,2)*px2*Power(py2,2)*py3-18*bx*dx*px2*Power(py2,2)*py3-3*Power(bx,2)*cy*px1*Power(py3,2)-6*ax*cx*cy*px1*Power(py3,2)-6*ax*bx*dy*px1*Power(py3,2)+18*cx*cy*Power(px1,2)*Power(py3,2)+18*bx*dy*Power(px1,2)*Power(py3,2)+3*Power(bx,2)*cy*px2*Power(py3,2)+6*ax*cx*cy*px2*Power(py3,2)+6*ax*bx*dy*px2*Power(py3,2)-9*cx*cy*px1*px2*Power(py3,2)-9*bx*dy*px1*px2*Power(py3,2)-Power(bx,2)*cy*px3*Power(py3,2)-2*ax*cx*cy*px3*Power(py3,2)-2*ax*bx*dy*px3*Power(py3,2)+9*Power(bx,2)*cx*py1*Power(py3,2)+9*ax*Power(cx,2)*py1*Power(py3,2)+18*ax*bx*dx*py1*Power(py3,2)-18*Power(cx,2)*px1*py1*Power(py3,2)-36*bx*dx*px1*py1*Power(py3,2)+6*Power(cx,2)*px2*py1*Power(py3,2)+12*bx*dx*px2*py1*Power(py3,2)-9*Power(bx,2)*cx*py2*Power(py3,2)-9*ax*Power(cx,2)*py2*Power(py3,2)-18*ax*bx*dx*py2*Power(py3,2)+3*Power(cx,2)*px1*py2*Power(py3,2)+6*bx*dx*px1*py2*Power(py3,2)+Power(bx,2)*cx*Power(py3,3)+ax*Power(cx,2)*Power(py3,3)+2*ax*bx*dx*Power(py3,3)+Power(by,2)*Power(tmp1,2)*(-(cy*(tmp1))+cx*tmp2)-by*(-54*ax*dy*Power(px1,2)*py1+108*ax*dy*px1*px2*py1-54*ax*dy*Power(px2,2)*py1-54*dy*px1*Power(px2,2)*py1-36*ax*dy*px1*px3*py1+36*dy*Power(px1,2)*px3*py1+36*ax*dy*px2*px3*py1+54*dy*px1*px2*px3*py1-6*ax*dy*Power(px3,2)*py1-36*dy*px1*Power(px3,2)*py1+6*dy*px2*Power(px3,2)*py1+54*ax*dx*px1*Power(py1,2)-54*ax*dx*px2*Power(py1,2)+27*dx*Power(px2,2)*Power(py1,2)+18*ax*dx*px3*Power(py1,2)-36*dx*px1*px3*Power(py1,2)-27*dx*px2*px3*Power(py1,2)+18*dx*Power(px3,2)*Power(py1,2)+54*ax*dy*Power(px1,2)*py2-108*ax*dy*px1*px2*py2+54*dy*Power(px1,2)*px2*py2+54*ax*dy*Power(px2,2)*py2+36*ax*dy*px1*px3*py2-108*dy*Power(px1,2)*px3*py2-36*ax*dy*px2*px3*py2+54*dy*px1*px2*px3*py2-18*dy*Power(px2,2)*px3*py2+6*ax*dy*Power(px3,2)*py2+12*dy*px1*Power(px3,2)*py2-108*ax*dx*px1*py1*py2+108*ax*dx*px2*py1*py2-36*ax*dx*px3*py1*py2+81*dx*px1*px3*py1*py2-27*dx*px2*px3*py1*py2-9*dx*Power(px3,2)*py1*py2+54*ax*dx*px1*Power(py2,2)-27*dx*Power(px1,2)*Power(py2,2)-54*ax*dx*px2*Power(py2,2)+18*ax*dx*px3*Power(py2,2)-27*dx*px1*px3*Power(py2,2)+18*dx*px2*px3*Power(py2,2)-18*ax*dy*Power(px1,2)*py3-36*dy*Power(px1,3)*py3+36*ax*dy*px1*px2*py3+54*dy*Power(px1,2)*px2*py3-18*ax*dy*Power(px2,2)*py3-54*dy*px1*Power(px2,2)*py3+18*dy*Power(px2,3)*py3-12*ax*dy*px1*px3*py3+36*dy*Power(px1,2)*px3*py3+12*ax*dy*px2*px3*py3-18*dy*px1*px2*px3*py3-2*ax*dy*Power(px3,2)*py3+36*ax*dx*px1*py1*py3+36*dx*Power(px1,2)*py1*py3-36*ax*dx*px2*py1*py3-81*dx*px1*px2*py1*py3+27*dx*Power(px2,2)*py1*py3+12*ax*dx*px3*py1*py3+3*dx*px2*px3*py1*py3-36*ax*dx*px1*py2*py3+27*dx*Power(px1,2)*py2*py3+36*ax*dx*px2*py2*py3+27*dx*px1*px2*py2*py3-18*dx*Power(px2,2)*py2*py3-12*ax*dx*px3*py2*py3-3*dx*px1*px3*py2*py3+6*ax*dx*px1*Power(py3,2)-18*dx*Power(px1,2)*Power(py3,2)-6*ax*dx*px2*Power(py3,2)+9*dx*px1*px2*Power(py3,2)+2*ax*dx*px3*Power(py3,2)-2*bx*(tmp1)*tmp2*(cy*(tmp1)-cx*tmp2)+2*ay*Power(tmp1,2)*(dy*(tmp1)-dx*tmp2))-ay*(tmp1)*(Power(cy,2)*Power(tmp1,2)-2*cx*cy*(tmp1)*tmp2+tmp2*(-2*bx*dy*(tmp1)+Power(cx,2)*tmp2+2*bx*dx*tmp2)));T c5=-3*(27*bx*Power(cy,2)*Power(px1,2)*py1+54*ax*cy*dy*Power(px1,2)*py1-54*bx*Power(cy,2)*px1*px2*py1-108*ax*cy*dy*px1*px2*py1+27*bx*Power(cy,2)*Power(px2,2)*py1+54*ax*cy*dy*Power(px2,2)*py1+54*cy*dy*px1*Power(px2,2)*py1+18*bx*Power(cy,2)*px1*px3*py1+36*ax*cy*dy*px1*px3*py1-36*cy*dy*Power(px1,2)*px3*py1-18*bx*Power(cy,2)*px2*px3*py1-36*ax*cy*dy*px2*px3*py1-54*cy*dy*px1*px2*px3*py1+3*bx*Power(cy,2)*Power(px3,2)*py1+6*ax*cy*dy*Power(px3,2)*py1+36*cy*dy*px1*Power(px3,2)*py1-6*cy*dy*px2*Power(px3,2)*py1-54*bx*cx*cy*px1*Power(py1,2)-54*ax*cy*dx*px1*Power(py1,2)-27*Power(bx,2)*dy*px1*Power(py1,2)-54*ax*cx*dy*px1*Power(py1,2)+54*bx*cx*cy*px2*Power(py1,2)+54*ax*cy*dx*px2*Power(py1,2)+27*Power(bx,2)*dy*px2*Power(py1,2)+54*ax*cx*dy*px2*Power(py1,2)-27*cy*dx*Power(px2,2)*Power(py1,2)-27*cx*dy*Power(px2,2)*Power(py1,2)-18*bx*cx*cy*px3*Power(py1,2)-18*ax*cy*dx*px3*Power(py1,2)-9*Power(bx,2)*dy*px3*Power(py1,2)-18*ax*cx*dy*px3*Power(py1,2)+36*cy*dx*px1*px3*Power(py1,2)+36*cx*dy*px1*px3*Power(py1,2)+27*cy*dx*px2*px3*Power(py1,2)+27*cx*dy*px2*px3*Power(py1,2)-18*cy*dx*Power(px3,2)*Power(py1,2)-18*cx*dy*Power(px3,2)*Power(py1,2)+27*bx*Power(cx,2)*Power(py1,3)+27*Power(bx,2)*dx*Power(py1,3)+54*ax*cx*dx*Power(py1,3)-36*cx*dx*px3*Power(py1,3)-27*bx*Power(cy,2)*Power(px1,2)*py2-54*ax*cy*dy*Power(px1,2)*py2+54*bx*Power(cy,2)*px1*px2*py2+108*ax*cy*dy*px1*px2*py2-54*cy*dy*Power(px1,2)*px2*py2-27*bx*Power(cy,2)*Power(px2,2)*py2-54*ax*cy*dy*Power(px2,2)*py2-18*bx*Power(cy,2)*px1*px3*py2-36*ax*cy*dy*px1*px3*py2+108*cy*dy*Power(px1,2)*px3*py2+18*bx*Power(cy,2)*px2*px3*py2+36*ax*cy*dy*px2*px3*py2-54*cy*dy*px1*px2*px3*py2+18*cy*dy*Power(px2,2)*px3*py2-3*bx*Power(cy,2)*Power(px3,2)*py2-6*ax*cy*dy*Power(px3,2)*py2-12*cy*dy*px1*Power(px3,2)*py2+108*bx*cx*cy*px1*py1*py2+108*ax*cy*dx*px1*py1*py2+54*Power(bx,2)*dy*px1*py1*py2+108*ax*cx*dy*px1*py1*py2-108*bx*cx*cy*px2*py1*py2-108*ax*cy*dx*px2*py1*py2-54*Power(bx,2)*dy*px2*py1*py2-108*ax*cx*dy*px2*py1*py2+36*bx*cx*cy*px3*py1*py2+36*ax*cy*dx*px3*py1*py2+18*Power(bx,2)*dy*px3*py1*py2+36*ax*cx*dy*px3*py1*py2-81*cy*dx*px1*px3*py1*py2-81*cx*dy*px1*px3*py1*py2+27*cy*dx*px2*px3*py1*py2+27*cx*dy*px2*px3*py1*py2+9*cy*dx*Power(px3,2)*py1*py2+9*cx*dy*Power(px3,2)*py1*py2-81*bx*Power(cx,2)*Power(py1,2)*py2-81*Power(bx,2)*dx*Power(py1,2)*py2-162*ax*cx*dx*Power(py1,2)*py2+54*cx*dx*px2*Power(py1,2)*py2+54*cx*dx*px3*Power(py1,2)*py2-54*bx*cx*cy*px1*Power(py2,2)-54*ax*cy*dx*px1*Power(py2,2)-27*Power(bx,2)*dy*px1*Power(py2,2)-54*ax*cx*dy*px1*Power(py2,2)+27*cy*dx*Power(px1,2)*Power(py2,2)+27*cx*dy*Power(px1,2)*Power(py2,2)+54*bx*cx*cy*px2*Power(py2,2)+54*ax*cy*dx*px2*Power(py2,2)+27*Power(bx,2)*dy*px2*Power(py2,2)+54*ax*cx*dy*px2*Power(py2,2)-18*bx*cx*cy*px3*Power(py2,2)-18*ax*cy*dx*px3*Power(py2,2)-9*Power(bx,2)*dy*px3*Power(py2,2)-18*ax*cx*dy*px3*Power(py2,2)+27*cy*dx*px1*px3*Power(py2,2)+27*cx*dy*px1*px3*Power(py2,2)-18*cy*dx*px2*px3*Power(py2,2)-18*cx*dy*px2*px3*Power(py2,2)+81*bx*Power(cx,2)*py1*Power(py2,2)+81*Power(bx,2)*dx*py1*Power(py2,2)+162*ax*cx*dx*py1*Power(py2,2)-54*cx*dx*px1*py1*Power(py2,2)-54*cx*dx*px3*py1*Power(py2,2)-27*bx*Power(cx,2)*Power(py2,3)-27*Power(bx,2)*dx*Power(py2,3)-54*ax*cx*dx*Power(py2,3)+18*cx*dx*px3*Power(py2,3)+9*bx*Power(cy,2)*Power(px1,2)*py3+18*ax*cy*dy*Power(px1,2)*py3+36*cy*dy*Power(px1,3)*py3-18*bx*Power(cy,2)*px1*px2*py3-36*ax*cy*dy*px1*px2*py3-54*cy*dy*Power(px1,2)*px2*py3+9*bx*Power(cy,2)*Power(px2,2)*py3+18*ax*cy*dy*Power(px2,2)*py3+54*cy*dy*px1*Power(px2,2)*py3-18*cy*dy*Power(px2,3)*py3+6*bx*Power(cy,2)*px1*px3*py3+12*ax*cy*dy*px1*px3*py3-36*cy*dy*Power(px1,2)*px3*py3-6*bx*Power(cy,2)*px2*px3*py3-12*ax*cy*dy*px2*px3*py3+18*cy*dy*px1*px2*px3*py3+bx*Power(cy,2)*Power(px3,2)*py3+2*ax*cy*dy*Power(px3,2)*py3-36*bx*cx*cy*px1*py1*py3-36*ax*cy*dx*px1*py1*py3-18*Power(bx,2)*dy*px1*py1*py3-36*ax*cx*dy*px1*py1*py3-36*cy*dx*Power(px1,2)*py1*py3-36*cx*dy*Power(px1,2)*py1*py3+36*bx*cx*cy*px2*py1*py3+36*ax*cy*dx*px2*py1*py3+18*Power(bx,2)*dy*px2*py1*py3+36*ax*cx*dy*px2*py1*py3+81*cy*dx*px1*px2*py1*py3+81*cx*dy*px1*px2*py1*py3-27*cy*dx*Power(px2,2)*py1*py3-27*cx*dy*Power(px2,2)*py1*py3-12*bx*cx*cy*px3*py1*py3-12*ax*cy*dx*px3*py1*py3-6*Power(bx,2)*dy*px3*py1*py3-12*ax*cx*dy*px3*py1*py3-3*cy*dx*px2*px3*py1*py3-3*cx*dy*px2*px3*py1*py3+27*bx*Power(cx,2)*Power(py1,2)*py3+27*Power(bx,2)*dx*Power(py1,2)*py3+54*ax*cx*dx*Power(py1,2)*py3+36*cx*dx*px1*Power(py1,2)*py3-108*cx*dx*px2*Power(py1,2)*py3+36*cx*dx*px3*Power(py1,2)*py3+36*bx*cx*cy*px1*py2*py3+36*ax*cy*dx*px1*py2*py3+18*Power(bx,2)*dy*px1*py2*py3+36*ax*cx*dy*px1*py2*py3-27*cy*dx*Power(px1,2)*py2*py3-27*cx*dy*Power(px1,2)*py2*py3-36*bx*cx*cy*px2*py2*py3-36*ax*cy*dx*px2*py2*py3-18*Power(bx,2)*dy*px2*py2*py3-36*ax*cx*dy*px2*py2*py3-27*cy*dx*px1*px2*py2*py3-27*cx*dy*px1*px2*py2*py3+18*cy*dx*Power(px2,2)*py2*py3+18*cx*dy*Power(px2,2)*py2*py3+12*bx*cx*cy*px3*py2*py3+12*ax*cy*dx*px3*py2*py3+6*Power(bx,2)*dy*px3*py2*py3+12*ax*cx*dy*px3*py2*py3+3*cy*dx*px1*px3*py2*py3+3*cx*dy*px1*px3*py2*py3-54*bx*Power(cx,2)*py1*py2*py3-54*Power(bx,2)*dx*py1*py2*py3-108*ax*cx*dx*py1*py2*py3+54*cx*dx*px1*py1*py2*py3+54*cx*dx*px2*py1*py2*py3-18*cx*dx*px3*py1*py2*py3+27*bx*Power(cx,2)*Power(py2,2)*py3+27*Power(bx,2)*dx*Power(py2,2)*py3+54*ax*cx*dx*Power(py2,2)*py3-18*cx*dx*px2*Power(py2,2)*py3-6*bx*cx*cy*px1*Power(py3,2)-6*ax*cy*dx*px1*Power(py3,2)-3*Power(bx,2)*dy*px1*Power(py3,2)-6*ax*cx*dy*px1*Power(py3,2)+18*cy*dx*Power(px1,2)*Power(py3,2)+18*cx*dy*Power(px1,2)*Power(py3,2)+6*bx*cx*cy*px2*Power(py3,2)+6*ax*cy*dx*px2*Power(py3,2)+3*Power(bx,2)*dy*px2*Power(py3,2)+6*ax*cx*dy*px2*Power(py3,2)-9*cy*dx*px1*px2*Power(py3,2)-9*cx*dy*px1*px2*Power(py3,2)-2*bx*cx*cy*px3*Power(py3,2)-2*ax*cy*dx*px3*Power(py3,2)-Power(bx,2)*dy*px3*Power(py3,2)-2*ax*cx*dy*px3*Power(py3,2)+9*bx*Power(cx,2)*py1*Power(py3,2)+9*Power(bx,2)*dx*py1*Power(py3,2)+18*ax*cx*dx*py1*Power(py3,2)-36*cx*dx*px1*py1*Power(py3,2)+12*cx*dx*px2*py1*Power(py3,2)-9*bx*Power(cx,2)*py2*Power(py3,2)-9*Power(bx,2)*dx*py2*Power(py3,2)-18*ax*cx*dx*py2*Power(py3,2)+6*cx*dx*px1*py2*Power(py3,2)+bx*Power(cx,2)*Power(py3,3)+Power(bx,2)*dx*Power(py3,3)+2*ax*cx*dx*Power(py3,3)-2*ay*(tmp1)*(cy*(tmp1)-cx*tmp2)*(dy*(tmp1)-dx*tmp2)+Power(by,2)*Power(tmp1,2)*(-(dy*(tmp1))+dx*tmp2)-by*(tmp1)*(Power(cy,2)*Power(tmp1,2)-2*cx*cy*(tmp1)*tmp2+tmp2*(-2*bx*dy*(tmp1)+Power(cx,2)*tmp2+2*bx*dx*tmp2)));T c6=+(Power(cy,3)*Power(tmp1,3)-162*by*cx*dy*Power(px1,2)*py1-81*ax*Power(dy,2)*Power(px1,2)*py1+324*by*cx*dy*px1*px2*py1+162*ax*Power(dy,2)*px1*px2*py1-162*by*cx*dy*Power(px2,2)*py1-81*ax*Power(dy,2)*Power(px2,2)*py1-81*Power(dy,2)*px1*Power(px2,2)*py1-108*by*cx*dy*px1*px3*py1-54*ax*Power(dy,2)*px1*px3*py1+54*Power(dy,2)*Power(px1,2)*px3*py1+108*by*cx*dy*px2*px3*py1+54*ax*Power(dy,2)*px2*px3*py1+81*Power(dy,2)*px1*px2*px3*py1-18*by*cx*dy*Power(px3,2)*py1-9*ax*Power(dy,2)*Power(px3,2)*py1-54*Power(dy,2)*px1*Power(px3,2)*py1+9*Power(dy,2)*px2*Power(px3,2)*py1+162*by*cx*dx*px1*Power(py1,2)+162*bx*cx*dy*px1*Power(py1,2)+162*ax*dx*dy*px1*Power(py1,2)-162*by*cx*dx*px2*Power(py1,2)-162*bx*cx*dy*px2*Power(py1,2)-162*ax*dx*dy*px2*Power(py1,2)+81*dx*dy*Power(px2,2)*Power(py1,2)+54*by*cx*dx*px3*Power(py1,2)+54*bx*cx*dy*px3*Power(py1,2)+54*ax*dx*dy*px3*Power(py1,2)-108*dx*dy*px1*px3*Power(py1,2)-81*dx*dy*px2*px3*Power(py1,2)+54*dx*dy*Power(px3,2)*Power(py1,2)-27*Power(cx,3)*Power(py1,3)-162*bx*cx*dx*Power(py1,3)-81*ax*Power(dx,2)*Power(py1,3)+54*Power(dx,2)*px3*Power(py1,3)+162*by*cx*dy*Power(px1,2)*py2+81*ax*Power(dy,2)*Power(px1,2)*py2-324*by*cx*dy*px1*px2*py2-162*ax*Power(dy,2)*px1*px2*py2+81*Power(dy,2)*Power(px1,2)*px2*py2+162*by*cx*dy*Power(px2,2)*py2+81*ax*Power(dy,2)*Power(px2,2)*py2+108*by*cx*dy*px1*px3*py2+54*ax*Power(dy,2)*px1*px3*py2-162*Power(dy,2)*Power(px1,2)*px3*py2-108*by*cx*dy*px2*px3*py2-54*ax*Power(dy,2)*px2*px3*py2+81*Power(dy,2)*px1*px2*px3*py2-27*Power(dy,2)*Power(px2,2)*px3*py2+18*by*cx*dy*Power(px3,2)*py2+9*ax*Power(dy,2)*Power(px3,2)*py2+18*Power(dy,2)*px1*Power(px3,2)*py2-324*by*cx*dx*px1*py1*py2-324*bx*cx*dy*px1*py1*py2-324*ax*dx*dy*px1*py1*py2+324*by*cx*dx*px2*py1*py2+324*bx*cx*dy*px2*py1*py2+324*ax*dx*dy*px2*py1*py2-108*by*cx*dx*px3*py1*py2-108*bx*cx*dy*px3*py1*py2-108*ax*dx*dy*px3*py1*py2+243*dx*dy*px1*px3*py1*py2-81*dx*dy*px2*px3*py1*py2-27*dx*dy*Power(px3,2)*py1*py2+81*Power(cx,3)*Power(py1,2)*py2+486*bx*cx*dx*Power(py1,2)*py2+243*ax*Power(dx,2)*Power(py1,2)*py2-81*Power(dx,2)*px2*Power(py1,2)*py2-81*Power(dx,2)*px3*Power(py1,2)*py2+162*by*cx*dx*px1*Power(py2,2)+162*bx*cx*dy*px1*Power(py2,2)+162*ax*dx*dy*px1*Power(py2,2)-81*dx*dy*Power(px1,2)*Power(py2,2)-162*by*cx*dx*px2*Power(py2,2)-162*bx*cx*dy*px2*Power(py2,2)-162*ax*dx*dy*px2*Power(py2,2)+54*by*cx*dx*px3*Power(py2,2)+54*bx*cx*dy*px3*Power(py2,2)+54*ax*dx*dy*px3*Power(py2,2)-81*dx*dy*px1*px3*Power(py2,2)+54*dx*dy*px2*px3*Power(py2,2)-81*Power(cx,3)*py1*Power(py2,2)-486*bx*cx*dx*py1*Power(py2,2)-243*ax*Power(dx,2)*py1*Power(py2,2)+81*Power(dx,2)*px1*py1*Power(py2,2)+81*Power(dx,2)*px3*py1*Power(py2,2)+27*Power(cx,3)*Power(py2,3)+162*bx*cx*dx*Power(py2,3)+81*ax*Power(dx,2)*Power(py2,3)-27*Power(dx,2)*px3*Power(py2,3)-54*by*cx*dy*Power(px1,2)*py3-27*ax*Power(dy,2)*Power(px1,2)*py3-54*Power(dy,2)*Power(px1,3)*py3+108*by*cx*dy*px1*px2*py3+54*ax*Power(dy,2)*px1*px2*py3+81*Power(dy,2)*Power(px1,2)*px2*py3-54*by*cx*dy*Power(px2,2)*py3-27*ax*Power(dy,2)*Power(px2,2)*py3-81*Power(dy,2)*px1*Power(px2,2)*py3+27*Power(dy,2)*Power(px2,3)*py3-36*by*cx*dy*px1*px3*py3-18*ax*Power(dy,2)*px1*px3*py3+54*Power(dy,2)*Power(px1,2)*px3*py3+36*by*cx*dy*px2*px3*py3+18*ax*Power(dy,2)*px2*px3*py3-27*Power(dy,2)*px1*px2*px3*py3-6*by*cx*dy*Power(px3,2)*py3-3*ax*Power(dy,2)*Power(px3,2)*py3+108*by*cx*dx*px1*py1*py3+108*bx*cx*dy*px1*py1*py3+108*ax*dx*dy*px1*py1*py3+108*dx*dy*Power(px1,2)*py1*py3-108*by*cx*dx*px2*py1*py3-108*bx*cx*dy*px2*py1*py3-108*ax*dx*dy*px2*py1*py3-243*dx*dy*px1*px2*py1*py3+81*dx*dy*Power(px2,2)*py1*py3+36*by*cx*dx*px3*py1*py3+36*bx*cx*dy*px3*py1*py3+36*ax*dx*dy*px3*py1*py3+9*dx*dy*px2*px3*py1*py3-27*Power(cx,3)*Power(py1,2)*py3-162*bx*cx*dx*Power(py1,2)*py3-81*ax*Power(dx,2)*Power(py1,2)*py3-54*Power(dx,2)*px1*Power(py1,2)*py3+162*Power(dx,2)*px2*Power(py1,2)*py3-54*Power(dx,2)*px3*Power(py1,2)*py3-108*by*cx*dx*px1*py2*py3-108*bx*cx*dy*px1*py2*py3-108*ax*dx*dy*px1*py2*py3+81*dx*dy*Power(px1,2)*py2*py3+108*by*cx*dx*px2*py2*py3+108*bx*cx*dy*px2*py2*py3+108*ax*dx*dy*px2*py2*py3+81*dx*dy*px1*px2*py2*py3-54*dx*dy*Power(px2,2)*py2*py3-36*by*cx*dx*px3*py2*py3-36*bx*cx*dy*px3*py2*py3-36*ax*dx*dy*px3*py2*py3-9*dx*dy*px1*px3*py2*py3+54*Power(cx,3)*py1*py2*py3+324*bx*cx*dx*py1*py2*py3+162*ax*Power(dx,2)*py1*py2*py3-81*Power(dx,2)*px1*py1*py2*py3-81*Power(dx,2)*px2*py1*py2*py3+27*Power(dx,2)*px3*py1*py2*py3-27*Power(cx,3)*Power(py2,2)*py3-162*bx*cx*dx*Power(py2,2)*py3-81*ax*Power(dx,2)*Power(py2,2)*py3+27*Power(dx,2)*px2*Power(py2,2)*py3+18*by*cx*dx*px1*Power(py3,2)+18*bx*cx*dy*px1*Power(py3,2)+18*ax*dx*dy*px1*Power(py3,2)-54*dx*dy*Power(px1,2)*Power(py3,2)-18*by*cx*dx*px2*Power(py3,2)-18*bx*cx*dy*px2*Power(py3,2)-18*ax*dx*dy*px2*Power(py3,2)+27*dx*dy*px1*px2*Power(py3,2)+6*by*cx*dx*px3*Power(py3,2)+6*bx*cx*dy*px3*Power(py3,2)+6*ax*dx*dy*px3*Power(py3,2)-9*Power(cx,3)*py1*Power(py3,2)-54*bx*cx*dx*py1*Power(py3,2)-27*ax*Power(dx,2)*py1*Power(py3,2)+54*Power(dx,2)*px1*py1*Power(py3,2)-18*Power(dx,2)*px2*py1*Power(py3,2)+9*Power(cx,3)*py2*Power(py3,2)+54*bx*cx*dx*py2*Power(py3,2)+27*ax*Power(dx,2)*py2*Power(py3,2)-9*Power(dx,2)*px1*py2*Power(py3,2)-Power(cx,3)*Power(py3,3)-6*bx*cx*dx*Power(py3,3)-3*ax*Power(dx,2)*Power(py3,3)-3*cx*Power(cy,2)*Power(tmp1,2)*tmp2+3*ay*(tmp1)*Power(dy*(tmp1)-dx*tmp2,2)+3*cy*(tmp1)*(2*by*(tmp1)*(dy*(tmp1)-dx*tmp2)+tmp2*(-2*bx*dy*(tmp1)+Power(cx,2)*tmp2+2*bx*dx*tmp2)));T c7=+3*(dy*(tmp1)-dx*tmp2)*(Power(cy,2)*Power(tmp1,2)-2*cx*cy*(tmp1)*tmp2+by*(tmp1)*(dy*(tmp1)-dx*tmp2)+tmp2*(-(bx*dy*(tmp1))+Power(cx,2)*tmp2+bx*dx*tmp2));T c8=+3*(cy*(tmp1)-cx*tmp2)*Power(dy*(tmp1)-dx*tmp2,2);T tmp=dy*(tmp1)-dx*tmp2;T c9=tmp*tmp*tmp;return{c0,c1,c2,c3,c4,c5,c6,c7,c8,c9};

            // clang-format on
        }

        template<typename T>
        void generate_coefficients(const CubicBezier &bezier, const CubicBezier &q, std::array<T, 10> &out) {
            T dx = q.p0.x();
            T dy = q.p0.y();
            T b0x = bezier.p0.x(), b1x = bezier.p1.x(), b2x = bezier.p2.x(), b3x = bezier.p3.x();
            T b0y = bezier.p0.y(), b1y = bezier.p1.y(), b2y = bezier.p2.y(), b3y = bezier.p3.y();
            T px1 = q.p1.x();
            T py1 = q.p1.y();
            T px2 = q.p2.x();
            T py2 = q.p2.y();
            T px3 = q.p3.x();
            T py3 = q.p3.y();

            b0x -= dx;
            b0y -= dy;
            b1x -= dx;
            b1y -= dy;
            b2x -= dx;
            b2y -= dy;
            b3x -= dx;
            b3y -= dy;

            px1 -= dx;
            py1 -= dy;
            px2 -= dx;
            py2 -= dy;
            px3 -= dx;
            py3 -= dy;

            out = generate_coefficients<T>(
                    b0x, b0y, b1x, b1y, b2x, b2y, b3x, b3y,
                    px1, py1, px2, py2, px3, py3
            );
            T absmax = out.front();
            using std::abs;
            for (std::size_t i = 1; i < out.size(); ++i) {
                if (abs(out[i]) > absmax) {
                    absmax = abs(out[i]);
                }
            }
            if (absmax > 1) {
                for (std::size_t i = 0; i < out.size(); ++i) {
                    out[i] /= absmax;
                }
            }
        }

        template<typename IntersectionConsumer, class numeric_t = double, typename CollinearFunc = IsCollinear>
        void get_intersection_detail(const CubicBezier &bezier, const CubicBezier &q,
                                     IntersectionConsumer &report,
                                     const numeric_t &abstol, const double interval_eps = 1e-14,
                                     const int maxiter = 25) {
            // q is the curve that is transformed in implicit form, bezier is the curve on which t values are found
            std::array<numeric_t, 10> poly;
            generate_coefficients(bezier, q, poly);

            std::array<double, 10> poly_d{double(poly[0]), double(poly[1]), double(poly[2]), double(poly[3]),
                                          double(poly[4]),
                                          double(poly[5]), double(poly[6]), double(poly[7]), double(poly[8]),
                                          double(poly[9])};
            polynomialsolver::PolynomialFunc<10, double> polyfunc_d(poly_d);
            polynomialsolver::PolynomialFunc<10, numeric_t> polyfunc(poly);

            auto inv = curveinverter<CollinearFunc>(q);

            auto handle_root = [&](const double t) {
                Point2d p = beziermap(bezier, t);
                double u = inv(p.x(), p.y());
                if (unit_inter(u)) {
                    report(t, u, p);
                }
            };

            double prev_root = -1.;
            using std::abs;
            auto add = [&](const std::pair<numeric_t, numeric_t> &interval) {
#ifdef DEBUG
                std::cout << "interval: " << double(interval.first) << ", " << double(interval.second) << '\n';
#endif
                if (interval.first == numeric_t(0)) {
#ifdef DEBUG
                    std::cout << "boundary val: " << polyfunc(numeric_t(0)) << '\n';
#endif
                    if (abs(polyfunc(numeric_t(0))) <= abstol) {
                        if (prev_root != 0.) {
                            prev_root = 0.;
                            if (!(bezier.p0 == q.p0 || bezier.p0 == q.p3)) {
                                handle_root(0.);
                            }
                        }
                        return;
                    }
                }

                if (interval.second == numeric_t(1)) {
#ifdef DEBUG
                    std::cout << "boundary val: " << polyfunc(numeric_t(1)) << '\n';
#endif
                    if (abs(polyfunc(numeric_t(1))) <= abstol) {
                        if (prev_root != 1.) {
                            prev_root = 1.;
                            if (!(bezier.p3 == q.p3 || bezier.p3 == q.p0)) {
                                handle_root(1.);
                            }
                        }
                        return;
                    }
                }

                auto maybe_root = polynomialsolver::itp_root_refine(polyfunc_d,
                                                                    double(interval.first), double(interval.second),
                                                                    interval_eps, maxiter);
                if (maybe_root) {
                    auto t = (double) *maybe_root;
#ifdef DEBUG
                    std::cout << "computed general root " << t << '\n';
#endif
                    if (prev_root != t) {
                        prev_root = t;
                        handle_root(t);
                    }
                }
            };
            polynomialsolver::rootbracket(poly, add, numeric_t(abstol));
        }

        // Intersects the curves b1 and b2. Intersections are passed to the IntersectionConsumer callback,
        // more precisely the parameters t, u, p are passed, where t is the parametric value of the first curve,
        // u the parametric value of the second curve, and p is the intersection point_.
        // preconditions: - curve must not be equivalent to a line -curve cannot intersect itself in another point_ than
        // the start-end points.
        template<typename IntersectionConsumer, class numeric_t = double, typename CollinearFunc = IsCollinear>
        void curve_curve_inter(const CubicBezier &b1, const CubicBezier &b2, IntersectionConsumer &report,
                               const numeric_t &abstol = (numeric_t) 1e-10) {
            if (b1 == b2) {
                return;
            }
            BBox box1{}, box2{};
            make_hull_bbox(b1, box1);
            make_hull_bbox(b2,box2);
            bool nointersect = !box1.strict_overlap(box2);
            if (nointersect) {
                return;
            }

            bool switch_order = !bezier_ordering(b1, b2);
            auto report_impl = [&](const double t, const double u, const Point2d &p) {
                if (switch_order) {
                    return report(u, t, p);
                } else {
                    return report(t, u, p);
                }
            };

            CubicBezier b1_ = b1, b2_ = b2;
            if (switch_order) std::swap(b1_, b2_);
            get_intersection_detail<decltype(report_impl), double, CollinearFunc>(b1_, b2_, report_impl, abstol);
        }
    }
#endif //CONTOURKLIP_BEZIER_UTILS_HPP