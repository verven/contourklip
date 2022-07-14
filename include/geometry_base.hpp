#ifndef CONTOURKLIP_GEOMETRY_BASE_HPP
#define CONTOURKLIP_GEOMETRY_BASE_HPP

#include <vector>
#include <array>
#include <cassert>
#include <optional>
#include <iostream>
#include <algorithm>
#include "direct_solvers.hpp"

namespace contourklip {
    struct Point2d {
    private:
        double _x;
        double _y;
    public:
        constexpr Point2d(double x, double y) : _x(x), _y(y) {}

        constexpr Point2d() : _x(0), _y(0) {}

        inline constexpr double x() const {
            return _x;
        }

        inline constexpr double y() const {
            return _y;
        }

        inline constexpr double& x() {
            return _x;
        }

        inline constexpr double& y() {
            return _y;
        }

        template<int idx>
        friend constexpr double get(const Point2d &p) {
            static_assert(idx == 0 || idx == 1);
            if constexpr(idx == 0) {
                return p.x();
            } else {
                return p.y();
            }
        }
    };

    template<int idx>
    constexpr double get(const Point2d &p);

    std::ostream &operator<<(std::ostream &o, const Point2d &p) {
        return o << "(" << p.x() << ", " << p.y() << ")";
    }

    inline bool operator==(const Point2d &p1, const Point2d &p2) {
        return (p1.x() == p2.x()) && (p1.y() == p2.y());
    }

    inline bool operator!=(const Point2d &p1, const Point2d &p2) { return !(p1 == p2); }

    constexpr auto increasing = [](const Point2d &a, const Point2d &b) {
        if (a.x() == b.x()) {
            return a.y() < b.y();
        }
        return a.x() < b.x();
    };

    bool operator<(const Point2d &a, const Point2d &b) { return increasing(a, b); }

    namespace detail {

        inline bool approx_equal(const Point2d &p1, const Point2d &p2, double eps) {
            return std::abs(p1.x() - p2.x()) < eps && std::abs(p1.y() - p2.y()) < eps;
        }

        inline double signed_area(const Point2d &p0, const Point2d &p1, const Point2d &p2) {
            using namespace directsolvers;
            return diff_of_products(p0.x() - p2.x(), p1.y() - p2.y(), p1.x() - p2.x(), p0.y() - p2.y());
        }
//
//        auto left_of_line = [](const Point2d &p0, const Point2d &p1, const Point2d &a) -> bool {
//            return signed_area(p0, p1, a) > 0;
//        };
//
//        auto is_collinear = [](const Point2d &a, const Point2d &b, const Point2d &p) -> bool {
//            if (a.x() == b.x()) {
//                return a.x() == p.x();
//            }
//            if (a.y() == b.y()) {
//                return a.y() == p.y();
//            }
//            return signed_area(a, b, p) == 0;
//        };

        struct LeftOfLine{
            bool operator()(const Point2d &p0, const Point2d &p1, const Point2d &a) const{
                return signed_area(p0, p1, a) > 0;
            }
        };

        struct IsCollinear{
            bool operator()(const Point2d &a, const Point2d &b, const Point2d &p) const{
                if (a.x() == b.x()) {
                    return a.x() == p.x();
                }
                if (a.y() == b.y()) {
                    return a.y() == p.y();
                }
                return signed_area(a, b, p) == 0;
            }
        };

        template<typename Orient2dFunc = LeftOfLine>
        inline bool above_line(const Point2d &a, const Point2d &b, const Point2d &p, const Orient2dFunc &on_left = {}) {
            return increasing(a, b) ? on_left(a, b, p) : on_left(b, a, p);
        }

        double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3) {
            return 0.5 * (x1 * (y2 - y3) +
                          x2 * (y3 - y1) +
                          x3 * (y1 - y2));
        }

        double quadri_area(const Point2d &a, const Point2d &b, const Point2d &c, const Point2d &d) {
            return triangle_area(a.x(), a.y(), b.x(), b.y(), c.x(), c.y())
                   + triangle_area(c.x(), c.y(), d.x(), d.y(), a.x(), a.y());
        }

        inline double sqdist(const Point2d &a, const Point2d &b) {
            double x = a.x() - b.x();
            double y = a.y() - b.y();
            return x * x + y * y;
        }

        inline double dist(const Point2d &a, const Point2d &b) {
            return sqrt(sqdist(a, b));
        }

        double sqdist_to(const Point2d &a, const Point2d &b, const Point2d &p) {
            double dx = b.x() - a.x();
            double dy = b.y() - a.y();
            double num = dy * p.x() - dx * p.y() + b.x() * a.y() - b.y() * a.x();
            double den = dy * dy + dx * dx;
            return num * num / den;
        }

        inline bool in_range(double x, double a, double b) {
            return a < x && x < b;
        }

        inline bool in_range_strict(double a, double b, double x) {
            return a < x && x < b;
        }

        inline bool in_range_closed(double a, double b, double x) {
            return a <= x && x <= b;
        }

        inline bool in_interval(double a, double b, double x) {
            return a < b ? in_range(x, a, b) : in_range(x, b, a);
        }

        inline bool in_box(const Point2d &a, const Point2d &b, const Point2d &p) {
            return in_interval(a.x(), b.x(), p.x()) && in_interval(a.y(), b.y(), p.y());
        }

        Point2d basic_intersection(const Point2d &p1, const Point2d &p2, const Point2d &p3, const Point2d &p4) {
            double den = ((p1.x() - p2.x()) * (p3.y() - p4.y()) - (p1.y() - p2.y()) * (p3.x() - p4.x()));
            double px = ((p1.x() * p2.y() - p1.y() * p2.x()) * (p3.x() - p4.x())
                         - (p1.x() - p2.x()) * (p3.x() * p4.y() - p3.y() * p4.x()))
                        / den;
            double py = ((p1.x() * p2.y() - p1.y() * p2.x()) * (p3.y() - p4.y()) -
                         (p1.y() - p2.y()) * (p3.x() * p4.y() - p3.y() * p4.x()))
                        / den;
            return {px, py};
        }

        class BBox {
        public:
            double min_x;
            double min_y;
            double max_x;
            double max_y;

            inline bool weak_contains(const Point2d &p) const {
                return min_x <= p.x() && p.x() <= max_x
                       && min_y <= p.y() && p.y() <= max_y;
            }

            inline bool strict_contains(const Point2d &p) const {
                return min_x < p.x() && p.x() < max_x
                       && min_y < p.y() && p.y() < max_y;
            }

            inline bool strict_contains_x(const double &x) const {
                return min_x < x && x < max_x;
            }

            inline bool strict_contains_y(const double &y) const {
                return min_y < y && y < max_y;
            }

            inline bool weak_contains_x(const double &x) const {
                return min_x <= x && x <= max_x;
            }

            inline bool weak_contains_y(const double &y) const {
                return min_y <= y && y <= max_y;
            }

            inline bool weak_overlap(const BBox &other) const {
                bool vertical = weak_contains_x(other.min_x)
                                || weak_contains_x(other.max_x)
                                || other.weak_contains_x(this->min_x)
                                || other.weak_contains_x(this->max_x);
                bool horizontal = weak_contains_y(other.min_y)
                                  || weak_contains_y(other.max_y)
                                  || other.weak_contains_y(this->min_y)
                                  || other.weak_contains_y(this->max_y);
                return vertical && horizontal;
            }

            inline bool strict_overlap(const BBox &other) const {
                if (weak_overlap(other)) {
                    bool vertical = strict_contains_x(other.min_x)
                                    || strict_contains_x(other.max_x)
                                    || other.strict_contains_x(this->min_x)
                                    || other.strict_contains_x(this->max_x);
                    bool horizontal = strict_contains_y(other.min_y)
                                      || strict_contains_y(other.max_y)
                                      || other.strict_contains_y(this->min_y)
                                      || other.strict_contains_y(this->max_y);
                    return vertical || horizontal;
                }
                return false;
            }


            friend std::ostream &operator<<(std::ostream &o, const BBox &bbox) {
                return o << "[" << Point2d{bbox.min_x, bbox.min_y}
                         << ", " << Point2d{bbox.max_x, bbox.max_y} << "]";
            }
        };

        struct Segment {
            Point2d first;
            Point2d second;

            friend std::ostream &operator<<(std::ostream &o, const Segment &p) {
                return o << "[" << p.first << ", " << p.second << "]";
            }
        };

        inline Point2d linear_map(const Point2d &first, const Point2d &second, double t) {
            return {first.x() + t * (second.x() - first.x()), first.y() + t * (second.y() - first.y())};
        }

        inline Point2d linear_map(const Segment &seg, double t) {
            return linear_map(seg.first, seg.second, t);
        }

        inline bool vertical(const Point2d &a, const Point2d &b) {
            return a.x() == b.x();
        }

        inline bool horizontal(const Point2d &a, const Point2d &b) {
            return a.y() == b.y();
        }

        template<typename T>
        T segment_tval(const T &a_x, const T &a_y, const T &b_x, const T &b_y, const T &p_x, const T &p_y) {
            return a_x;
        }

        double segment_tval(const Segment &seg, const Point2d &p) {
            if (p == seg.first) return 0;
            if (p == seg.second) return 1;
            if (seg.second.x() == seg.first.x()) {
                return (p.y() - seg.first.y()) / (seg.second.y() - seg.first.y());
            }
            return (p.x() - seg.first.x()) / (seg.second.x() - seg.first.x());
        }

        BBox make_bbox(const Segment &seg) {
            double min_x = std::min(seg.first.x(), seg.second.x());
            double min_y = std::min(seg.first.y(), seg.second.y());
            double max_x = std::max(seg.first.x(), seg.second.x());
            double max_y = std::max(seg.first.y(), seg.second.y());
            return {min_x, min_y, max_x, max_y};
        }

        struct SegInter {
            double t1;
            double t2;
            Point2d p;
        };

        std::ostream &operator<<(std::ostream &o, const SegInter &q) {
            return o << "[" << q.p << " " << q.t1 << " " << q.t2 << "]";
        }

        bool operator==(const SegInter &a, const SegInter &b) {
            return a.t1 == b.t1 && a.t2 == b.t2 && a.p == b.p;
        }

        // returns a SegInter if the following 3 conditions hold:
        // a) the segments do not share any endpoint
        // b) the segments are not parallel
        // c) at least one of the segments is split in 2 new segments by the other segment
        // note that points are passed by value
        template<typename collinearF = IsCollinear>
        std::optional<SegInter> intersect_segments_detail(Point2d a1,
                                                          Point2d a2,
                                                          Point2d b1,
                                                          Point2d b2) {
            if (a1 == b1
                || a2 == b2
                || a1 == b2
                || a2 == b1
                    ) {
                return {};
            }
            collinearF collinear;
            bool b1_on_a = collinear(a1, a2, b1);
            bool b2_on_a = collinear(a1, a2, b2);
            if (b1_on_a && b2_on_a) {
                return {};
            }

            auto make_inter = [&](double t, double u, const Point2d &p) -> SegInter {
                double a = t, b = u;
                return SegInter{a, b, p};
            };

            double x1 = a1.x(), y1 = a1.y();
            double x2 = a2.x(), y2 = a2.y();
            double x3 = b1.x(), y3 = b1.y();
            double x4 = b2.x(), y4 = b2.y();

            using namespace directsolvers;
            double den = diff_of_products((x1 - x2), (y3 - y4), (y1 - y2), (x3 - x4));
            double num1 = diff_of_products((x1 - x3), (y3 - y4), (y1 - y3), (x3 - x4));
            double num2 = diff_of_products((x2 - x1), (y1 - y3), (y2 - y1), (x1 - x3));

            bool a_notinrange = num1 * den < 0 || std::abs(num1) > std::abs(den);
            bool b_notinrange = num2 * den < 0 || std::abs(num2) > std::abs(den);
            if (b1_on_a) {
                if (a_notinrange) {
                    return {};
                }
                return make_inter(num1 / den, 0, b1);
            }
            if (b2_on_a) {
                if (a_notinrange) {
                    return {};
                }
                return make_inter(num1 / den, 1, b2);
            }
            bool a1_on_b = collinear(b1, b2, a1);
            bool a2_on_b = collinear(b1, b2, a2);
            if (a1_on_b) {
                if (b_notinrange) {
                    return {};
                }
                return make_inter(0, num2 / den, a1);
            }
            if (a2_on_b) {
                if (b_notinrange) {
                    return {};
                }
                return make_inter(1, num2 / den, a2);
            }
            if (a_notinrange || b_notinrange) {
                return {};
            }
            return make_inter(num1 / den, num2 / den, linear_map(a1, a2, num1 / den));
        }

        template<typename collinearF = IsCollinear>
        std::optional<SegInter> intersect_segments(Point2d a1, Point2d a2, Point2d b1, Point2d b2) {
            if (a1 == b1
                || a2 == b2
                || a1 == b2
                || a2 == b1
                    ) {
                return {};
            }
            collinearF collinear;
            bool b1_on_a = collinear(a1, a2, b1);
            bool b2_on_a = collinear(a1, a2, b2);
            if (b1_on_a && b2_on_a) {
                return {};
            }
            bool a_reversed = false, b_reversed = false;
            bool segments_swapped = false;

            if ((a_reversed = !increasing(a1, a2))) {
                std::swap(a1, a2);
            }
            if ((b_reversed = !increasing(b1, b2))) {
                std::swap(b1, b2);
            }
            double dx1 = a2.x() - a1.x(), dy1 = a2.y() - a1.y();
            double dx2 = b2.x() - b1.x(), dy2 = b2.y() - b1.y();

            if (!increasing({dx1, dy1}, {dx2, dy2})) {
                std::swap(a1, b1);
                std::swap(a2, b2);
                segments_swapped = true;
            }
            if (auto ret = intersect_segments_detail(a1, a2, b1, b2)) {
                if (segments_swapped) {
                    std::swap(ret->t1, ret->t2);
                }
                if (a_reversed) {
                    ret->t1 = 1 - ret->t1;
                }
                if (b_reversed) {
                    ret->t2 = 1 - ret->t2;
                }
                return ret;
            }
            return {};
        }

        template<typename collinearF = IsCollinear>
        std::optional<SegInter> intersect_segments(const Segment &a, const Segment &b) {
            return intersect_segments(a.first, a.second, b.first, b.second);
        }

        struct CubicBezier {
            Point2d p0;
            Point2d p1;
            Point2d p2;
            Point2d p3;

            CubicBezier() = default;

            CubicBezier(const Point2d &p0, const Point2d &p1, const Point2d &p2, const Point2d &p3) : p0(p0), p1(p1),
                                                                                                      p2(p2), p3(p3) {}

            explicit CubicBezier(std::array<Point2d, 4> &in) {
                p0 = in[0];
                p1 = in[1];
                p2 = in[2];
                p3 = in[3];
            }

            constexpr std::array<Point2d, 4> as_array() const {
                return {p0, p1, p2, p3};
            }

            friend std::ostream &operator<<(std::ostream &o, const CubicBezier &p) {
                return o << "[" << p.p0 << " " << p.p1 << " " << p.p2 << " " << p.p3 << "]";
            }

            friend bool operator==(const CubicBezier &a, const CubicBezier &b) {
                return a.as_array() == b.as_array();
            }
        };

        inline void make_hull_bbox(const CubicBezier &c, BBox &out) {
            double min_x = std::min(c.p0.x(), c.p3.x());
            double min_y = std::min(c.p0.y(), c.p3.y());
            double max_x = std::max(c.p0.x(), c.p3.x());
            double max_y = std::max(c.p0.y(), c.p3.y());
            out.min_x = std::min(min_x, std::min(c.p1.x(), c.p2.x()));
            out.min_y = std::min(min_y, std::min(c.p1.y(), c.p2.y()));
            out.max_x = std::max(max_x, std::max(c.p1.x(), c.p2.x()));
            out.max_y = std::max(max_y, std::max(c.p1.y(), c.p2.y()));
        }
    }
    class Contour;

    enum ComponentType {
        LINE = 0,
        CUBIC_BEZIER = 1,
    };

    /// \brief a simple struct to represent a segment of a path.
    struct ContourComponent {
        friend class Contour;
    private:
        ComponentType component_type_;
        Point2d c_1_;
        Point2d c_2_;
        Point2d point_;
    public:
        /// \brief constructs an instance this representing a line segment of a Contour.
        /// If given a Contour c, adding this to it will represent the segment [c.back_point(), this->point()]
        /// \param pLast the point_ Point2d representing the end point
        explicit ContourComponent(const Point2d &p) : component_type_(LINE), point_(p) {}

        /// \brief constructs an instance representing a cubic bezier segment of a contour.
        /// given some first point p, it will represent the bezier [p, c1, c2, point].
        /// \param p_1 the first Point2d control point
        /// \param p_2 the second Point2d control point
        /// \param p_last the Point2d endpoint
        ContourComponent(const Point2d &c_1,
                         const Point2d &c_2,
                         const Point2d &p) :
                component_type_(CUBIC_BEZIER), c_1_(c_1), c_2_(c_2), point_(p) {}

        Point2d c1() const {
#ifdef DEBUG
            assert(segment_shape() == CUBIC_BEZIER);
#endif
            return c_1_;
        }

        Point2d& c1() {
            return c_1_;
        }

        Point2d c2() const {
#ifdef DEBUG
            assert(segment_shape() == CUBIC_BEZIER);
#endif
            return c_2_;
        }

        Point2d& c2() {
            return c_2_;
        }

        Point2d point() const {
            return point_;
        }

        Point2d &point() {
            return point_;
        }

        bool bcurve() const {
            return segment_shape() == CUBIC_BEZIER;
        }

        /// \brief returns the shape type tag associated with this instance.
        /// \return the shape type ComponentType
        ComponentType segment_shape() const {
            return component_type_;
        }

        friend bool operator==(const ContourComponent &a, const ContourComponent &b) {
            if (a.component_type_ != b.component_type_) {
                return false;
            }
            if (a.bcurve()) {
                return a.c1() == b.c1()
                       && a.c2() == b.c2()
                       && a.point() == b.point();
            }
            return a.point() == b.point();
        }

    private:
        /// \brief reverses the control points associated with this, irrespective of the shape type.
        void reverse_controlp() {
            Point2d temp = c_1_;
            c_1_ = c_2_;
            c_2_ = temp;
        }

    };

    std::ostream &operator<<(std::ostream &o, const ContourComponent &comp) {
        switch (comp.segment_shape()) {
            case CUBIC_BEZIER:
                return o << comp.c1() << " " << comp.c2() << " " << comp.point();
            case LINE:
                return o << comp.point();
        }
        return o;
    }

    class Contour {
    private:
        using ContainerType = std::vector<ContourComponent>;
        ContainerType container{};
    public:
        Contour() = default;

        explicit Contour(const Point2d &start) {
            push_back(start);
        }

        Contour(const Point2d &p0, const Point2d &p1, const Point2d &p2, const Point2d &p3) {
            push_back(p0);
            push_back(p1, p2, p3);
        }

        Contour(const Point2d &p0, const Point2d &p1) {
            push_back(p0);
            push_back(p1);
        }

        void push_back(const Point2d &p) {
            push_back(ContourComponent(p));
        }

        void push_back(const Point2d &p2,
                       const Point2d &p3,
                       const Point2d &p) {
            push_back(
                    ContourComponent{p2, p3, p}
            );
        }

        void push_back(const ContourComponent &start) {
            container.push_back(start);
        }

        ContourComponent operator[](const std::size_t idx) const {
            return container[idx];
        }

        ContourComponent& operator[](const std::size_t idx) {
            return container[idx];
        }

        std::size_t size() const {
            return container.size();
        }

        Point2d front_point() const {
            return container.front().point();
        }

        Point2d back_point() const {
            return container.back().point();
        }


        ContourComponent &front() {
            return container.front();
        }

        ContourComponent &back() {
            return container.back();
        }

        ContourComponent front() const {
            return container.front();
        }

        ContourComponent back() const {
            return container.back();
        }

        bool is_closed() const {
            return front_point() == back_point();
        }

        void close() {
            if (!is_closed()) {
                this->push_back(front_point());
            }
        }

        void reverse() {
            std::reverse(container.begin(), container.end());
            for (std::size_t i = 0; i < container.size() - 1; ++i) {
                container[i].point() = container[i + 1].point();
                if (container[i].segment_shape() == CUBIC_BEZIER) {
                    container[i].reverse_controlp();
                }
            }
            //rotate to the right so that we start with a simple point.
            std::rotate(container.rbegin(), container.rbegin() + 1, container.rend());
        }

        template<ComponentType T, typename Consumer>
        void forward_segments(Consumer &out) const {
            if (container.empty()) { return; }
            auto it = container.begin();
            for (auto prev = it++; it != container.end(); prev++, it++) {
                auto pair = std::make_pair(*prev, *it);
                if (pair.second.segment_shape() == T) {
                    if constexpr (T == LINE) {
                        out(pair.first.point(), pair.second.point());
                    } else {
                        out(pair.first.point(),
                            pair.second.c1(),
                            pair.second.c2(),
                            pair.second.point());
                    }
                }
            }
        }

        auto begin() const {
            return container.begin();
        }

        auto begin() {
            return container.begin();
        }

        auto end() const {
            return container.end();
        }

        auto end() {
            return container.end();
        }

        friend bool operator==(const Contour &a, const Contour &b) {
            return a.container == b.container;
        }

        friend std::ostream &operator<<(std::ostream &o, const Contour &c) {
            for (const auto &seg: c) {
                std::cout << seg << '\n';
            }
            return o;
        }
    };


    namespace detail {
        std::tuple<double, double, double, double> contourbbox(const Contour &a) {
            double min_x = a.front_point().x();
            double max_x = min_x;
            double min_y = a.front_point().y();
            double max_y = min_y;
            for (const auto &seg: a) {
                min_x = std::min(min_x, seg.point().x());
                min_y = std::min(min_y, seg.point().y());
                max_x = std::max(max_x, seg.point().x());
                max_y = std::max(max_y, seg.point().y());
                switch (seg.segment_shape()) {
                    case LINE:
                        continue;
                    case CUBIC_BEZIER:
                        min_x = std::min(min_x, std::min(seg.c1().x(), seg.c2().x()));
                        min_y = std::min(min_y, std::min(seg.c1().y(), seg.c2().y()));
                        max_x = std::max(max_x, std::max(seg.c1().x(), seg.c2().x()));
                        max_y = std::max(max_y, std::max(seg.c1().y(), seg.c2().y()));
                        continue;
                }
            }
            return {min_x, min_y, max_x, max_y};
        }

        double bezier_area(const Point2d &p0, const Point2d &p1,
                           const Point2d &p2, const Point2d p3) {
            double x0 = p0.x(), y0 = p0.y(), x1 = p1.x(), y1 = p1.y(),
                    x2 = p2.x(), y2 = p2.y(), x3 = p3.x(), y3 = p3.y();
            return (x0 * (-2 * y1 - y2 + 3 * y3)
                    + x1 * (2 * y0 - y2 - y3)
                    + x2 * (y0 + y1 - 2 * y3)
                    + x3 * (-3 * y0 + y1 + 2 * y2)
                   ) * 3. / 20.;
        }

        double contour_area(const Contour &c) {
            if (c.size() < 2) return 0.;
            double area = 0.0;
            for (std::size_t i = 0; i < c.size() - 1; ++i) {
                area += c[i].point().x() * c[i + 1].point().y()
                        - c[i + 1].point().x() * c[i].point().y();
                if (c[i + 1].segment_shape() == CUBIC_BEZIER) {
                    double t = bezier_area(c[i].point(), c[i + 1].c1(),
                                           c[i + 1].c2(), c[i + 1].point());
                    //mult. by 2 since at the end we div by 2.
                    area -= 2 * t;
                }
            }
            if (!c.is_closed()) {
                // we only have an implicit line segment
                area += c.back_point().x() * c.front_point().y()
                        - c.back_point().x() * c.front_point().y();
            }
            return 0.5 * area;
        }

        double multipolygon_area(const std::vector<Contour> &poly) {
            double area = 0;
            for (const auto &c: poly) {
                area += contour_area(c);
            }
            return area;
        }
    }
}
#endif //CONTOURKLIP_GEOMETRY_BASE_HPP