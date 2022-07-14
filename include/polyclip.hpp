#ifndef CONTOURKLIP_POLYCLIP_HPP
#define CONTOURKLIP_POLYCLIP_HPP

#include <map>
#include <deque>
#include <algorithm>
#include "bezier_utils.hpp"
#include "geometry_base.hpp"
#include "sweep_point.hpp"
#include "contour_postprocessing.hpp"
#ifdef DEBUG
#include "svg_io.hpp"
#endif

namespace contourklip {
    using namespace detail;

#ifdef DEBUG
#define CONTOURKLIP_IF_NOT(stmt, on_fail) assert((stmt) && (on_fail)); if(false)
#else
#define CONTOURKLIP_IF_NOT(stmt, msg) if( ! (stmt) )
#endif

    template<typename T = Segment>
    struct ContourSegment{
        T seg;
        std::size_t contourid;
        detail::PolygonType ptype;
    };

    struct SegIndent {
        detail::PolygonType ptype;
        bool curve;
        std::size_t idx;
        std::size_t contour_id;
    };

    bool operator==(const SegIndent &a, const SegIndent &b) {
        return a.idx == b.idx && a.curve == b.curve && a.ptype == b.ptype;
    }

    struct IntersectionValue {
        double t_val;
        Point2d p{};
        SegIndent id;
    };

    bool operator<(const IntersectionValue &a, const IntersectionValue &b) {
        if (a.id.ptype != b.id.ptype) {
            return a.id.ptype < b.id.ptype;
        }
        if (a.id.curve != b.id.curve) {
            return a.id.curve < b.id.curve;
        }
        if (a.id.idx != b.id.idx) {
            return a.id.idx < b.id.idx;
        }
        if(a.p == b.p) {
            return false;
        }
        return a.t_val < b.t_val;
    }

    std::ostream &operator<<(std::ostream &o, const IntersectionValue &v) {
        return o << "[t:" << v.t_val << ", " << v.p << ", "
                 << "ptype: " << v.id.ptype << ", c:" << v.id.curve << ", i:" << v.id.idx << "]";
    }

    struct Config{
        bool postprocess = true;
        bool postprocess_collinear = true;
        bool fail_on_approx_equal = true;
        double approx_equal_tol = 1e-6;
    };

    template<typename Orient2dFunc = LeftOfLine, typename CollinearFunc = IsCollinear>
    class PolyClip {
    private:
        using LineSegmentIt= std::vector<ContourSegment<Segment>>::const_iterator;
        using CurveSegmentIt = std::vector<ContourSegment<CubicBezier>>::const_iterator;
        enum BooleanOpTypeImpl {
            IMPL_UNION_ = 0,
            IMPL_INTERSECTION_ = 1,
            IMPL_DIFFERENCE_2_ = 2,
            IMPL_DIFFERENCE_ = 3,
            IMPL_XOR_ = 4
        };

        BooleanOpType initial_clippingop;
        LineSegmentIt a_begin, a_end;
        CurveSegmentIt b_begin, b_end;

        std::vector<Contour> &resultpoly;

        BooleanOpTypeImpl used_clippingop_;
        std::set<IntersectionValue> inters_;
        std::vector<ContourSegment<Segment>> lines_{};
        std::vector<ContourSegment<CubicBezier>> curves_{};

        std::deque<SweepPoint> resource_holder_{};
        std::vector<SweepPoint *> queue_{};
        detail::SComp<Orient2dFunc, CollinearFunc> sline_comp_{};
        std::set<SweepPoint *, decltype(sline_comp_)> sline_{sline_comp_};
        using sline_iterator_t = typename decltype(sline_)::iterator;
        std::vector<sline_iterator_t> other_iters_;

        bool bad_ = false;
        Config config_;
    public:
#ifdef DEBUG
        BasicPathWriter *w = nullptr;
        bool connectedges_verbose = false;
        bool initqueue_verbose = false;
        bool computefields_verbose = false;
        bool inresult_verbose = false;
        bool compute_verbose = false;
        bool computecontours_verbose = false;
        bool sweep_verbose = false;
#endif
PolyClip(const std::vector<Contour> &a, const std::vector<Contour> &b,
                                       std::vector<Contour> &result, BooleanOpType clippingop, Config c = {}) :
                initial_clippingop(clippingop), resultpoly(result), config_(c) {
            collect_segments(a, b);
            init_op_enum();
        }

        /// \brief indicates if the clipping operation succeded. If it is false, the Multipolygon this.resultpoly
        /// which was supposed to store the result may contain anything.
        /// \return true if the clipping operation succeded, false otherwise.
        bool success() noexcept {
            return !this->bad_;
        }
        /// \brief computes the clipping operation of the stored input
        void compute() noexcept {
            //phase 1
            inters_ = std::set<IntersectionValue>{};
            if (!intersect_all_segments(lines_.begin(), lines_.end(),
                                    curves_.begin(), curves_.end(), config_, inters_)){
                this->bad_ = true;
                return;
            }

#ifdef DEBUG
            if(compute_verbose) {
                std::cout << "intersections: " << '\n';
                for (const auto &interpoint: inters_) {
                    std::cout << interpoint << '\n';
                }
            }
#endif

            std::size_t num_estimate = 2 * inters_.size() - 2 * lines_.size() + curves_.size();
            queue_.reserve(num_estimate);
            init_queue();
            std::sort(queue_.begin(), queue_.end(), detail::queue_comp<Orient2dFunc, CollinearFunc>);
            other_iters_ = std::vector<sline_iterator_t>(queue_.size());
#ifdef DEBUG
            if(compute_verbose) {
                std::cout << "\n\nsorted the queue\n\n";
            }
#endif
#ifdef DEBUG
            if(compute_verbose) {
                std::cout << "queue:\n";
                for (const auto &item: queue_) {
                    std::cout << *item << '\n';
                }
            }
            if( w != nullptr) {
                for (const auto &curr: queue_) {
                    w->push_circle(curr->point.x(), curr->point.y());
                }
                for (const auto &curr: queue_) {
                    if (curr->left) {
                        w->push_path_str(
                                (!curr->curve ? segment_to_path_str(curr->point, curr->other_point->point)
                                              : bezier_to_path_str(curr->point, curr->controlp, curr->other_point->controlp,
                                                                   curr->other_point->point)),
                                "none", "lightgrey");
                    }
                }
                w->write_to("debug_allpoints_alledges.svg");
            }
                std::cout << std::flush;
#endif
            queue_correctness();
            if (!success()) {
#ifdef DEBUG
                std::cout << "queue correctness failed\n";
#endif
                return;
            }

            //phase 2
            sweep();
            if (!success()) {
#ifdef DEBUG
                std::cout << "sweep failed\n";
#endif
                return;
            }

            connect_edges(0);

            if (initial_clippingop == XOR){
                connect_edges(IMPL_DIFFERENCE_2_);
            }

            if (initial_clippingop == DIVIDE){
                connect_edges(IMPL_DIFFERENCE_2_);
                connect_edges(IMPL_INTERSECTION_);
            }
        }

    private:
        void collect_segments(const std::vector<Contour> &a, const std::vector<Contour> &b) {
            std::size_t contour_idx = 0;
            detail::PolygonType currpolygon = detail::SUBJECT;

            auto add_line_segment = [this, &contour_idx, &currpolygon](const Point2d& a, const Point2d b) {
                bool inc = increasing(a, b);
                if (inc) {
                    this->lines_.push_back({{a, b}, contour_idx, currpolygon});
                } else {
                    this->lines_.push_back({{b, a}, contour_idx, currpolygon});
                }
            };

            auto add_curve_segment = [this, &contour_idx, &currpolygon](const Point2d& p0,
                    const Point2d& p1, const Point2d& p2, const Point2d& p3) {
                bool inc = bezier_direction({p0, p1, p2, p3});
                if (inc) {
                    this->curves_.push_back({{p0, p1, p2, p3}, contour_idx, currpolygon});
                } else {
                    this->curves_.push_back({{p3, p2, p1, p0}, contour_idx, currpolygon});
                }
            };

            for (auto it = a.begin(); it != a.end(); ++it, ++contour_idx) {
                if(it->size() < 1) continue;
                it->template forward_segments<contourklip::LINE>(add_line_segment);
                it->template forward_segments<contourklip::CUBIC_BEZIER>(add_curve_segment);
                //the input should be const, hence we don't actually close it.
                if (!it->is_closed()) {
                    add_line_segment(it->front_point(), it->back().point());
                }
            }

            currpolygon = detail::CLIPPING;
            contour_idx = 0;
            for (auto it = b.begin(); it != b.end(); ++it, ++contour_idx) {
                if(it->size() < 1) continue;
                it->template forward_segments<contourklip::LINE>(add_line_segment);
                it->template forward_segments<contourklip::CUBIC_BEZIER>(add_curve_segment);
                if (!it->is_closed()) {
                    add_line_segment(it->front_point(), it->back().point());
                }
            }
        }

        bool intersect_all_segments(
                const LineSegmentIt &lines_begin, const LineSegmentIt &lines_end,
                const CurveSegmentIt &curves_begin, const CurveSegmentIt &curves_end,
                const Config &config,
                std::set<IntersectionValue> &out) {
            bool success = true;
            // lines against lines
            std::size_t i = 0;
            auto add = [&out](const IntersectionValue &t) {
                out.insert(t);
            };
            for (auto it = lines_begin; it != lines_end; ++it, ++i) {
                IntersectionValue curr_start{0., (*it).seg.first, {(*it).ptype, false, i}};
                IntersectionValue curr_end{1., (*it).seg.second, {(*it).ptype, false, i}};
                add(curr_start);
                add(curr_end);

                std::size_t j = i + 1;
                auto jt = it;
                ++jt;
                for (; jt != lines_end; ++jt, ++j) {
                    if (auto val = intersect_segments((*it).seg, (*jt).seg)) {
                        Point2d mapped = val->p;
                        double t1 = val->t1;
                        double t2 = val->t2;

                        if (t1 > 0 && t1 < 1.) {
                            SegIndent curr{(*it).ptype, false, i};
                            add(IntersectionValue{t1, mapped, curr});
                        }
                        if (t2 > 0 && t2 < 1.) {
                            SegIndent curr{(*jt).ptype, false, j};
                            add(IntersectionValue{t2, mapped, curr});
                        }
                    }
                }
            }
            i = 0;
            auto a_it = lines_begin;
            // curves against lines
            for (auto it = curves_begin; it != curves_end; ++it, ++i, ++a_it) {
                IntersectionValue curr_start{0., (*it).seg.p0, {(*it).ptype, true, i}};
                IntersectionValue curr_end{1., (*it).seg.p3, {(*it).ptype, true, i}};
                add(curr_start);
                add(curr_end);
                std::size_t j = 0;
                for (auto jt = lines_begin; jt != lines_end; ++jt, ++j) {
                    auto add_val = [&](double t, double u, Point2d p) {
                        bool parametric_ok = !std::isnan(u) || !std::isinf(u);
                        bool intersection_ok = approx_equal(beziermap((*it).seg, u),
                                                            linear_map((*jt).seg, t),
                                                            1e-4);
#ifdef DEBUG
                        if (std::isnan(u) || std::isinf(u)) {
                            std::cout << "t values not valid\n";
                        }
                        if (!intersection_ok) {
                            std::cout << "issue with intersection point\n";
                            std::cout << (*it).seg << '\n' << (*jt).seg
                                      << " " << t << " " << u << " " << p << '\n';
                            assert(false);
                        }
#endif
                        success = parametric_ok && intersection_ok;
                        IntersectionValue curr_line{t, p, {(*jt).ptype, false, j}};
                        IntersectionValue curr_curve{u, p, {(*it).ptype, true, i}};
                        add(curr_line);
                        add(curr_curve);
                    };
                    line_bezier_inter((*jt).seg, (*it).seg, add_val);
                    if (!success) return false;
                }
            }
            i = 0;
            // curves against curves
            for (auto it = curves_begin; it != curves_end; ++it, ++i) {
                std::size_t j = i + 1;
                auto jt = it;
                ++jt;
                for (; jt != curves_end; ++jt, ++j) {
                    auto add_val = [&](double t, double u, Point2d p) {
                        Point2d interp = beziermap((*it).seg, t);
                        Point2d interp2 = beziermap((*jt).seg, u);
                        bool parametric_ok = !std::isnan(u) && !std::isinf(u);
                        bool intersection_ok = approx_equal(interp, interp2, 1e-4);
#ifdef DEBUG
                        if (!approx_equal(interp, interp2, 1e-4)) {
                            std::cout << "issue with intersection point\n";
    //                        std::cout << curves.at(i).seg << '\n' << curves.at(j).seg
    //                                  << " " << t << " " << u << '\n';
                            std::abort();
                        }
#endif
                        success = parametric_ok && intersection_ok;
                        // we use inter. point of first curve
                        IntersectionValue a{t, p, {(*it).ptype, true, i}};
                        IntersectionValue b{u, p, {(*jt).ptype, true, j}};
                        add(a);
                        add(b);
                    };
                    curve_curve_inter<decltype(add_val), double, CollinearFunc>((*it).seg, (*jt).seg, add_val);
                    if (!success) return false;
                }
            }
            return true;
        }

        void init_queue() noexcept {
            std::size_t k = 1;
            auto prev_it = inters_.begin();
            auto it = prev_it;
            ++it;
            while (it != inters_.end() && k < inters_.size()) {
                while (it != inters_.end() && k < inters_.size() && (*prev_it).id == (*it).id) {
#ifdef DEBUG
                    //                std::cout << "added(\n" << inters_[k - 1] << ",\n" << inters_[k] << "\n)";
#endif
                    IntersectionValue prev = (*prev_it);
                    IntersectionValue curr = (*it);
                    if (!curr.id.curve) {
                        // although the intersections are sorted with t values,
                        // it does not guarantee that the second is on the right of the first.
                        if (!increasing(prev.p, curr.p)) {
                            std::swap(prev, curr);
                        }
                        auto *a = add_sweep_point(true, prev.p, nullptr);
                        auto *b = add_sweep_point(false, curr.p, a);
                        //setting the correct fields
                        a->other_point = b;
                        a->ptype = prev.id.ptype;
                        b->ptype = curr.id.ptype;
                        a->contourid = b->contourid = lines_[curr.id.idx].contourid;
                        //setting to some useful value. this is important when comparing sweeppoints.
                        a->controlp = b->point;
                        b->controlp = a->point;
                        queue_.push_back(a);
                        queue_.push_back(b);
                    } else {
                        // we need to split curves into monotonic segments
#ifdef DEBUG
                        if(initqueue_verbose) {
                            std::cout << "splitting curve " << curves_[curr.id.idx].seg
                                      << "\n    at " << prev << " " << curr << '\n';
                        }
#endif
                        auto maybebezier = sub_bezier(curves_[curr.id.idx].seg,
                                                      prev.t_val, curr.t_val);
                        if (!maybebezier) {
                            this->bad_ = true;
                            return;
                        }
                        CubicBezier currbezier = *maybebezier;
                        // this is for consistency
                        currbezier.p0 = prev.p;
                        currbezier.p3 = curr.p;

                        double t_prev = 0.;
                        SweepPoint *first = nullptr;

                        Extremity_Direction d_prev;

                        auto process_f = [&](double t, Extremity_Direction d) {
                            auto maybemonoton = sub_bezier(currbezier, t_prev, t);
                            if (!maybemonoton) {
                                this->bad_ = true;
                                return;
                            }
                            CubicBezier monoton = *maybemonoton;
                            auto *a = add_sweep_point(false, monoton.p0, nullptr);
                            auto *b = add_sweep_point(false, monoton.p3, a);
                            a->controlp = monoton.p1;
                            b->controlp = monoton.p2;

                            // consistency regarding inner control points. tangents are either horizontal or vertical.
                            if (t != 1.) {
                                if (d == X_Extremity) {
                                    b->controlp = {b->point.x(), b->controlp.y()};
                                } else if (d == Y_Extremity) {
                                    b->controlp = {b->controlp.x(), b->point.y()};
                                }
                            }
                            if (t_prev != 0.) {
                                if (d_prev == X_Extremity) {
                                    a->controlp = {a->point.x(), a->controlp.y()};
                                } else if (d_prev == Y_Extremity) {
                                    a->controlp = {a->controlp.x(), a->point.y()};
                                }
                            }
                            a->initial_controlp = currbezier.p1;
                            b->initial_controlp = currbezier.p2;

                            a->curve = b->curve = true;
                            a->other_point = b;
                            a->ptype = b->ptype = curr.id.ptype;
                            a->contourid = b->contourid = curves_[curr.id.idx].contourid;
                            a->set_if_left();

                            if (t_prev == 0.) {
                                a->islast = true;
                                first = a;
                            }
                            if (t == 1.) {
                                b->islast = true;
                            }
                            a->start = first;
                            b->start = first;
                            queue_.push_back(a);
                            queue_.push_back(b);
#ifdef DEBUG
                            assert(a->point != a->other_point->point);
                            assert(b->point != b->other_point->point);
#endif
                            t_prev = t;
                            d_prev = d;
                        };
                        bezier_monotonic_split(currbezier, process_f);
                        // last segment
                        process_f(1, X_Extremity);
                    }
                    it++;
                    prev_it++;
                    k++;
                }
                if (it == inters_.end()) {
                    break;
                }
                it++;
                prev_it++;
                k++;
            }
        }

        SweepPoint *add_sweep_point(bool left, const Point2d &p, SweepPoint *other) {
            return &resource_holder_.emplace_back(left, p, other);
        }

        void compute_fields(SweepPoint *curr,
                            const sline_iterator_t& prev) {
            if (prev != sline_.end()) {
                if (curr->ptype == (*prev)->ptype) {
#ifdef DEBUG
                    //                std::cout << "same polygon ";
#endif
                    curr->in_out = !(*prev)->in_out;
                    curr->other_in_out = (*prev)->other_in_out;

                } else {
#ifdef DEBUG
                    //                std::cout << "different polygon ";
#endif
                    curr->in_out = !(*prev)->other_in_out;
                    curr->other_in_out = (*prev)->vertical() ? !(*prev)->in_out : (*prev)->in_out;

                    if (overlapping(curr, (*prev))) {
#ifdef DEBUG
                        std::cout << "dupe edge ";
#endif
                        if (curr->in_out == (*prev)->in_out) {
                            curr->edgetype = detail::SAME_TRANSITION;
                            (*prev)->edgetype = detail::SAME_TRANSITION;
                        } else {
                            curr->edgetype = detail::DIFFERENT_TRANSITION;
                            (*prev)->edgetype = detail::DIFFERENT_TRANSITION;
                        }
                        // a duplicate edge is added at most once.
                        // Therefore, we set the previous in result to false, and it is then
                        // decided if the current is in the result.
#ifdef DEBUG
                        //                        std::cout << "previous in result, now " <<(*prev)->in_result ;
#endif

                        (*prev)->in_result[0] = false;
                        (*prev)->in_result[IMPL_INTERSECTION_] = false;
                        (*prev)->in_result[IMPL_DIFFERENCE_2_] = false;
                    }
                }

                curr->prev_in_result[0] = (!(*prev)->in_result[0]
                                        || (*prev)->vertical()) ? (*prev)->prev_in_result[0] : *prev;
                curr->prev_in_result[IMPL_INTERSECTION_] = (!(*prev)->in_result[IMPL_INTERSECTION_]
                                        || (*prev)->vertical()) ? (*prev)->prev_in_result[IMPL_INTERSECTION_] : *prev;
                curr->prev_in_result[IMPL_DIFFERENCE_2_] = (!(*prev)->in_result[IMPL_DIFFERENCE_2_]
                                        || (*prev)->vertical()) ? (*prev)->prev_in_result[IMPL_DIFFERENCE_2_] : *prev;
            }
                //special case: no predecessor
            else {
                curr->in_out = false;
                curr->other_in_out = true;
            }
#ifdef DEBUG
            //        std::cout << "computed fields inout, otherinout " << curr->inOut << " " << curr->other_in_out;
#endif
            curr->in_result[0] = in_result(curr, used_clippingop_);
            curr->in_result[IMPL_INTERSECTION_] = in_result(curr, IMPL_INTERSECTION_);
            curr->in_result[IMPL_DIFFERENCE_2_] = in_result(curr, IMPL_DIFFERENCE_2_);
        }

        bool in_result(SweepPoint *curr, PolyClip::BooleanOpTypeImpl clippingop) {
            // We have to check if the edge lies inside or outside the other polygon.
            // Then, it depends on the boolean operation:
            //
            // intersection:    we select the edges which are inside the other polygon.
            // union:           we select the edges which are outside the other polygon.
            // subtract A - B:  we select the outside edges from A and the inside edges from B.
            // xor:             every edge is in the result.
            //
            // If the closest other edge is an outside inside transition, then current edge is inside
            // the other polygon, otherwise outside. We also have to treat the special case with duplicate edges.

            switch (curr->edgetype) {
                case detail::NORMAL:
                    switch (clippingop) {
                        case IMPL_INTERSECTION_:
                            return !curr->other_in_out;
                        case IMPL_UNION_:
                            return curr->other_in_out;
                        case IMPL_DIFFERENCE_:
                            return (curr->ptype == detail::SUBJECT && curr->other_in_out) ||
                                   (curr->ptype == detail::CLIPPING && !curr->other_in_out);
                        case IMPL_DIFFERENCE_2_:
                            return (curr->ptype == detail::SUBJECT && !curr->other_in_out) ||
                                   (curr->ptype == detail::CLIPPING && curr->other_in_out);
                        case IMPL_XOR_:
                            return true;
                    }
                case detail::SAME_TRANSITION:
                    return clippingop == IMPL_INTERSECTION_ || clippingop == IMPL_UNION_;
                case detail::DIFFERENT_TRANSITION:
                    return clippingop == IMPL_DIFFERENCE_ || clippingop == IMPL_DIFFERENCE_2_;
                default:
                    return false;
            }
            return true;
        }

        void connect_edges(int optype_idx=0) {
            std::vector<SweepPoint *> sorted_result{};

            int k = 0;
            for (auto sp: queue_) {
                if (sp->in_result[optype_idx]) {
                    sorted_result.push_back(sp);
                    sp->pos = k;
                    k++;
                }
            }

            if (sorted_result.size() <= 2) return;
#ifdef DEBUG
            if(connectedges_verbose) {
                std::cout << "sorted result: \n";
                for (const auto &i: sorted_result) {
                    std::cout << *i << '\n';
                }
                std::cout << '\n';
            }
#endif
            std::vector<bool> processed(sorted_result.size(), false);
            std::vector<std::size_t> depth{};
            std::size_t contourId;
            for (std::size_t i = 0; i < sorted_result.size(); ++i) {
                if (processed[i]) {
                    continue;
                }
#ifdef DEBUG
                if(connectedges_verbose) {
                    std::cout << "\n--------computing new contour, starting with " << *sorted_result[i] << '\n';
                }
#endif
                contourId = depth.size();
                depth.push_back(0); // we first assume depth[result_contour_id] is an outer contour
                Contour c{};

                compute_contour(sorted_result[i], c, contourId, processed, sorted_result);
                if (!success()) { return; }

                if (c.size() < 3) {
#ifdef DEBUG
                    if(connectedges_verbose) {
                        std::cout << "contour is only 2 points";
                    }
#endif
                    continue;
                }
                SweepPoint* previnresult = sorted_result[i]->prev_in_result[optype_idx];
                if (previnresult) {
                    std::size_t lowercontourId = previnresult->result_contour_id;
#ifdef DEBUG
                    assert(lowercontourId < depth.size() && "issue with result_contour_id");
#endif
                    if (!previnresult->result_in_out) {
                        depth[contourId] = depth[lowercontourId] + 1;
                    }
                }

                bool reverse = depth[contourId] % 2 == 1;
#ifdef DEBUG
                if(connectedges_verbose) {
                    std::cout << "resulting contour\n";
                    std::cout << "result_contour_id " << contourId << " depth of contour " << depth[contourId] << '\n';
                    std::cout << "reverse " << reverse << " contourid " << contourId << " depth " << depth[contourId]
                              << '\n';
                }
                assert(contourId < depth.size());
                assert(depth[contourId] < queue_.size());
#endif
                if (config_.postprocess) {
                    auto add = [&](Contour &t) {
                        add_contour(t, ((contour_area(t) > 0) == reverse));
                    };
                    postprocess_contour<decltype(add), CollinearFunc>(c, add, config_.postprocess_collinear);
                } else {
                    add_contour(c, reverse);
                }
            }
        }

        void compute_contour(SweepPoint *sp, Contour &c, std::size_t contourId,
                             std::vector<bool> &processed,
                             const std::vector<SweepPoint *> &sorted_result) {
            std::size_t currpos;
            Point2d startpoint = sp->point;
            c.push_back(startpoint);
            processed[sp->pos] = true;
            sp->result_contour_id = contourId;
            sp->other_point->result_contour_id = contourId;

            currpos = sp->other_point->pos;
            CONTOURKLIP_IF_NOT(currpos > 0 && currpos < sorted_result.size(), "currpos should be valid index") {
                this->bad_ = true;
                return;
            }
            processed[currpos] = true;
            //note that we are necessarily starting with the lowest left edge
            sp->result_in_out = false;
            sp->other_point->result_in_out = false;

            std::size_t k = 0;
            bool cycle = false;

            while (sorted_result[currpos]->point != startpoint && !cycle && k < sorted_result.size()) {
                if (sorted_result[currpos]->curve) {
                    if (sorted_result[currpos]->islast) {
                        c.push_back(sorted_result[currpos]->other_point->initial_controlp,
                                    sorted_result[currpos]->initial_controlp, sorted_result[currpos]->point);
                    }
#ifdef DEBUG
                    if(computecontours_verbose) {
                        std::cout << "added a curve\n";
                    }
#endif
                } else {
                    c.push_back(sorted_result[currpos]->point);
                }
                sorted_result[currpos]->result_contour_id = contourId;
                sorted_result[currpos]->other_point->result_contour_id = contourId;
                //if we're traversing from left to right, we're at a right point, hence we
                // don't have an inout transition.
                sorted_result[currpos]->result_in_out
                        = sorted_result[currpos]->other_point->result_in_out = sorted_result[currpos]->left;

#ifdef DEBUG
                if(computecontours_verbose) {
                    std::cout << "added idx:" << currpos << " " << sorted_result[currpos]->point << '\n';
                    if (w != nullptr) {
                        w->push_path_str(path_to_str(c, false), "none", "black", 2);
                        w->write_to("debug_contour_" + std::to_string(contourId) + "_" + std::to_string(k) + ".svg");
                    }
                }
#endif
                std::size_t npos;
                if (auto maybe_npos = nextpos(currpos, processed, sorted_result)) {
                    npos = sorted_result[*maybe_npos]->other_point->pos;
                } else {
#ifdef DEBUG
                    std::cout << "issue with nextpos\n";
                    assert(false);
#endif
                    this->bad_ = true;
                    return;
                }

                // we need to detect cycles, but we cannot break because
                // we have to properly set the processed[] entries.
                cycle = (npos == currpos);
                currpos = npos;
#ifdef DEBUG
                if(computecontours_verbose) {
                    std::cout << '\n';
                    if (cycle) {
                        std::cout << "2-cycle detected\n";
                    }
                }
#endif
                processed[sorted_result[currpos]->pos] = true;
                //if we have processed left, we also have processed right
                processed[sorted_result[currpos]->other_point->pos] = true;
#ifdef DEBUG
                if(computecontours_verbose) {
                    std::cout << "\nnext position: " << currpos;
                    std::cout << "\n\n";
                }
#endif
                k++;
            }
            sorted_result[currpos]->result_in_out
                    = sorted_result[currpos]->other_point->result_in_out = true;

            // close the contour.
            if (sorted_result[currpos]->curve) {
                //special case: the first and last components (ie their monotonic segments) are from the same curve
                //since it has already been included, we only need to change the starting point
                if (sp->curve && sp->start == sorted_result[currpos]->start) {
#ifdef DEBUG
                    if(computecontours_verbose) {
                        std::cout << "last curve is same as start\n";
                        std::cout << c.front().point() << " " << sp->start->point << " "
                                  << sorted_result[currpos]->point << " " << sorted_result[currpos]->other_point->point
                                  << '\n';
                    }
#endif

                    c.front().point() = c.back().point();
                } else {
                    c.push_back(sorted_result[currpos]->other_point->initial_controlp,
                                sorted_result[currpos]->initial_controlp, sorted_result[currpos]->point);
                }
            } else {
                c.push_back(sorted_result[currpos]->point);
            }
#ifdef DEBUG
            if(computecontours_verbose) {
                std::cout << "added idx:" << currpos << " " << sorted_result[currpos]->point << '\n';
                if (w != nullptr) {
                    w->push_path_str(path_to_str(c, false), "none", "black", 2);
                    w->write_to("debug_contour_" + std::to_string(contourId) + "_" + std::to_string(k) + ".svg");
                }
            }
#endif
        }

        std::optional<std::size_t> nextpos(std::size_t currpos, const std::vector<bool> &processed,
                                                            const std::vector<SweepPoint *> &sorted_result) {
            //we first try to find one which is from another contour
            std::size_t npos = currpos;
            npos++;

            while (npos < sorted_result.size()
                   && sorted_result.at(npos)->point == sorted_result.at(currpos)->point) {
                if (!processed.at(npos)) {
                    return npos;
                }
                npos++;
            }
            npos = currpos;
            npos--;

            while (npos > 0
                   && sorted_result.at(npos)->point ==
                      sorted_result.at(currpos)->point) {
                if (!processed.at(npos)) {
                    return npos;
                }
                npos--;
            }

            npos = currpos;
            npos++;
            while (npos < sorted_result.size()
                   && sorted_result.at(npos)->point == sorted_result.at(currpos)->point) {
                if (!processed.at(npos)) {
                    return npos;
                }
                npos++;
            }
            npos = currpos;
            npos--;
            while (npos > 0
                   && sorted_result.at(npos)->point == sorted_result.at(currpos)->point) {
                if (!processed.at(npos)) {
                    return npos;
                }
                npos--;
            }

            CONTOURKLIP_IF_NOT(true, "issue finding next position") {
                this->bad_ = true;
            }

            return {};
        }

        void add_contour(Contour &c, bool reverse) {
            if (reverse) {
                c.reverse();
            }
            resultpoly.push_back(c);
        }

        void queue_correctness() {
            int count = 1;
            for (std::size_t i = 1; i < queue_.size(); ++i) {
                CONTOURKLIP_IF_NOT(queue_.at(i)->point.x() >= queue_.at(i - 1)->point.x(),
                                 "queue_ needs to be sorted") {
                    this->bad_ = true;
                    return;
                }

                if (config_.fail_on_approx_equal) {
                    CONTOURKLIP_IF_NOT(!approx_equal(queue_.at(i)->point, queue_.at(i)->other_point->point,
                                                     config_.approx_equal_tol),
                                     "other point should not be approx. equal to current point") {
                        this->bad_ = true;
                        return;
                    }

                    CONTOURKLIP_IF_NOT((queue_.at(i)->point == queue_.at(i - 1)->point
                                      || !approx_equal(queue_.at(i)->point, queue_.at(i - 1)->point,
                                                       config_.approx_equal_tol)),
                                     "points should not be almost the same, otherwise possible numerical issue"
                    ) {
                        this->bad_ = true;
                        return;
                    }
                }

                CONTOURKLIP_IF_NOT(queue_.at(i)->point.x() >= queue_.at(i)->other_point->point.x() || queue_.at(i)->left,
                                 "right point of the same segment comes after left point") {
                    this->bad_ = true;
                    return;
                }
                if (queue_.at(i)->point == queue_.at(i - 1)->point) {
                    CONTOURKLIP_IF_NOT(((!queue_.at(i - 1)->left) || queue_.at(i)->left),
                                     "right points should come before left points"
                    ) {
                        this->bad_ = true;
                        return;
                    }
                    count++;
                } else {
                    CONTOURKLIP_IF_NOT(count % 2 == 0, "same points should occur an even number of times") {
                        this->bad_ = true;
                        return;
                    }
                    count = 1;
                }
            }
        }

        void sweep() {
            for (std::size_t k = 0; k < queue_.size(); ++k) {
                SweepPoint *curr = queue_.at(k);
#ifdef DEBUG
                if(sweep_verbose) {
                    std::cout << '\n' << k << "----------------------------------------\n";
                    std::cout << "sline_ before insert: ";
                    for (const auto &i: sline_) {
                        std::cout << (*i) << '\n';
                    }
                    std::cout << "\n\n";
                    std::cout << " curr " << *curr << '\n' << std::flush;
                }
#endif

                if (curr->left) {
                    auto it = sline_.insert(curr);

                    if (!it.second) {
#ifdef DEBUG
                        std::cerr << "the insertion must work since there shouldn't be duplicates\n";
                        std::cerr << "curr: " << *curr << '\n';
                        assert(false);
#endif
                        this->bad_ = true;
                        return;
                    }

                    other_iters_[k] = it.first;
                    curr->other_point->other_iter_pos = k;
                    auto prev = sline_.end();
                    //there's a valid predecessor
                    if (it.first != sline_.begin()) {
                        prev = it.first;
                        prev--;
#ifdef DEBUG
                        //                    std::cout << "predecessor is" << (*prev)->point << ", " << (*prev)->other_point->point << std::flush;
#endif
#ifdef DEBUG
                        if (curr->point == (*prev)->point && curr->other_point->point == (*prev)->other_point->point) {
                            //                        std::cout << "duplicate edge " << std::flush;
                        }
#endif
                    }
                    compute_fields(curr, prev);
#ifdef DEBUG
                    if (curr->in_result[0]) {
                        if(w != nullptr) {
                            if(sweep_verbose) {
                                std::cout << k << " in result_candidates " << std::flush;
                            }
                            w->push_path_str((!curr->curve
                                              ? segment_to_path_str(curr->point, curr->other_point->point)
                                              : bezier_to_path_str(curr->point, curr->controlp,
                                                                   curr->other_point->controlp, curr->other_point->point)),
                                             "none", "#555555");
                            w->write_to("debug_sweep" + std::to_string(k) + ".svg");
                    }
                    }
#endif
                } else {
                    auto it2 = other_iters_[curr->other_iter_pos];
                    auto it = sline_.find(curr->other_point);
                    if (it != it2) {
#ifdef DEBUG
                        std::cerr << "the other iterator must be present\n";
                        std::cerr << "curr: " << *curr << '\n';
                        std::cerr << "other: " << *(curr->other_point) << '\n';
                        std::abort();
#endif
                        this->bad_ = true;
                        return;
                    }
                    sline_.erase(it);

                    curr->in_result[0] = curr->other_point->in_result[0];
                    curr->in_result[IMPL_INTERSECTION_] = curr->other_point->in_result[IMPL_INTERSECTION_];
                    curr->in_result[IMPL_DIFFERENCE_2_] = curr->other_point->in_result[IMPL_DIFFERENCE_2_];
                }
            }
        }

        void init_op_enum() {
            switch (initial_clippingop) {
                case INTERSECTION:
                    used_clippingop_ = IMPL_INTERSECTION_;
                    break;
                case UNION:
                    used_clippingop_ = IMPL_UNION_;
                    break;
                case DIFFERENCE:
                    used_clippingop_ = IMPL_DIFFERENCE_;
                    break;
                case XOR:
                    used_clippingop_ = IMPL_DIFFERENCE_;
                    break;
                case DIVIDE:
                    used_clippingop_ = IMPL_DIFFERENCE_;
                    break;
            }
        }
    };

#undef CONTOURKLIP_IF_NOT

    bool clip(const std::vector<Contour> &a, const std::vector<Contour> &b,
              std::vector<Contour> &result, BooleanOpType clippingop) {
        PolyClip c(a, b, result, clippingop);
        c.compute();
        return c.success();
    }

    bool clip(const Contour &a, const Contour &b,
              std::vector<Contour> &result, BooleanOpType clippingop) {
        std::vector poly1{a}, poly2{b};
        PolyClip c(poly1, poly2, result, clippingop);
        c.compute();
        return c.success();
    }
}
#endif //CONTOURKLIP_POLYCLIP_HPP