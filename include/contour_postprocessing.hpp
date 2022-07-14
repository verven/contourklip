#ifndef CONTOURKLIP_CONTOUR_POSTPROCESSING_HPP
#define CONTOURKLIP_CONTOUR_POSTPROCESSING_HPP

#include <map>
#include "geometry_base.hpp"

namespace contourklip::detail {
        template<typename outF, typename CollinearF = IsCollinear>
        void postprocess_contour(Contour &c, outF &report_contour, bool remove_collinear = true,
                                 CollinearF collinear_f = {}) {
            if (c.size() <= 3) {
                report_contour(c);
                return;
            }

            std::map<Point2d, std::size_t> visited_p{};
            std::vector<std::size_t> skip(c.size(), 0);
            std::fill(skip.begin(), skip.end(), 0);

            std::size_t prev;
            for (std::size_t i = 0; i < c.size(); ++i) {
                auto v = visited_p.find(c[i].point());
                if (v != visited_p.end()) {
                    prev = v->second;
                    // update the value to be the last updated
                    v->second = i;
                    // it could also be that we have an overl. point of a sub-contour that is already deleted.
                    if (skip[prev] != 0) {
                        continue;
                    }

                    Contour subcontour{};
                    //we still need to keep one of the 2 points on the contour
                    for (std::size_t j = prev; j < i; ++j) {
                        if (remove_collinear && subcontour.size() >= 2
                            && !subcontour.back().bcurve() && !c[j].bcurve()
                            && collinear_f(subcontour[subcontour.size() - 2].point(),
                                           subcontour.back_point(),
                                           c[j].point())
                                ) {
                            subcontour.back().point() = c[j].point();
                        } else {
                            subcontour.push_back(c[j]);
                        }
                        if (skip[j]) {
                            j = skip[j];
                        }
                        skip[j] = i;
                    }
                    subcontour.push_back(c[i]);
                    if (subcontour.size() <= 2) {
                        continue;
                    }
                    report_contour(subcontour);
                }
                    // current point is new
                else {
                    visited_p.insert({c[i].point(), i});
                }
            }
        }
    }
#endif //CONTOURKLIP_CONTOUR_POSTPROCESSING_HPP