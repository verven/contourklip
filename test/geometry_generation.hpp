#ifndef CONTOURKLIP_GEOMETRY_GENERATION_HPP
#define CONTOURKLIP_GEOMETRY_GENERATION_HPP

#include <random>
#include "geometry_base.hpp"

namespace geometrygen{
    struct UniformF{
        std::mt19937 gen;
        std::uniform_real_distribution<> d;
        UniformF(double a, double b, int seed = 1): d(a, b){
            gen =  std::mt19937(seed);
        }
        double operator()(void){
            return d(gen);
        }
    };

    template<typename P>
    class PointGenerator{
        UniformF h, v;
     public:
        PointGenerator(double xa=0, double xb=1, double ya=0,
                       double yb=0, int seed=1) : h(xa, xb, seed), v(ya, yb, seed+1) {}

        P operator()(void){
            return P{h(), v()};
        }
    };

    std::vector<double> random_partition(double a, double b, int n, int seed=0){
        UniformF dis(a, b, seed);

        std::vector<double> v;
        for (int i = 0; i < n-1; ++i) {
            v.push_back(dis());
        }
        std::sort(v.begin(), v.end());
        return v;
    }

    std::vector<double> random_partition2(double a, double b, int n, int seed=0){
        double interval = (b-a)/(double(n));

        UniformF dis(0, interval, seed);

        std::vector<double> v{};
        v.reserve(n-1);
        for (int i = 1; i < n; ++i) {
            v.push_back(dis() + double(i*interval));
        }
        return v;
    }

    void generate_contour(int n, double a, double b, double c, int seed, contourklip::Contour& out, contourklip::Point2d offset = {0, 0}, bool withcurves = true){
        UniformF rand_dist(a, b, seed);
        UniformF rand_offset(std::max(0., b-c), b+c, seed+1);

        auto partition = random_partition2(0, 2*M_PI, n, seed+2);
        double prev_angle =0;
        double d = rand_dist();
        out.push_back({d + offset.x(), 0+offset.y()});

        auto frompolar = [&offset](double r, double angle) -> contourklip::Point2d{
            return {r*std::cos(angle) + offset.x(), r * std::sin(angle) + offset.y()};
        };
        int k =0;
        for (const auto &angle: partition) {
            d = rand_dist();
            if (withcurves && (rand_dist() < a + 0.5*(b-a))){
                UniformF rand_angle(prev_angle, angle, seed);
                double a1 = rand_angle(), a2 = rand_angle();
                if(a1 > a2) std::swap(a1, a2);
                double d1 = rand_offset(), d2 = rand_offset();
                out.push_back(frompolar(d1, a1),
                              frompolar(d2, a2), frompolar(d, angle));
            }else{
                out.push_back(frompolar(d, angle));
            }
            prev_angle = angle;
            k++;
        }
        out.close();
    }

    void generate_rectangle(const contourklip::Point2d& a, const contourklip::Point2d& b, contourklip::Contour& out){
        out.push_back(a);
        out.push_back({a.x(), b.y()});
        out.push_back(b);
        out.push_back({b.x(), a.y()});
        out.close();
    }

};
#endif //CONTOURKLIP_GEOMETRY_GENERATION_HPP