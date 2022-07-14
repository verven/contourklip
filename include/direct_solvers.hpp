#ifndef CONTOURKLIP_DIRECT_SOLVERS_HPP
#define CONTOURKLIP_DIRECT_SOLVERS_HPP
#include <cmath>
namespace directsolvers {
    double clip(double val, double lower, double upper) {
        return std::max(lower, std::min(val, upper));
    }

    /*
      see: https://stackoverflow.com/questions/63665010

      diff_of_products() computes a*b-c*d with a maximum error <= 1.5 ulp

      Claude-Pierre Jeannerod, Nicolas Louvet, and Jean-Michel Muller,
      "Further Analysis of Kahan's Algorithm for the Accurate Computation
      of 2x2 Determinants". Mathematics of Computation, Vol. 82, No. 284,
      Oct. 2013, pp. 2245-2264
    */
    double diff_of_products (double a, double b, double c, double d) {
        double w = d * c;
        double e = std::fma (-d, c, w);
        double f = std::fma (a, b, -w);
        return f + e;
    }


    int solve_quadratic(double a_0, double a_1, double a_2, std::pair<double, double> &r) {
        if (std::abs(a_2) < 1e-15) {
            if (std::abs(a_1) < 1e-15) {
                return 0;
            }
            r.first = r.second = -a_0 / a_1;
            return 1;
        }
        double d = diff_of_products(a_1, a_1, 4.0*a_2, a_0);
        if (d < 0) {
            return 0;
        }
        double sqd = sqrt(d);
        double u = 1.0 / a_2;
        if (a_1 >= 0.0) {
            double t = 0.5 * (-a_1 - sqd) * u;
            r.first = t;
            r.second = u * a_0 / t;
        } else {
            double t = 0.5 * (-a_1 + sqd) * u;
            r.first = u * a_0 / t;
            r.second = t;
        }
        return 2;
    }

    template<typename RootConsumer>
    void solve_cubic_real(double a_0, double a_1, double a_2, double a_3, RootConsumer &c, double tol) {
        //special case: not a cubic, fall back to quadratic
        if (std::abs(a_3) < tol) {
            std::pair<double, double> r{};
            int t = solve_quadratic(a_0, a_1, a_2, r);
            if (t == 1) {
                c(r.first);
            }
            if (t == 2) {
                c(r.first);
                c(r.second);
            }
            return;
        }
        //normalize so that highest coefficient is 1
        a_2 /= a_3;
        a_1 /= a_3;
        a_0 /= a_3;

        double a2_2 = a_2 * a_2;
        double a_2over3 = a_2 / 3.;
        double q = a_1 / 3.0 - a2_2 / 9.0;
        double r = (a_1 * a_2 - 3. * a_0) / 6.0 - (a2_2 * a_2) / 27.0;
        double rr = r * r;
        double q3 = q * q * q;
        double check = rr + q3;
        // case: three real solutions
        if (check <= 0 || check < tol) {
            double theta = 0;
            if (!(abs(q) < tol)) {
                double temp = clip(r / sqrt(-q3), -1, 1);
                theta = acos(temp);
            }
            double angle1 = theta / 3.;
            double angle2 = angle1 - 2. * M_PI / 3.;
            double angle3 = angle1 + 2. * M_PI / 3.;
            double sq = 2. * std::sqrt(-q);
            double r1, r2, r3;
            r1 = sq * cos(angle3) - a_2over3;
            r2 = sq * cos(angle2) - a_2over3;
            r3 = sq * cos(angle1) - a_2over3;
            // it holds that r1 <= r2 <= r3
            c(r1);
            c(r2);
            c(r3);
            return;

        }
        //only one real solution
        double sq = sqrt(check);
        double u = cbrt(r + sq);
        double v = cbrt(r - sq);
        double r1 = u + v - a_2over3;
        c(r1);
    }
}
#endif //CONTOURKLIP_DIRECT_SOLVERS_HPP