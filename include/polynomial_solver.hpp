#ifndef CONTOURKLIP_POLYNOMIAL_SOLVER_HPP
#define CONTOURKLIP_POLYNOMIAL_SOLVER_HPP
#include <array>
#include <cmath>
#include <optional>
namespace polynomialsolver {

    template<typename U, typename T>
    std::optional<T> itp_root_refine(U &func, T a, T b, T eps, int maxiter) {
        //leveraging argument dependent lookup
        using std::log2;
        using std::abs;
        using std::exp2;
        T a_start = a;
        T b_start = b;
        T k1 = (T) 0.2 / (b - a);
        int nmax =  int(log2((b - a) / (2 * eps))) + 2;
        int i = 0;
        T fa = func(a);
        T fb = func(b);
        while (b - a > 2 * eps && i < maxiter) {
            if(fa ==(T)0){
                return a;
            }
            if(fb==(T)0){
                return b;
            }
            //safety check in case interval is degenerate
            if(a < a_start || b > b_start) return {};
            if(fa * fb >0) return {};
            T x_mid = (a + b) / 2.0;
            T r = eps * (exp2(nmax - i)) - (b - a) / 2.0;
            T delta = k1 * (b - a) * (b - a);
            T x_f = (fb * a - fa * b) / (fb - fa);
            T si = x_mid - x_f;
            si = si < 0 ? -1 : 1;
            T x_t = (delta <= abs(x_mid - x_f)) ? x_f + si * delta : x_mid;
            T x_itp = (abs(x_t - x_mid) <= r) ? x_t : x_mid - si * r;
            T f_x = func(x_itp);
            if (f_x * fa < 0) {
                b = x_itp;
                fb = f_x;
            } else {
                a = x_itp;
                fa = f_x;
            }
            i++;
        }
        return (T)0.5*(a+b);
    }

    template<typename T>
    constexpr T linearinter(T p0, T p1, T t) {
        return (1 - t) * p0 + t * p1;
    }

    // casteljau subdivision using O(N) memory
    template<std::size_t N, typename T>
    void casteljau_subdiv(const std::array<T, N> &coeffs, T t, std::array<T, N> &res_first, std::array<T, N> &res_second) {
        std::array<std::array<T, N>, 2> table{};
        std::size_t curr = 0, prev = 0;
        for (std::size_t i = 0; i < N; ++i) {
            table[curr][i] = coeffs[i];
        }
        curr = 1;
        for (std::size_t i = 1; i < N; ++i) {
            for (std::size_t j = 0; j < N - i; ++j) {
                table[curr][j] = linearinter(table[prev][j], table[prev][j + 1], t);
            }
            res_first[i] = table[curr][0];
            res_second[N - i - 1] = table[curr][N - i - 1];
            prev = curr;
            curr = (1 - curr);
        }
        res_first[0] = coeffs[0];
        res_second[N - 1] = coeffs[N - 1];
    }

    // converts a polynomial from the monomial basis given by coeffs to the bezier basis by
    // storing the result in out. O(N^2) operation with O(N) memory
    template<std::size_t N, typename T>
    constexpr void basis_conversion(const std::array<T, N> &coeffs, std::array<T, N> &out) {
        long bin_c = 1;
        std::array<std::array<T, N>, 2> table{};
        for (std::size_t i = 0; i < N; ++i) {
            table[0][i] = coeffs[i] / ((T) bin_c);
            //careful about operation order
            bin_c = (bin_c * ((N - 1) - i)) / (i + 1);
        }

        std::size_t curr = 0, prev = 1;
        for (std::size_t i = 0; i < N; ++i) {
            out[i] = table[curr][0];
            prev = curr;
            curr = 1-curr;
            for (std::size_t j = 0; j < N - i -1; ++j) {
                table[curr][j] = table[prev][j] + table[prev][j +1];
            }
        }
    }

    // returns the number of sign changes in the coefficients. Numerically zero
    // coefficients  as given by the functor are ignored.
    template<std::size_t N, typename T, typename zeroF>
    constexpr int sign_changes(const std::array<T, N> &coeffs, zeroF& is_zero) {
        std::size_t start = 0, end = N - 1;
        while (start < N - 1 && is_zero(coeffs[start])) {
            start++;
        }
        // start ==  N -1 || c >0;
        while (end > 1 && is_zero(coeffs[end])) {
            end--;
        }
        T prev = coeffs[start];
        int out = 0;
        for (std::size_t i = start; i <= end; ++i) {
            if (is_zero(coeffs[i]) ) {
                continue;
            }
            if (coeffs[i] * prev < 0) {
                out++;
            }
            prev = coeffs[i];
        }
        return out;
    }

    // brackets the roots of the polynomial defined by the coeffs array, passes the interval to the callback
    template<std::size_t N, typename T, typename OutF>
    void rootbracket_bezier(const std::array<T, N> &coeffs, T a, T b, OutF &process, T abstol) {
        using std::abs;

        if (abs(b-a) <= abstol){
            return;
        }
        auto iszero = [&abstol](T d){
            return abs(d) <= abstol;
        };
        bool boundary1 = iszero(coeffs.front());
        bool boundary2 = iszero(coeffs.back());

        switch (sign_changes(coeffs, iszero)) {
            case 0:
                if(boundary1 && boundary2) break;
                if( boundary2 || (a ==0 && boundary1) ){
                    process({a, b});
                }
                return;
            case 1:
                if (!boundary1 && !boundary2) {
                    process({a, b});
                    return;
                }
        }
        T mid = 0.5 * (a + b);
        std::array<T, N> leftcoeffs{};
        std::array<T, N> rightcoeffs{};
        casteljau_subdiv(coeffs, (T) 0.5, leftcoeffs, rightcoeffs);
        rootbracket_bezier(leftcoeffs, a, mid, process, abstol);
        rootbracket_bezier(rightcoeffs, mid, b, process, abstol);
    }

    // Horner polynomial evaluation scheme.
    template<std::size_t N, typename T>
    constexpr T polyval(const std::array<T, N> &coeffs, T t) {
        T out = coeffs[N - 1];
        for (int i = N - 2; i >= 0; i--) {
            out = coeffs[i] + t * out;
        }
        return out;
    }

    // derivative of polynomial
    template<std::size_t N, typename T>
    constexpr void poly_der(const std::array<T, N> &coeffs, std::array<T, N - 1> &out) {
        for (int i = 1; i < N; ++i) {
            out[i - 1] = coeffs[i] * i;
        }
    }

    template<std::size_t N, typename T, bool BezierRepr = false>
    struct PolynomialFunc {
        std::array<T, N> coeffs;
        explicit PolynomialFunc(const std::array<T, N> &poly) : coeffs(poly) {}

        PolynomialFunc<N - 1, T, BezierRepr> derivative(){
            std::array<T, N> derivative_coeffs;
            if constexpr(BezierRepr){
                poly_der_bezier(coeffs, derivative_coeffs);
            } else{
                poly_der(coeffs, derivative_coeffs);
            }
            return {derivative_coeffs};
        }
        T operator()(T x) const{
            if constexpr(BezierRepr){
                return polyeval_bezier(coeffs, x);
            }else{
                return polyval(coeffs, x);
            }
        }
    };

    template<std::size_t N, typename T, typename IntervalConsumer>
    void rootbracket(const std::array<T, N> &poly, IntervalConsumer &out, T abstol = 1e-15) {
        using std::abs;
        T acc = 0;
        for (const auto &c: poly) {
            acc += abs(c);
        }
        if(abs(acc)  <= abstol){
            return;
        }
        std::array<T, N> bezier_coeffs{};
        basis_conversion(poly, bezier_coeffs);
        rootbracket_bezier(bezier_coeffs, (T) 0, (T) 1, out, abstol);
    }

    // numerically computes roots of the input coeffs, given in monomial basis form.
    template<std::size_t N, typename T, typename RootConsumer>
    void getpolyroots(const std::array<T, N> &poly, RootConsumer &out,
                      int niter = 25, T abstol = 1e-15, T interval_eps = 1e-20) {
        //Note that intervals are added in increasing order
        std::size_t numroots =0;
        std::array<std::pair<T, T>, N> root_intervals{};
        auto add = [&](const std::pair<T, T>& interval){
            if(numroots < N) root_intervals[numroots++] = interval;
        };
        using std::abs;
        rootbracket(poly,  add, abstol);
        PolynomialFunc poly_function{poly};
        for (std::size_t i = 0; i < numroots; ++i) {
            auto [a, b] = root_intervals[i];
            // by assumption each interval contains exactly 1 root
            if(abs(poly_function(a))<abstol) {
                out(a);
                continue;
            }
            if(abs(poly_function(b))<abstol) {
                out(b);
                continue;
            };
            if (auto curr = itp_root_refine(poly_function, a, b, (T) interval_eps, niter)){
                out(*curr);
            }
        }
    }

    // convenience function that returns the monomial representation of the polynomial
    // defined by the N input roots. The leading coefficient is not normalized.
    template<std::size_t N, typename T>
    std::array<T, N + 1> poly_from_roots(const std::array<T, N> roots) {
        std::array<T, N + 1> out{};
        out[0] = -roots[0];
        out[1] = 1;
        for (std::size_t i = 1; i < roots.size(); ++i) {
            for (int j = i; j >= 0; --j) {
                out[j + 1] = out[j];
            }
            out[0] = 0;
            for (int j = 0; j < i + 1; ++j) {
                out[j] += out[j + 1] * -roots[i];
            }
        }
        return out;
    }
}
#endif //CONTOURKLIP_POLYNOMIAL_SOLVER_HPP