#ifndef CONTOURKLIP_TEST_UTILITIES_HPP
#define CONTOURKLIP_TEST_UTILITIES_HPP

#include "svg_io.hpp"
#include "geometry_generation.hpp"
#include "bezier_utils.hpp"
#include "polyclip.hpp"

#include "doctest.h"

std::string TESTCASE_ROOT_DIR = "./../testcases";
std::string TESTCASE_DIR = std::string{TESTCASE_ROOT_DIR} + "/geometry" ;
std::string TESTCASE_OTHERS_DIR = std::string{TESTCASE_ROOT_DIR} + "/various" ;
constexpr double AREA_EPS = 1;

void set_cout_precision(int prec = 10){
    std::cout.setf(std::ios::unitbuf);
    std::cout.precision(prec);
    std::cout << std::fixed;
}

std::tuple<double, double, double, double> load_testcase(const std::string& case_name,
                                                         std::vector<contourklip::Contour> &a,
                                                         std::vector<contourklip::Contour> &b){
    std::string in_path = TESTCASE_DIR +"/" + case_name +"/in.svg";
    BasicPathLoader pl(in_path);
    REQUIRE(pl.paths.size() >=2);
    multipolygon_from_str(pl.paths[0], a);
    multipolygon_from_str(pl.paths[1], b);
    return pl.get_dims();
}

bool parametric_v(double t){
    bool numeric= !std::isnan(t) && !std::isinf(t);
    bool tinrange =0<=t && t<=1;
    return numeric && tinrange;
}

void sweeppoint_init(const contourklip::Segment &seg, contourklip::detail::SweepPoint *a, contourklip::detail::SweepPoint *b) {

    if (contourklip::increasing(seg.first, seg.second)) {
        a->point = seg.first;
        b->point = seg.second;
    } else {
        b->point = seg.first;
        a->point = seg.second;
    }
    a->left = true;
    b->left = false;
    a->other_point = b;
    b->other_point = a;
    a->controlp = a->other_point->point;
    b->controlp = b->other_point->point;
    a->curve = b->curve = false;
}

void sweeppoint_init(const contourklip::CubicBezier &c, contourklip::detail::SweepPoint *a, contourklip::detail::SweepPoint *b) {

    if (contourklip::increasing(c.p0, c.p3)) {
        a->point = c.p0;
        b->point = c.p3;
        a->controlp = c.p1;
        b->controlp = c.p2;
    } else {
        a->point = c.p3;
        b->point = c.p0;
        b->controlp = c.p1;
        a->controlp = c.p2;
    }
    a->left = true;
    b->left = false;
    a->other_point = b;
    b->other_point = a;
    a->curve = b->curve = true;
}

std::tuple<double, double, double, double> load_testcase(const std::string &case_name, contourklip::Contour &a,
                                                         contourklip::Contour &b) {
    std::vector<contourklip::Contour> p1{};
    std::vector<contourklip::Contour> p2{};
    auto dims = load_testcase(case_name, p1, p2);
    REQUIRE(!p1.empty());
    REQUIRE(!p2.empty());
    a = p1.front();
    b = p2.front();
    return dims;
}

std::tuple<double, double, double, double> load_contour(const std::string &in_path, contourklip::Contour &a) {

    std::vector<contourklip::Contour> p1{};
    BasicPathLoader pl(in_path);
    multipolygon_from_str(pl.paths.front(), p1);
    REQUIRE(!p1.empty());
    a = p1.front();

    return pl.get_dims();
}

void save_contour(const std::string &out_path, contourklip::Contour &a){
    auto [min_x, min_y, max_x, max_y] = contourklip::contourbbox(a);
    BasicPathWriter w{min_x, min_y, (max_x-min_x), (max_y-min_y)};
    w.push_path_str(multipolygon_to_str({a}));
    w.write_to(out_path);
}

void save_contours(const std::string &out_path, std::vector<contourklip::Contour> &a){
    if(a.empty()) return;
    if(a.front().size() ==0 ) return;
    auto [min_x, min_y , max_x , max_y] = contourklip::contourbbox(a.front());
    for (const auto &contour: a) {
        auto [cmin_x, cmin_y, cmax_x, cmax_y] = contourklip::contourbbox(contour);
        min_x = std::min(cmin_x, min_x);
        min_y = std::min(cmin_y, min_y);
        max_x = std::max(cmax_x, max_x);
        max_y = std::max(cmax_y, max_y);
    }
    double offset = 5;
    min_x -= offset;
    max_x += offset;
    min_y -= offset;
    max_y += offset;
    BasicPathWriter w{min_x, min_y, (max_x-min_x), (max_y-min_y)};
    w.push_path_str(multipolygon_to_str(a));
    w.write_to(out_path);
}

void save_output(const std::string &case_name, BasicPathWriter &out, const contourklip::BooleanOpType& op) {
    std::ostringstream out_path{};
    out_path << TESTCASE_DIR << "/" << case_name << "/" << op <<  "_out.svg";
//    std::string out_path = std::string{TESTCASE_DIR} + "/" + case_name + "/out.svg";
    out.write_to(out_path.str());
}

void reverse_contours(std::vector<contourklip::Contour>& poly){
    for (auto &contour: poly) {
        contour.reverse();
    }
}

int check_monotonic_split(const contourklip::CubicBezier& c){

    auto horizontalorvertical
    = [](const contourklip::Point2d& a, const contourklip::Point2d& b){
        return (std::abs(a.y()-b.y()) < 1e-8) || (std::abs(a.x()-b.x()) < 1e-8);
    };

    double t =0;
    int num =0;
    auto dosplit = [&](double u, contourklip::Extremity_Direction d){
        REQUIRE(u != t);
        auto maybesubsegment = contourklip::sub_bezier(c, t, u);
        REQUIRE(maybesubsegment);
        contourklip::CubicBezier subsegment = *maybesubsegment;
        bool tmp = horizontalorvertical(subsegment.p0, subsegment.p1)
                   || horizontalorvertical(subsegment.p2, subsegment.p3);
        CHECK(tmp);
        t = u;
        num++;
    };
    contourklip::bezier_monotonic_split(c, dosplit);
    CHECK(num <= 4);
    if(num >0){
        dosplit(1, contourklip::X_Extremity);
        //n inner points lead to n+1 intervals
        num--;
    }
    return num;
}

int check_line_curve_inter(const contourklip::Segment& seg, const contourklip::CubicBezier& c){
    int num =0;
    auto verify = [&seg , &c, &num](double t, double u, contourklip::Point2d p) {
        CHECK(parametric_v(t));
        CHECK(parametric_v(u));
        contourklip::Point2d interp = linear_map(seg, t);
        contourklip::Point2d interp2 = beziermap(c, u);
        CHECK(contourklip::detail::approx_equal(interp, interp2, 1e-4));
#ifdef DEBUG
        std::cout << "intersection at " << t << ", " << u << " with point " << p << '\n';
#endif
        num++;
    };

    contourklip::line_bezier_inter(seg, c, verify);
    return num;
}

int check_curve_curve_inter(const contourklip::CubicBezier& a, const contourklip::CubicBezier& b){

    int num = 0;
    std::vector<contourklip::Point2d> inters{};
    std::vector<contourklip::Point2d> inters1{};
    bool switched = false;

    #ifdef DEBUG
    std::cout.precision(20);
    std::cout << std::fixed;
    #endif
    auto verify = [&](double t, double u, contourklip::Point2d p) {
        CHECK(parametric_v(t));
        CHECK(parametric_v(u));
        if (switched) std::swap(t, u);
        contourklip::Point2d interp = beziermap(a, t);
        contourklip::Point2d interp2 = beziermap(b, u);
        #ifdef DEBUG
        std::cout << "intersection at " << t << " " << u << interp << " " << interp2 << '\n';
        #endif
        CHECK(contourklip::detail::approx_equal(interp, interp2, 1e-4));
        CHECK(contourklip::detail::approx_equal(interp, p, 1e-4));
        inters.push_back(p);
    };

    contourklip::curve_curve_inter(a, b, verify, (double) 1e-9);
#ifdef DEBUG
    std::cout << "------------------------\n";
#endif
    CHECK(inters.size()<=9);
    num = (int)inters.size();
    switched= true;
    inters1 = inters;
    inters.clear();

    contourklip::curve_curve_inter(b, a, verify, (double) 1e-9);

    CHECK(inters.size()<=9);
    CHECK(inters1.size() == inters.size());
    CHECK(inters1 == inters);
  #ifdef DEBUG
    for (int i = 0; i < inters.size(); ++i) {
        std::cout << inters1[i] << " "<< inters[i] << '\n';
    }
  #endif

    return num;
}


bool area_check(const std::vector<contourklip::Contour> &a, const std::vector<contourklip::Contour> &b){
    return (contourklip::detail::multipolygon_area(a) - contourklip::detail:: multipolygon_area(b)) < AREA_EPS;
}

bool area_check_inter(const std::vector<contourklip::Contour> &u_ab, const std::vector<contourklip::Contour> &xor_ab, const std::vector<contourklip::Contour> &i_ab){

    double r1 = contourklip::detail::multipolygon_area(u_ab);
    double r2 = contourklip::detail::multipolygon_area(xor_ab);
    double r3 =  contourklip::detail::multipolygon_area(i_ab);

    return std::abs(r1-(r2+r3)) < AREA_EPS;
}


bool area_check_diff(const std::vector<contourklip::Contour> &a_sub_b, const std::vector<contourklip::Contour> &b_sub_a, const std::vector<contourklip::Contour> &xor_ab){

    double r1 = contourklip::detail::multipolygon_area(xor_ab);
    double r2 = contourklip::detail::multipolygon_area(a_sub_b);
    double r3 = contourklip::detail::multipolygon_area(b_sub_a);
    #ifdef DEBUG
    std::cout << "xor, a-b, b-a areas: " << r1 << ", " << r2 << ", " << r3 << '\n';
    #endif
    return std::abs(r1 - (r2+r3)) < AREA_EPS;
}

void test_op_unary(const std::vector<contourklip::Contour> &a){
    std::vector<contourklip::Contour> a_inter_a{};
    std::vector<contourklip::Contour> a_union_a{};
    std::vector<contourklip::Contour> a_diff_a{};
    std::vector<contourklip::Contour> a_xor_a{};

    std::vector<contourklip::Contour> a_union_none{};
    std::vector<contourklip::Contour> none_union_a{};
    std::vector<contourklip::Contour> a_inter_none{};
    std::vector<contourklip::Contour> none_inter_a{};

   #ifdef DEBUG
    std::cout << "testing a inter a\n";
   #endif
    CHECK(contourklip::clip(a, a, a_inter_a, contourklip::INTERSECTION));
    #ifdef DEBUG
    std::cout << "testing a union a\n";

    #endif
    CHECK(contourklip::clip(a, a, a_union_a, contourklip::UNION));

    #ifdef DEBUG
    std::cout << "testing a inter none\n";
    #endif
    CHECK(contourklip::clip(a, {}, a_inter_none, contourklip::INTERSECTION));
   #ifdef DEBUG
    std::cout << "testing none union a\n";
   #endif
    CHECK(contourklip::clip({}, a, none_union_a, contourklip::INTERSECTION));
    #ifdef DEBUG
    std::cout << "testing a union none\n";
    #endif
    CHECK(contourklip::clip(a, {}, a_union_none, contourklip::UNION));

}

auto test_op_binary(const std::vector<contourklip::Contour> &a, const std::vector<contourklip::Contour> &b,
                    bool check_area = true){
    std::vector<contourklip::Contour> a_inter_b{};
    std::vector<contourklip::Contour> a_union_b{};
    std::vector<contourklip::Contour> b_inter_a{};
    std::vector<contourklip::Contour> b_union_a{};
    std::vector<contourklip::Contour> a_sub_b{};
    std::vector<contourklip::Contour> b_sub_a{};
    std::vector<contourklip::Contour> a_xor_b{};
    std::vector<contourklip::Contour> b_xor_a{};

    std::vector<contourklip::Contour> a_dissolve_b{};
    std::vector<contourklip::Contour> b_dissolve_a{};

    #ifdef DEBUG
    std::cout << "testing intersection a b\n";
    #endif
    CHECK(contourklip::clip(a, b, a_inter_b, contourklip::INTERSECTION));
    #ifdef DEBUG
    std::cout << "testing intersection b a\n";
    #endif
    CHECK(contourklip::clip(b, a, b_inter_a, contourklip::INTERSECTION));
    #ifdef DEBUG
    std::cout << "testing union a b\n";
    #endif
    CHECK(contourklip::clip(a, b, a_union_b, contourklip::UNION));
    #ifdef DEBUG
    std::cout << "testing union b a\n";
    #endif
    CHECK(contourklip::clip(b, a, b_union_a, contourklip::UNION));
    #ifdef DEBUG
    std::cout << "testing difference a b\n";
    #endif
    CHECK(contourklip::clip(a, b, a_sub_b, contourklip::DIFFERENCE));
    #ifdef DEBUG
    std::cout << "testing difference b a\n";
    #endif
    CHECK(contourklip::clip(b, a, b_sub_a, contourklip::DIFFERENCE));
    #ifdef DEBUG
    std::cout << "testing xor a b\n";
    #endif
    CHECK(contourklip::clip(a, b, a_xor_b, contourklip::XOR));
    #ifdef DEBUG
    std::cout << "testing xor b a\n";
    #endif
    CHECK(contourklip::clip(b, a, b_xor_a, contourklip::XOR));

    #ifdef DEBUG
    std::cout << "testing dissolve a, b\n";
    #endif
    CHECK(contourklip::clip(a, b, a_dissolve_b, contourklip::DIVIDE));
    #ifdef DEBUG
    std::cout << "testing dissolve b, a\n";
    #endif
    CHECK(contourklip::clip(b, a, b_dissolve_a, contourklip::DIVIDE));

    if(check_area) {
        CHECK(area_check(a_inter_b, b_inter_a));
        CHECK(area_check(a_union_b, b_union_a));
        CHECK(area_check(a_xor_b, b_xor_a));
        CHECK(area_check(a_dissolve_b, b_dissolve_a));

        CHECK(area_check_inter(a_union_b, a_xor_b, a_inter_b));
        CHECK(area_check_diff(a_sub_b, b_sub_a, a_xor_b));
    }

    return std::array<double, 4>{contourklip::detail::multipolygon_area(a_inter_b),
                                 contourklip::detail::multipolygon_area(a_union_b),
                                 contourklip::detail::multipolygon_area(a_sub_b),
                                 contourklip::detail::multipolygon_area(b_sub_a)};
}

void test_ops(const std::vector<contourklip::Contour> &a, const std::vector<contourklip::Contour> &b,
              bool check_area = true){
    std::vector<contourklip::Contour> a_rev = a;
    std::vector<contourklip::Contour> b_rev = b;
    reverse_contours(a_rev);
    reverse_contours(b_rev);

    test_op_unary(a);
    test_op_unary(a_rev);
    test_op_unary(b);
    test_op_unary(b_rev);

    test_op_binary(a, b, check_area);
    test_op_binary(a, b_rev, check_area);
    test_op_binary(a_rev, b, check_area);
    test_op_binary(a_rev, b_rev, check_area);
    test_op_binary(a, a_rev, check_area);
    test_op_binary(b, b_rev, check_area);
}

auto do_test = [](const std::string& tcase){
    std::vector<contourklip::Contour> a;
    std::vector<contourklip::Contour> b;
    load_testcase(tcase, a, b);
    test_ops(a, b);
};
#endif //CONTOURKLIP_TEST_UTILITIES_HPP