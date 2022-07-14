#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <utility>

#include "doctest.h"

#include "polyclip.hpp"
#include "svg_io.hpp"
#include "test_utilities.hpp"
#include "bezier_utils.hpp"

using namespace contourklip;

TEST_CASE("contour segments"){
    contourklip::Contour c1{{0, 0}};
    c1.push_back({1, 0});
    c1.push_back({2, 0}, {2, 1}, {1, 1});
    c1.push_back({0, 1});
    c1.close();
    int num_lines=0, num_curves= 0;

    auto lines = [&](const Point2d& a, const Point2d& b){
        CHECK(a != b);
        num_lines++;
    };
    auto curves = [&](const Point2d& p0, const Point2d& p1, const Point2d& p2, const Point2d& p3){
        CHECK(p0 != p1);
        num_curves++;
    };
    c1.forward_segments<contourklip::LINE>(lines);
    c1.forward_segments<contourklip::CUBIC_BEZIER>(curves);
    CHECK(num_lines == 3);
    CHECK(num_curves == 1);
}

TEST_SUITE("linear intersections"){

bool check_linear_inter(const Segment &a, const Segment &b) {

    if(!intersect_segments(a, b)){
        return false;
    }

    Segment a1 = {a.second, a.first};
    Segment b1 = {b.second, b.first};

    std::vector<std::pair<Segment, Segment>> cases{
            {a, b},
            {b, a},
            {a1, b1},
            {b1, a1},
            {a1, b},
            {b, a1},
            {a, b1},
            {b1, a}
    };

    std::vector<std::optional<SegInter>> result{};
    for (const auto &seg_pair: cases) {
        result.push_back(intersect_segments(seg_pair.first, seg_pair.second));
        REQUIRE(result.back());
    }

    Point2d p = result.front()->p;
    for (std::size_t i = 0; i < result.size(); ++i) {
        CHECK(result[i]->p == p);
        auto p1 = linear_map(cases[i].first, result[i]->t1);
        auto p2 = linear_map(cases[i].second, result[i]->t2);
        CHECK(approx_equal(p1, p, 1e-5));
        CHECK(approx_equal(p2, p, 1e-5));
    }
    return true;
}

TEST_CASE("linear intersection commutativity")
{
    Segment a{{54.4, 288.5},
              {500, 509.5}};
    Segment b{{500, 288.5},
              {104.4, 500}};

    CHECK(check_linear_inter(a, b));
}

TEST_CASE("linear intersection commutativity 2")
{

    Segment a{{1, 1}, {4, 1}};
    Segment b{{2, 1}, {3, 4}};

    CHECK(check_linear_inter(a, b));
}

TEST_CASE("linear intersection on segment"){

    Segment a{{1, 1}, {5, 1}};
    Segment b{{0, 5}, {3, 1}};
    Segment c{{8, 2}, {3, 1}};

    SegInter u1 = *intersect_segments(a, b);
    SegInter u2 = *intersect_segments(a, c);

    CHECK(u1.p == u2.p);
    }
}

TEST_SUITE("splitting of beziers") {
    TEST_CASE("test casteljau subcurve split 2") {
        CubicBezier c{
                {2089.13, 1077.93},
                {1720.02, 1453.7},
                {1720.02, 1453.7},
                {3411.43, 855.997}
        };
        CubicBezier subcurve = *sub_bezier(c, 0.1, 0.25);
        CHECK(beziermap(c, 0.1) == subcurve.p0);
        CHECK(beziermap(c, 0.25) ==  subcurve.p3);
    }

    TEST_CASE("test bezier monotonic 4 points") {

        CubicBezier c{{5,  -8},
                      {1,  -15},
                      {12, -5},
                      {11, -8}};
        int num = check_monotonic_split(c);
        CHECK(num == 4);
    }

    TEST_CASE("test bezier monotonic split 3 points") {

        CubicBezier c{{5,  -8},
                      {1,  -12},
                      {12, -14},
                      {11, -5}};
        int num = check_monotonic_split(c);
        CHECK(num == 3);
    }

    TEST_CASE("test monotonic bezier already monotonic") {

        CubicBezier c{{1845.298655314596090, 1408.195851834146197},
                      {1769.666501760373194, 1436.156138345120098},
                      {1720.019999999999982, 1453.700000000000045},
                      {1720.019999999999982, 1453.700000000000045}};

        int num = check_monotonic_split(c);
        CHECK(num ==0);
    }

    TEST_CASE("test monotonic bezier with overlapping points") {

        CubicBezier c{
                {1037.349999999999909, 169.037000000000006},
                {1037.349999999999909, 169.037000000000006},
                {865.234133861356099,  313.209452692659681},
                {967.979899428897852,  412.444395647281397}};

        int num = check_monotonic_split(c);
        CHECK(num ==1);
    }

    TEST_CASE("bezier monotonic split 1 point"){

        CubicBezier c{{2517.160194232299091, 1028.619576515788594},
                      {2485.342247547211628, 1028.020143901298070},
                      {2443.735516136445767, 1030.236455866035385},
                      {2392.340000000000600, 1035.268512410000312}};

        int num = check_monotonic_split(c);
        CHECK(num == 1);
    }

    TEST_CASE("bezier monotonic split close extrema"){
        CubicBezier c {{2.52351, -2.82288}, {2.52393, -1.92694},
                       {2.44519, -0.6975}, {2.30787, 0.791758}};

        int num = check_monotonic_split(c);
        CHECK(num == 1);
    }
}

TEST_SUITE("line-curve intersection") {
    TEST_CASE("line curve intersections basic") {
        Segment seg{{1, 2},
                    {5, 3}};
        CubicBezier c{{1, 5},
                      {3, -2},
                      {4, -5},
                      {5, 8}};

        int num = check_line_curve_inter(seg, c);
        CHECK(num == 2);
    }

    TEST_CASE("line curve intersections aligned") {
        Segment seg{{1, 0},
                    {1, 5}};
        CubicBezier c{{1,  1},
                      {-3, 2},
                      {5,  3},
                      {1,  4}};

        int num = check_line_curve_inter(seg, c);
        CHECK(num == 3);
    }

    TEST_CASE("line curve intersections connected2") {
        Segment seg{{1, 0},
                    {9, 4}};
        CubicBezier c{{5.1, 2.05},
                      {-3,  2},
                      {12,  -3},
                      {10,  -5}};

        int num = check_line_curve_inter(seg, c);
        CHECK(num == 2);
    }

    TEST_CASE("line curve intersections vertical connected") {
        Segment seg{{1, 0},
                    {1, 8}};
        CubicBezier c{{1,   0},
                      {-4,  3},
                      {1.7, 4},
                      {1,   8}};
        int num = check_line_curve_inter(seg, c);
        CHECK(num == 1);
    }

    TEST_CASE("line curve intersections vertical") {
        Segment seg{{2392.34, 917.237},
                    {2392.34, 1077.93}};
        CubicBezier c{{2089.13, 1077.93},
                      {3411.43, 855.997},
                      {1720.02, 1453.7},
                      {1720.02, 1453.7}};
        int num = check_line_curve_inter(seg, c);
        CHECK(num ==1);
    }

    TEST_CASE("line curve intersections vertical 2"){

        CubicBezier c{{410.707,862.071}, {330.302,862.071},{265.023,796.792}, {265.023,716.388}};
        Segment seg{{265.023,425.02}, {265.023,716.388}};

        int num = check_line_curve_inter(seg, c);
        CHECK(num == 0);

    }

    TEST_CASE("line curve intersection adjacent") {

        std::cout.precision(25);
        std::cout << std::fixed;

        Segment seg1{{487.548, 643.546},
                     {487.548, 352.179}};
        Segment seg2{{196.181, 643.546},
                     {487.548, 643.546}};
        Segment seg3{{487.548,  352.179},
                     {196.1810, 352.179}};
        CubicBezier c1{{633.231, 497.863},
                       {633.231, 578.268},
                       {567.953, 643.546},
                       {487.548, 643.546}};
        CubicBezier c2{{487.548, 643.546},
                       {407.143, 643.546},
                       {341.864, 578.268},
                       {341.864, 497.863}};
        CubicBezier c3{{341.863999999999976, 497.863000000000000},
                       {341.863999999999976, 417.458000000000027},
                       {407.142999999999972, 352.179},
                       {487.548, 352.179}};
        int num;
        num = check_line_curve_inter(seg1, c1);
        CHECK(num == 0);
        num = check_line_curve_inter(seg1, c2);
        CHECK(num == 0);
        num = check_line_curve_inter(seg2, c1);
        CHECK(num == 0);
        num = check_line_curve_inter(seg3, c2);
        CHECK(num == 0);
        num = check_line_curve_inter(seg3, c1);
        CHECK(num == 0);
        num = check_line_curve_inter(seg3, c3);
        CHECK(num == 0);
        num = check_line_curve_inter(seg2, c3);
        CHECK(num == 0);
    }

    TEST_CASE("line curve intersection almost on line"){

        CubicBezier c{{2.06669, 0}, {7.59173, 0.158415}, {7.20307, 3.47477}, {-1.03651, 3.2766}};
        Segment seg {{-4.73344, 0.0258641}, {3.36403, 0}};

        int num = check_line_curve_inter(seg, c);
        CHECK( num ==  1 );
    }

    TEST_CASE("line curve intersection fuzzy"){
        geometrygen::PointGenerator<Point2d> gen{0, 10, 0, 10};

        auto gen_segment= [&gen](){
            return Segment{gen(), gen()};
        };
        auto gen_bezier= [&gen](){
            return CubicBezier{gen(), gen(), gen(), gen()};
        };
        int num = 0;
        for (int i = 0; i < 10'000; ++i) {
            #ifdef DEBUG
            std::cout << "processing "<< i << "\n";
            #endif
            num += check_line_curve_inter(gen_segment(), gen_bezier());
        }
    }
}

TEST_SUITE("curve curve intersection") {
    TEST_CASE("curve curve intersection already axis aligned") {

        CubicBezier b{{-2, 5},
                      {1,  -2},
                      {3,  -1},
                      {5,  4}};
        CubicBezier q{{0,   0},
                      {0.2, 1.2},
                      {0.5, 1.9},
                      {1.0, 0}};

        int num = check_curve_curve_inter(b, q);
        CHECK(num == 2);
    }

    TEST_CASE("test curve curve inter basic") {
        CubicBezier B1{{1, 2},
                       {3, -4},
                       {4, 10},
                       {5, 8}};//0
        CubicBezier B2{{1,  2},
                       {3,  -14},
                       {4,  20},
                       {10, -10}}; //1
        CubicBezier Q{{-2, 4},
                      {1,  -8},
                      {8,  -5},
                      {10, 9}};
        int num = check_curve_curve_inter(B1, Q);
        CHECK(num == 0);
        num = check_curve_curve_inter(B2, Q);
        CHECK(num == 1);
    }

    TEST_CASE("test curve curve inter circle approximations") {

        CubicBezier B{{562.098, 517.77},
                      {546.036, 628.361},
                      {443.209, 705.107},
                      {332.618, 689.045}};
        CubicBezier Q{{636.808, 555.905},
                      {525.489, 565.729},
                      {427.136, 483.328},
                      {417.312, 372.009}};
        int num = check_curve_curve_inter(B, Q);
        CHECK(num == 1);
    }

    TEST_CASE("bezier intersection start-end") {
        CubicBezier b1{{1196.340, 617.296},
                       {2440.450, 617.296},
                       {1196.340, 1748.460},
                       {2361.720, 1452.770}};
        CubicBezier b2{{1196.340, 617.296},
                       {3321.060, 617.296},
                       {236.996,  1452.770},
                       {2361.720, 1452.770}};
        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 2);
    }

    TEST_CASE("bezier intersection start-end self intersections") {
        Point2d p1 = {528.31, 12.1};
        Point2d p2 = {801, 49.5};

        CubicBezier b1{p1,
                       {2440.450, 617.296},
                       {1196.340, 1748.460},
                       p2};
        CubicBezier b2{p1,
                       {3321.060, 617.296},
                       {236.996,  1452.770},
                       p2};

        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 4);
    }

    TEST_CASE("bezier intersection start-end 4") {

        Point2d p1 = {528.31, 12.1};
        Point2d p2 = {801, 49.5};

        CubicBezier b1{p1,
                       {1196.340, 1748.460},
                       {2440.450, 617.296},
                       p2};
        CubicBezier b2{p1,
                       {236.996,  1452.770},
                       {3321.060, 617.296},
                       p2};

        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 2);
    }

    TEST_CASE("bezier bezier intersections circle_connected") {

        CubicBezier c1{{633.231, 497.863},
                       {633.231, 578.268},
                       {567.953, 643.546},
                       {487.548, 643.546}};
        CubicBezier c2{{487.548, 643.546},
                       {407.143, 643.546},
                       {341.864, 578.268},
                       {341.864, 497.863}};

        int num = check_curve_curve_inter(c1, c2);
        CHECK(num == 0);
    }

    TEST_CASE("bezier bezier intersections quarter circle") {
        CubicBezier c1{{50, 100},{77.5, 100}, {100, 77.5}, {100, 50}};
        CubicBezier c2{{50, 75},{50, 102.5}, {72.3, 125}, {100, 125}};
        int num = check_curve_curve_inter(c1, c2);
        CHECK(num == 1);
    }

    TEST_CASE("bezier intersection start-start horizontal"){
        CubicBezier b1{{500, 500}, {525, 550}, {530, 480}, {550, 500}};
        CubicBezier b2{{500, 500}, {530, 490}, {535, 525}, {540, 550}};
        int num = check_curve_curve_inter(b1, b2);
        CHECK(num ==1);
    }

    TEST_CASE("curve-curve 5 intersections"){
        CubicBezier b1{{308.1, 725.4}, {129.4, 93.6}, {792.2, 941.1}, {594.6, 433.6}};
        CubicBezier b2{{304.2, 734.6}, {168.1, 70.2}, {687, 968.8}, {603.8, 409.5}};
        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 5);
    }

    TEST_CASE("curve-curve 9 intersections"){
        CubicBezier b1{{385.1, 514.3}, {802.1, 514.3}, {99.8, 412.8}, {547.4, 412.8}};
        CubicBezier b2{{394.9, 528.3}, {394.9, 169}, {522.7, 739.3}, {522.7, 398.2}};
        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 9);
    }

    TEST_CASE("curve-curve 3 intersections"){
        CubicBezier b1{{385.1, 514.3}, {802.1, 514.3}, {99.8, 412.8}, {547.4, 412.8}};
        CubicBezier b2{{385.1, 512.3}, {802.1, 512.3}, {99.8, 414.8}, {547.4, 414.8}};

        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 3);
    }

    TEST_CASE("curve-curve 1 intersection"){
        CubicBezier b1{{464.7, 372.8}, {452.9, 438.5}, {171.1, 533.4}, {121.9, 488.2}};
        CubicBezier b2{{500, 580.6}, {467.8, 580.6},{356.7, 469.5}, {356.7, 437.3} };
        int num = check_curve_curve_inter(b1, b2);
        CHECK(num == 1);
    }
}

TEST_SUITE("segment ordering") {
    TEST_CASE("bezier ordering visual check") {

        std::string dir = TESTCASE_OTHERS_DIR + "/bezier_ordering/";
        BasicPathLoader pl(dir + "in.svg");

        auto beziercurve_below = [](const CubicBezier &a, const CubicBezier &b) -> bool {

            detail::SweepPoint a1{a.p0};
            detail::SweepPoint a2{a.p3};
            detail::SweepPoint b1{b.p0};
            detail::SweepPoint b2{b.p3};
            a1.controlp = a.p1;
            a2.controlp = a.p2;
            b1.controlp = b.p1;
            b2.controlp = b.p2;
            a1.other_point = &a2;
            a2.other_point = &a1;
            b1.other_point = &b2;
            b2.other_point = &b1;

            a1.set_if_left();
            b1.set_if_left();

            return curve_below(&a1, &b1);
        };

        auto output_order = [&](bool  reflected ) {
            BasicPathWriter out(pl.get_dims());
            std::vector<CubicBezier> v{};
            for (const auto &path: pl.paths) {
                v.push_back(bezier_from_str(path));
                double offset = +v.back().p0.x();
                if(reflected) {
                    auto mirrorf = [&](const Point2d &p) -> Point2d { return {-p.x() + 5 * offset, p.y()}; };
                    v.back() = transform(v.back(), mirrorf);
                }
                out.push_path_str(bezier_to_path_str(v.back()), "none");
            }

            std::sort(v.begin(), v.end(), beziercurve_below);

            int idx = 0;
            for (const auto &curve: v) {
                out.push_text(beziermap(curve, 0.55), std::to_string(idx));
                idx++;
            }
            out.write_to(dir + (reflected ? "out_reflected.svg" : "out.svg") );
        };
        output_order(false);
        output_order(true);
    }

    TEST_CASE("line-bezier ordering") {
        Segment seg{{1166.7, 350},
                    {1352.3, 505.9}};

        CubicBezier c{{1166.7, 350},
                      {1166.7, 350},
                      {1176.7, 148.3},
                      {1232.1, 148.3}};

        detail::SweepPoint s1;
        detail::SweepPoint s2;
        detail::SweepPoint c1;
        detail::SweepPoint c2;
        sweeppoint_init(seg, &s1, &s2);
        sweeppoint_init(c, &c1, &c2);

        CHECK(curve_below(&c1, &s1));
        CHECK(detail::queue_comp<LeftOfLine, IsCollinear>(&c1, &s1));
    }
    TEST_CASE("line-bezier ordering horizontal-vertical") {

        Segment seg{{269.023, 425.02}, {560.39, 425.021}};

        CubicBezier c{{269.023, 425.02}, {269.023, 401.929},
                      {274.408, 380.083},{283.989, 360.674}};
        detail::SweepPoint s1;
        detail::SweepPoint s2;
        detail::SweepPoint c1;
        detail::SweepPoint c2;
        sweeppoint_init(seg, &s1, &s2);
        sweeppoint_init(c, &c1, &c2);

        CHECK(curve_below(&c1, &s1));
        CHECK(detail::queue_comp(&c1, &s1));
    }

    TEST_CASE("line-bezier ordering 2") {
        CubicBezier c{{2.523510524871271, -2.822884461255850}, {2.523511628233836, -2.820518727440074},
                {2.523511628233836, -2.820518727440074}, {2.523512179725838, -2.815780285295673}};

        Segment seg{{2.523510524871271, -2.822884461255850}, {6.795631221097975, 0.000000000000000}};
        detail::SweepPoint s1;
        detail::SweepPoint s2;
        detail::SweepPoint c1;
        detail::SweepPoint c2;
        sweeppoint_init(seg, &s1, &s2);
        sweeppoint_init(c, &c1, &c2);
        CHECK( curve_below(&s1, &c1));
    }

    TEST_CASE("line-bezier ordering overlapping"){

        Segment seg{{15.068416403647555, 12.765976081320677}, {15.070287752788607, 12.757544852322230}};
        CubicBezier c{seg.first, {15.069040488591890, 12.763164135717986},
                      {15.069664271610140, 12.760353726150072}, seg.second};

        detail::SweepPoint s1;
        detail::SweepPoint s2;
        detail::SweepPoint c1;
        detail::SweepPoint c2;
        sweeppoint_init(seg, &s1, &s2);
        sweeppoint_init(c, &c1, &c2);

        s1.ptype = s2.ptype = c1.ptype = c2.ptype = contourklip::detail::CLIPPING;
        CHECK(detail::curve_below(&s1, &c1) !=  detail::curve_below(&c1, &s1));
    }
}

TEST_SUITE("postprocess contour"){
    TEST_CASE("postprocess contour 1") {

    contourklip::Contour c;
    c.push_back({0, 0});
    c.push_back({0, 1});
    c.push_back({0.5, 0.5});
    c.push_back({1, 0});
    c.push_back({1, 1});
    c.push_back({0.5, 0.5});
    c.close();

    std::vector<contourklip::Contour> out{};
    auto process = [&](const contourklip::Contour &c) {
        out.push_back(c);
    };
    postprocess_contour(c, process);
    REQUIRE(out.size()==2);
    CHECK(out.front().size() == 4);
    CHECK(out.back().size() == 4);
    }

    TEST_CASE("postprocess contour 2") {

        contourklip::Contour c;
        c.push_back({0.5, 0.5});
        c.push_back({0, 0});
        c.push_back({0, 1});
        c.push_back({0.5, 0.5});
        c.push_back({1, 1});
        c.push_back({1, 0});
        c.close();

        std::vector<contourklip::Contour> out{};
        auto process = [&](const contourklip::Contour &c) {
            out.push_back(c);
        };
        postprocess_contour(c, process);
        REQUIRE(out.size()==2);
        CHECK(out.front().size() == 4);
        CHECK(out.back().size() == 4);
    }

    TEST_CASE("postprocess contour 3") {
        contourklip::Contour c;
        c.push_back({0, 2});
        c.push_back({1, 2});
        c.push_back({2, 2});
        c.push_back({3, 2});
        c.push_back({3, 1});
        c.push_back({4, 1});
        c.push_back({4, 3});
        c.push_back({2, 3});
        c.push_back({2, 2});
        c.push_back({2, 1});
        c.push_back({3, 1});
        c.push_back({3, 0});
        c.push_back({0, 0});
        c.close();

        std::vector<contourklip::Contour> out{};
        auto process = [&](const contourklip::Contour &c) {
            out.push_back(c);
        };
        postprocess_contour(c, process);
        REQUIRE(out.size()==2);
        CHECK(out.front().size() == 7);
        CHECK(out.back().size() == 7);
    }
}

TEST_SUITE("geometry generation") {
    TEST_CASE("random partition") {
        double a = 0, b = 2 * M_PI_2;
        for (int i = 2; i < 50; ++i) {
            for (int seed = 0; seed < 10; ++seed) {
                std::vector<double> out = geometrygen::random_partition(a, b, i, 1);
                CHECK(out.size() == i - 1);
                CHECK(std::is_sorted(out.begin(), out.end()));
            }
        }
    }

    TEST_CASE("point generation"){
        geometrygen::PointGenerator<Point2d> getp{0, 1, 0, 1, 5};
        std::set<Point2d> generated{};
        for (int i = 0; i <25; ++i) {
            Point2d p{getp()};
            CHECK(p.x() <= 1);
            CHECK(p.x() >= 0);
            CHECK(p.y() <= 1);
            CHECK(p.y() >= 0);
            CHECK(generated.find(p) == generated.end());
            generated.insert(p);
        }
    }

    TEST_CASE("random contour generation"){
        Contour out;
        geometrygen::generate_contour(30, 2., 5., 1.5, 15, out, {}, true);
        save_contour("debug_random_contour.svg", out);
    }
}

TEST_CASE("specific testcase") {
    std::string case_name = "curves03";
    std::string t_dir = std::string{TESTCASE_DIR} + "/" + case_name + "/";

    BooleanOpType op = contourklip::INTERSECTION;
    std::vector<Contour> a;
    std::vector<Contour> b;
    std::vector<Contour> a_rev;
    std::vector<Contour> b_rev;

    auto dims = load_testcase(case_name, a, b);
    a_rev = a; b_rev = b;
    reverse_contours(a_rev);
    reverse_contours(b_rev);

    std::vector<Contour> res;
    BasicPathWriter dbg{dims};
    dbg.path_prefx = t_dir + "/";

    Config config;
    PolyClip t(a, b, res, op, config);
#ifdef DEBUG
    t.w = &dbg;
    t.compute_verbose = true;
#endif
    t.compute();
    CHECK(t.success());
    if(!res.empty()) CHECK(contour_area(res.front()) > 0);
    BasicPathWriter out{dims};
    if (op == DIVIDE) {
        for (const auto &c: res) {
            //this way, each contour has its own path
            out.push_path_str(multipolygon_to_str(std::vector<Contour>{c}));
        }
    } else {
        out.push_path_str(multipolygon_to_str(res));
    }

    #ifdef DEBUG
    std::cout << "total area of output " << multipolygon_area(res) << "\n";
    std::cout << "consisting of areas:\n";
    for (const auto &c: res) {
        std::cout << contour_area(c) << "\n";
    }
    save_output(case_name, out, op);
    #endif
}