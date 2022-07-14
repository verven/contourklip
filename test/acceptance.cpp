#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "polyclip.hpp"
#include "test_utilities.hpp"

TEST_CASE("tiny contour fuzzy only lines"){
    for (int seed = 1; seed < 10000; ++seed) {
        contourklip::Contour c1, c2;
        geometrygen::generate_contour(10, 2., 8, 1.5, seed, c1, {}, false);
        geometrygen::generate_contour(10, 2., 8, 1.5, seed + 1, c2, {0.25, 1}, false);
        test_ops(std::vector{c1}, std::vector{c2}, true);
    }
}

TEST_CASE("tiny contour fuzzy tests"){
    int k = 100;
    for (int seed = 1; seed < k; ++seed) {
        for (int seed2 = seed+1; seed2 < k; ++seed2) {
            contourklip::Contour c1, c2;
            geometrygen::generate_contour(5, 2., 8, 1.5, seed, c1, {}, true);
            geometrygen::generate_contour(5, 2., 8, 1.5, seed2, c2, {0.25, 1}, true);
            test_ops(std::vector{c1}, std::vector{c2}, true);
        }
    }
}

TEST_CASE("disallowed contour"){
    contourklip::Contour c1{};
    c1.push_back({0, 1});
    c1.push_back({0, 2});
    c1.push_back({-1, 2});
    c1.push_back({2, 2});
    c1.close();

    contourklip::Contour c2{};
    geometrygen::generate_rectangle({-5, 5}, {1, -2}, c2);
    std::vector<contourklip::Contour> out{};
    bool res = contourklip::clip(c2, c1, out, contourklip::DIFFERENCE);

    CHECK( ! res);
}

TEST_CASE("disallowed contour 2"){
    contourklip::Contour c1{};
    c1.push_back({0, 0});
    c1.push_back({0, 1});
    c1.push_back({2, 2}, {0, 2}, {2, 1});
    c1.push_back({2, 0});
    c1.close();
    std::vector<contourklip::Contour> out{};
    bool res = contourklip::clip(c1, {}, out, contourklip::DIFFERENCE);

    CHECK( ! res);
}