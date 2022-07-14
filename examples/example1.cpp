#include <iostream>
#include "contourklip.hpp"

int main() {
    contourklip::Contour contour1{{0, 100}};
    contour1.push_back({50, 100});
    contour1.push_back({77.5, 100}, {100, 77.5}, {100, 50});
    contour1.push_back({100, 22.5}, {77.5, 0}, {50, 0});
    contour1.push_back({0, 0});
    contour1.close();

    contourklip::Contour contour2{{150, 25}};
    contour2.push_back({100, 25});
    contour2.push_back({72.3, 25}, {50, 47.3}, {50, 75});
    contour2.push_back({50, 102.5}, {72.3, 125}, {100, 125});
    contour2.push_back({150, 125});
    contour2.close();

    std::vector<contourklip::Contour> shape1{contour1};
    std::vector<contourklip::Contour> shape2{contour2};
    std::vector<contourklip::Contour> result{};

    if(contourklip::clip(shape1, shape2, result, contourklip::INTERSECTION)){
        std::cout << "clipping operation succeeded\n";
    }

    for (const auto &contour: result) {
        std::cout << "contour:\n";
        std::cout << contour;
    }
    return 0;
}