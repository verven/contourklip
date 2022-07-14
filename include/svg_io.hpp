#ifndef CONTOURKLIP_SVG_IO_HPP
#define CONTOURKLIP_SVG_IO_HPP

#include <array>
#include <sstream>
#include <regex>
#include <filesystem>
#include <iostream>
#include <fstream>
#include "geometry_base.hpp"

/*
 * some functionality to instantiate contours from a svg file, and write back contours
 * to a svg file. It works for svg files which have the same structure as the
 * testcases, and it is therefore quite limited and not very robust.
 */


void path_to_str(const contourklip::Contour &c, std::string &out, bool close_path = true) {
    if (c.size() == 0) return;

    std::ostringstream os;
    os << std::setprecision(10);
    os << "M";

    auto svgf = [&os](const contourklip::Point2d &p) {
        os << p.x() << ", " << p.y();
    };

    svgf(c.front_point());

    for (std::size_t i = 1; i < c.size(); ++i) {
        if (c[i].segment_shape() == contourklip::LINE) {
            os << (i == 0 ? " " : " L");
            svgf(c[i].point());
        } else {
            os << (i == 0 ? " " : " C");
            svgf(c[i].c1());
            os << " ";
            svgf(c[i].c2());
            os << " ";
            svgf(c[i].point());

        }
    }
    if (close_path) os << "Z";
    out += os.str();
}

std::string path_to_str(const contourklip::Contour &c, bool close_path = true) {
    std::string out{};
    path_to_str(c, out, close_path);
    return out;
}

void multipolygon_to_str(const std::vector<contourklip::Contour> &p, std::string &out) {
    for (const auto &c: p) {
        path_to_str(c, out);
    }
}

std::string multipolygon_to_str(const std::vector<contourklip::Contour> &p) {
    std::string out;
    multipolygon_to_str(p, out);
    return out;
}

std::string segment_to_path_str(const contourklip::Point2d &p0, const contourklip::Point2d p1) {

    contourklip::Contour p{p0, p1};
    std::string out;
    path_to_str(p, out, false);
    return out;
}

std::string segment_to_path_str(const contourklip::detail::Segment &seg) {
    return segment_to_path_str(seg.first, seg.second);
}

std::string
bezier_to_path_str(const contourklip::Point2d &p0, const contourklip::Point2d p1, const contourklip::Point2d p2,
                   const contourklip::Point2d p3) {
    contourklip::Contour p{p0, p1, p2, p3};
    std::string out;
    path_to_str(p, out, false);
    return out;
}

std::string bezier_to_path_str(const contourklip::detail::CubicBezier &c) {
    return bezier_to_path_str(c.p0, c.p1, c.p2, c.p3);
}

std::string bezier_to_path_str(const std::array<double, 4> &Bx, const std::array<double, 4> &By) {
    return bezier_to_path_str({Bx[0], By[0]}, {Bx[1], By[1]},
                              {Bx[2], By[2]}, {Bx[3], By[3]});
}

void str_replace(std::string &str, const std::string &toreplace, const std::string &with) {
    auto length = toreplace.size();
    auto length_new = with.size();
    std::size_t index = 0;
    while ((index = str.find(toreplace, index)) != std::string::npos) {
        str.replace(index, length, with);
        index += length_new;
    }
}

bool path_from_svg_string(std::istringstream &path_str, contourklip::Contour &out, bool verbose = false) {

    if (path_str.get() != 'M') {
        return 0;
    }

    auto parsepoint = [&path_str]() -> contourklip::Point2d {
        double curr_x, curr_y;
        path_str >> curr_x;
        if (path_str.peek() == ',') {
            char t;
            path_str >> t; // consume the separator
        }
        path_str >> curr_y;

        return contourklip::Point2d{curr_x, curr_y};
    };
    contourklip::Point2d start = parsepoint();
    if (verbose) std::cout << "starting point_: " << start << '\n';

    out.push_back(start);

    while (path_str.peek() != EOF) {
        char curr;
        path_str >> curr;
        if (verbose) std::cout << "curr " << curr << '\n';
        if (path_str.bad()) {
            return false;
        }
        switch (curr) {
            case 'Z':
                if (verbose) std::cout << "parse done\n";
                return true;
            case 'L': {
                contourklip::Point2d p = parsepoint();
                out.push_back(p);
                if (verbose) std::cout << "parsed L" << p << '\n';
            }
                continue;
            case 'C': {
                contourklip::Point2d p1 = parsepoint();
                contourklip::Point2d p2 = parsepoint();
                contourklip::Point2d pLast = parsepoint();
                out.push_back(p1, p2, pLast);
                if (verbose) std::cout << "parsed C " << p1 << " " << p2 << pLast << '\n';
            }
                continue;
            case 'H': {
                double x_curr;
                path_str >> x_curr;
                contourklip::Point2d p1{x_curr, out.back_point().y()};
                if (verbose) std::cout << "parsed H" << p1 << '\n';
                out.push_back(p1);
            }
                continue;

            case 'V': {
                double y_curr;
                path_str >> y_curr;
                contourklip::Point2d p1{out.back_point().x(), y_curr};
                if (verbose) std::cout << "parsed V" << p1 << '\n';
                out.push_back(p1);
            }
                continue;
            default:
                if (verbose) std::cout << "bad character encountered";
                return false;
        }
    }
    return true;
}

void multipolygon_from_str(const std::string &in, std::vector<contourklip::Contour> &out, bool verbose = false) {

    if (in.size() <= 1) {
        return;
    }
    std::string copy = in;
    std::transform(copy.begin(), copy.end(), copy.begin(), ::toupper);
    //we do this so that parsing with istringstream is much easier. Of course, this is not efficient.
    str_replace(copy, "C", " C");
    str_replace(copy, "L", " L");
    str_replace(copy, "Z", " Z");
    str_replace(copy, "H", " H");
    str_replace(copy, "V", " V");

    if (verbose) std::cout << "complete: \n" << copy << "\n\n";

    std::istringstream path_str(copy);

    while (path_str.peek() != EOF) {
        contourklip::Contour curr{};
        if (path_from_svg_string(path_str, curr)) {
            out.push_back(curr);
        }
    }
}

contourklip::detail::CubicBezier bezier_from_str(const std::string &path) {
    std::vector<contourklip::Contour> poly{};
    multipolygon_from_str(path, poly);
    assert(!poly.empty() && poly.at(0).size() > 1);
    auto seg = poly.at(0)[1];
    assert(seg.segment_shape() == contourklip::CUBIC_BEZIER);
    return {poly.at(0).front_point(), seg.c1(), seg.c2(), seg.point()};
}

class BasicPathLoader {
    std::string svg_in;
public:
    explicit BasicPathLoader(const std::string &filepath) {
        std::ifstream ifs{filepath};
        if (!ifs) {
            std::cout << "file opening unsuccessful. requested path:\n    ";
            std::cout << filepath;
            std::cout << "\ncurrent path is\n    ";
            std::cout << std::filesystem::current_path().string() << std::endl << std::flush;
            std::cout << '\n';
            assert(false);
        }
        svg_in.assign((std::istreambuf_iterator<char>(ifs)),
                      (std::istreambuf_iterator<char>()));

        preprocess();
        parse_dims();
        parse_paths();
    }

    std::tuple<double, double, double, double> get_dims() {
        return {start_x, start_y, width, height};
    }

    std::vector<std::string> paths;
    double start_x{}, start_y{}, width{}, height{};

private:
    void preprocess() {
        str_replace(svg_in, "  ", " ");
        str_replace(svg_in, "\n", "");
    }

    void parse_paths() {

        const std::regex path_regex{R"(<path(?:[\s\S]*?)d=\"([\s\S]*?)\"(?:[\s\S]*?)\/>)"};
        auto it = std::sregex_iterator(svg_in.begin(), svg_in.end(), path_regex);
        for (auto &match = it; match != std::sregex_iterator(); match++) {
            if (match->size() > 1) {
                paths.emplace_back(*(match->begin() + 1));
            }
        }
    }

    void parse_dims() {
        const std::regex viewbox_regex{
                R"(svg(?:[\s\S]*?)viewBox=\"((-?[\d.]+)? ([-?\d.]+)? ([-?\d.]+)? ([-?\d.]+)?)\"(?:[\s\S]*?))"
        };
        auto m = std::sregex_iterator(svg_in.begin(), svg_in.end(), viewbox_regex);
        if (m->size() >= 5) {
            try {
                start_x = std::stoi(*(m->begin() + 2));
                start_y = std::stoi(*(m->begin() + 3));
                width = std::stoi(*(m->begin() + 4));
                height = std::stoi(*(m->begin() + 5));
            }
            catch (std::invalid_argument &) {
                std::cout << "could not parse height and width";
                std::abort();
            }
        } else {
            std::cout << "could not parse dimensions data";
            std::abort();
        }
    }
};

class BasicPathWriter {
public:

    std::string path_prefx;

    BasicPathWriter(double width, double height) : BasicPathWriter{0, 0, width, height} {}

    BasicPathWriter(double start_x, double start_y, double width, double height) :
            start_x(start_x), start_y(start_y), width(width), height(height) {
        out = std::ostringstream{};
        put_prefix();
    }

    explicit BasicPathWriter(const std::tuple<double, double, double, double> &dims) :
            BasicPathWriter{std::get<0>(dims), std::get<1>(dims),
                            std::get<2>(dims), std::get<3>(dims)} {}

    void push_path_str(const std::string &path, const std::string &fillcolor = "lightgrey",
                       const std::string &stroke = "black", double factor = 1., double opacity = 0.6) {
        out << R"(<path d=")" << path << "\" "
            << R"(style="fill:)" << fillcolor << ";"
            << R"(stroke:)" << stroke << ";"
            << R"(stroke-width:)" << factor * (std::min(height, width) / 1000.0) << R"(px;)"
            << R"(opacity:)" << opacity << ";"
            << R"("/>)";
        out << '\n';
    }

    void push_circle(const double cx, const double cy, const double r, const std::string &color) {
        out << "<circle cx=\"" << cx << "\" cy=\""
            << cy << "\" r=\"" << r << "\" fill=\"" << color << "\"/>";
        out << '\n';
    }

    void push_circle(const double cx, const double cy) {
        push_circle(cx, cy, std::min(height, width) / 400, "black");
    }

    void push_circle(const double cx, const double cy, const std::string &color) {
        push_circle(cx, cy, std::min(height, width) / 400, color);
    }

    void push_text(const contourklip::Point2d &p, const std::string &text, double r = 0.01) {
        out << R"(<text x=")" << p.x() << R"(px" y=")" << p.y()
            << R"( px" style="font-family:'Helvetica';font-size:)"
            << r * std::min(height, width)
            << R"(px; fill:)" << "grey" << R"(;">)" << text << R"(</text>)";
    }

    void write_to(const std::string &filename) {
//        std::cout << "writing to file, contents: " << out.rdbuf();
        std::ofstream outfile{path_prefx + filename};
        outfile << out.rdbuf()->str();
        outfile << get_suffix();
    }

private:
    void put_prefix() {
        std::string t = R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)";
        std::string t2 = R"(<svg width="100%" height="100%" viewBox=")";
        std::string t3 = R"(" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linecap:round;">)";
        out << t << '\n' << t2;
        out << start_x << " " << start_y << " " << width << " " << height << t3 << '\n';
    }

    static std::string get_suffix() {
        return "</svg>";
    }

    double start_x, start_y, width, height;
    std::ostringstream out{};
};

#endif //CONTOURKLIP_SVG_IO_HPP