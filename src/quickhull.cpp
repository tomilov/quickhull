#include "quickhull.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <cstdlib>
#include <cstdio>

int main()
{
    using G = double;
    using H = convex_hull< G >;
    using point = typename H::point_type;

    /*{
        size_type const COUNT = 100;
        std::stringstream ss;
        ss << "bash -c 'rbox n D3 " << COUNT << " | tail -" << COUNT << " &> points.txt'";
        std::system(ss.str().c_str());
    }*/
    std::ifstream ifs_;
    ifs_.open("points.txt"); // rbox n D3 100 > points.txt
    if (!ifs_.is_open()) {
        std::cerr << "file is not open" << std::endl;
        return EXIT_FAILURE;
    }
    std::string line_;
    if (!std::getline(ifs_, line_)) {
        std::cerr << "no dim at first line" << std::endl;
        return EXIT_FAILURE;
    }
    size_type const dim_ = std::stoll(line_);
    if (!std::getline(ifs_, line_)) {
        std::cerr << "no count at second line" << std::endl;
        return EXIT_FAILURE;
    }
    size_type const count_ = std::stoll(line_);
    std::deque< point > points_;
    while (std::getline(ifs_, line_)) {
        point point_(dim_);
        std::istringstream iss(line_);
        for (size_type i = 0; i < dim_; ++i) {
            if (!(iss >> point_[i])) {
                std::cerr << "bad value faced at line " << points_.size() << " of data" << std::endl;
                return EXIT_FAILURE;
            }
        }
        points_.push_back(std::move(point_));
    }
    if (count_ != points_.size()) {
        std::cerr << "input file format error" << std::endl;
        return EXIT_FAILURE;
    }
    std::ofstream ofs_;
    ofs_.open("script.txt");
    if (!ofs_.is_open()) {
        std::cerr << "output file cannot be truncated" << std::endl;
        return EXIT_FAILURE;
    }
    ofs_ << "set view equal xyz; set view 0,0; set xrange [-0.5:0.5]; set yrange [-0.5:0.5]; set zrange [-0.5:0.5]; set xyplane at -0.5" << std::endl;
    ofs_ << "splot '-' w p, '-' w l" << std::endl;
    for (point const & point_ : points_) {
        for (G const & coordinate_ : point_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << std::endl;
    }
    ofs_ << 'e' << std::endl;
    if (!ofs_.is_open()) {
        std::cerr << "output file cannot be truncated" << std::endl;
        return EXIT_FAILURE;
    }
    H convex_hull_(points_.cbegin(), points_.cend());
    convex_hull_.create_simplex();
    for (auto const & facet_ : convex_hull_.facets_) {
        for (point const & point_ : facet_.second.vertices_) {
            for (G const & coordinate_ : point_) {
                ofs_ << coordinate_ << ' ';
            }
            ofs_ << std::endl;
        }
        point const & first_point_ = facet_.second.vertices_.front();
        for (G const & coordinate_ : first_point_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << std::endl << std::endl;
    }
    ofs_ << 'e' << std::endl;
    return EXIT_SUCCESS;
}
