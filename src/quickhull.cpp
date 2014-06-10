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
    using point_type = typename H::point_type;

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
    std::deque< point_type > points_;
    while (std::getline(ifs_, line_)) {
        points_.emplace_back(dim_);
        point_type & point_ = points_.back();
        std::istringstream iss(line_);
        for (size_type i = 0; i < dim_; ++i) {
            if (!(iss >> point_[i])) {
                std::cerr << "bad value faced at line " << points_.size() << " of data" << std::endl;
                return EXIT_FAILURE;
            }
        }
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
    ofs_ << "reset" << std::endl;
    switch (dim_) {
    case 1 :
    case 2 : {
        ofs_ << "plot";
        break;
    }
    case 3 : {
        ofs_ << "set view equal xyz; set view 0,0; set xrange [-0.5:0.5]; set yrange [-0.5:0.5]; set zrange [-0.5:0.5]; set xyplane at -0.5" << std::endl;
        ofs_ << "splot";
        break;
    }
    default : {
        std::cerr << "dimensionality value (" << dim_ << ") is out of supported range" << std::endl;
        return EXIT_FAILURE;
    }
    }
    ofs_ << " '-' w p, '-' w l" << std::endl;
    for (point_type const & point_ : points_) {
        for (G const & coordinate_ : point_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << std::endl;
    }
    ofs_ << 'e' << std::endl;
    H convex_hull_(points_.cbegin(), points_.cend());
    convex_hull_.create_simplex();
    for (auto const & facet_ : convex_hull_.facets_) {
        auto const & vertices_ = facet_.second.vertices_;
        for (point_type const & point_ : vertices_) {
            for (G const & coordinate_ : point_) {
                ofs_ << coordinate_ << ' ';
            }
            ofs_ << std::endl;
        }
        point_type const & first_point_ = vertices_.front();
        for (G const & coordinate_ : first_point_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << std::endl << std::endl;
    }
    ofs_ << 'e' << std::endl;
    return EXIT_SUCCESS;
}
