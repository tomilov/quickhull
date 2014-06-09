#include "quickhull.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <cstdlib>

int main()
{
    using G = double;
    using H = convex_hull< G >;
    using point = typename H::point_type;

    std::ifstream ifs_;
    ifs_.open("../10.txt");
    if (ifs_.is_open()) {
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
        H convex_hull_(points_.cbegin(), points_.cend());
        convex_hull_.create_simplex();
        for (auto const & facet_ : convex_hull_.facets_) {
            for (point const & point_ : facet_.second.vertices_) {
                for (G const & coordinate_ : point_) {
                    std::cout << coordinate_ << ' ';
                }
                std::cout << std::endl;
            }
            point const & first_point_ = facet_.second.vertices_.front();
            for (G const & coordinate_ : first_point_) {
                std::cout << coordinate_ << ' ';
            }
            std::cout << std::endl << std::endl;
        }
    } else {
        std::cerr << "file is not open" << std::endl;
    }
    return EXIT_SUCCESS;
}
