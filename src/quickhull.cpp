#include "quickhull.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <valarray>
#include <chrono>

#include <cstdlib>
#include <cstdio>

#include <iostream>

#include <cstdlib>
#include <cassert>

int
main(int argc, char * argv[])
{
    if (argc < 2) {
        std::cerr << "error: argc == "  << argc << std::endl;
        return EXIT_FAILURE;
    }
    using size_type = std::size_t;

    std::ifstream ifs_;
    ifs_.open(argv[1]); // rbox n D3 s 100 > points.txt
    if (!ifs_.is_open()) {
        std::cerr << "file is not open" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "read file: " << argv[1] << std::endl;
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
    using G = long double;
    using point_type = std::valarray< G >;
    using points_type = std::valarray< point_type >;
    size_type const count_ = std::stoll(line_);
    points_type points_(count_);
    std::istringstream iss;
    for (size_type i = 0; i < count_; ++i) {
        if (!std::getline(ifs_, line_)) {
            std::cerr << "io: line count error" << std::endl;
            return EXIT_FAILURE;
        }
        point_type & point_ = points_[i];
        point_.resize(dim_);
        iss.str(line_);
        std::copy_n(std::istream_iterator< G >(iss), dim_, std::begin(point_));
        if (!iss) {
            std::cerr << "bad value at line " << points_.size() << " of data" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::cout.rdbuf()->pubsetbuf(nullptr, 0);
    std::cout << "D = " << dim_ << std::endl;
    std::cout << "N = " << count_ << std::endl;
    using H = quick_hull< points_type >;
    H quick_hull_(dim_, points_);
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        {
            steady_clock::time_point const start = steady_clock::now();
            size_type const basis_size_ = quick_hull_.create_simplex().size();
            steady_clock::time_point const end = steady_clock::now();
            std::cout << "simplex time = " << duration_cast< microseconds >(end - start).count() << "us" << std::endl;
            if (basis_size_ != dim_ + 1) {
                std::cerr << "cant create a simplex" << std::endl;
                return EXIT_FAILURE;
            }
        }
        {
            steady_clock::time_point const start = steady_clock::now();
            quick_hull_.create_convex_hull();
            steady_clock::time_point const end = steady_clock::now();
            std::cout << "qh time = " << duration_cast< microseconds >(end - start).count() << "us" << std::endl;
        }
    }
    auto const & facets_ = quick_hull_.facets_;
    size_type const facets_count_ = facets_.size();
    std::cout << "number of facets: " << facets_count_ << std::endl;
#if 1
    std::ofstream ofs_;
    ofs_.open("script.txt"); // gnuplot> load 'script.txt'
    if (!ofs_.is_open()) {
        std::cerr << "output file cannot be truncated" << std::endl;
        return EXIT_FAILURE;
    }
    ofs_ << "clear" << std::endl;
    ofs_ << "set autoscale" << std::endl;
    switch (dim_) {
    case 1 : {
        ofs_ << "plot";
        break;
    }
    case 2 : {
        ofs_ << "plot";
        break;
    }
    case 3 : {
        ofs_ << "splot";
        break;
    }
    default : {
        std::cerr << "dimensionality value (" << dim_ << ") is out of supported range" << std::endl;
        return EXIT_FAILURE;
    }
    }
    ofs_ << " '-' with points notitle, '-' with labels offset character 0,character 1 notitle";
    for (std::size_t i = 0; i < facets_count_; ++i) {
        ofs_ << ", '-' with lines notitle";
    }
    ofs_ << ';' << std::endl;
    for (std::size_t i = 0; i < count_; ++i) {
        point_type const & point_ = points_[i];
        for (G const & coordinate_ : point_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << std::endl;
    }
    ofs_ << 'e' << std::endl;
    for (std::size_t i = 0; i < count_; ++i) {
        point_type const & point_ = points_[i];
        for (G const & coordinate_ : point_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << i << std::endl;
    }
    ofs_ << 'e' << std::endl;
    for (size_type i = 0; i < facets_count_; ++i) {
        auto const & vertices_ = facets_[i].vertices_;
        for (size_type const vertex_ : vertices_) {
            for (G const & coordinate_ : points_[vertex_]) {
                ofs_ << coordinate_ << ' ';
            }
            ofs_ << std::endl;
        }
        point_type const & first_vertex_ = points_[vertices_.front()];
        for (G const & coordinate_ : first_vertex_) {
            ofs_ << coordinate_ << ' ';
        }
        ofs_ << std::endl;
        ofs_ << 'e' << std::endl;
    }
#endif
    return EXIT_SUCCESS;
}
