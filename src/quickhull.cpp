#include "quickhull.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <valarray>
#include <chrono>

#include <cstdlib>
#include <cstdio>
#include <cassert>

int
main(int argc, char * argv[])
{
    using size_type = std::size_t;

    std::ifstream ifs_;
    if (!(argc < 2)) {
        ifs_.open(argv[1]); // rbox D3 s 100 > points.txt
        if (!ifs_.is_open()) {
            std::cerr << "file is not open" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::istream & in_ = (argc < 2) ? std::cin : ifs_;
    std::cout << "#read file: " << ((argc < 2) ? "stdin" : argv[1]) << '\n';
    std::string line_;
    if (!std::getline(in_, line_)) {
        std::cerr << "bad file format" << std::endl;
        return EXIT_FAILURE;
    }
    size_type dim_;
    std::istringstream iss(line_);
    iss >> dim_;
    {
        using char_type = typename std::string::value_type;
        std::cout << "#command line:";
        std::istreambuf_iterator< char_type > const ibeg(iss), iend;
        std::copy(ibeg, iend, std::ostreambuf_iterator< char_type >(std::cout));
        std::cout << '\n';
    }
    if (!std::getline(in_, line_)) {
        std::cerr << "no count at second line" << std::endl;
        return EXIT_FAILURE;
    }
    using G = float;
    using point_type = std::valarray< G >;
    using points_type = std::valarray< point_type >;
    size_type const count_ = std::stoll(line_);
    points_type points_(count_);
    for (size_type i = 0; i < count_; ++i) {
        if (!std::getline(in_, line_)) {
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
    //std::cout.rdbuf()->pubsetbuf(nullptr, 0);
    std::cout << "#D = " << dim_ << '\n';
    std::cout << "#N = " << count_ << '\n';
    using quick_hull_type = quick_hull< points_type >;
    quick_hull_type quick_hull_(dim_, points_);
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        {
            steady_clock::time_point const start = steady_clock::now();
            size_type const basis_size_ = quick_hull_.create_simplex().size();
            steady_clock::time_point const end = steady_clock::now();
            std::cout << "#simplex time = " << duration_cast< microseconds >(end - start).count() << "us\n";
            if (basis_size_ != dim_ + 1) {
                std::cerr << "cannot create a simplex" << std::endl;
                return EXIT_FAILURE;
            }
        }
        {
            steady_clock::time_point const start = steady_clock::now();
            quick_hull_.create_convex_hull();
            steady_clock::time_point const end = steady_clock::now();
            std::cout << "#quickhull time = " << duration_cast< microseconds >(end - start).count() << "us\n";
        }
    }
    auto const & facets_ = quick_hull_.facets_;
    size_type const facets_count_ = facets_.size();
    std::cout << "#number of facets: " << facets_count_ << std::endl;
    std::ostream & os_ = std::cout;
    os_ << "clear\n";
    os_ << "set autoscale\n";
    switch (dim_) {
    case 1 : {
        os_ << "plot";
        break;
    }
    case 2 : {
        os_ << "plot";
        break;
    }
    case 3 : {
        os_ << "splot";
        break;
    }
    default : {
        std::cerr << "dimensionality value (" << dim_ << ") is out of supported range: cannot generate output" << std::endl;
        return EXIT_FAILURE;
    }
    }
    os_ << " '-' with points notitle, '-' with labels offset character 0, character 1 notitle";
    for (std::size_t i = 0; i < facets_count_; ++i) {
        os_ << ", '-' with lines notitle";
    }
    os_ << ";\n";
    for (std::size_t i = 0; i < count_; ++i) {
        point_type const & point_ = points_[i];
        for (G const & coordinate_ : point_) {
            os_ << coordinate_ << ' ';
        }
        os_ << '\n';
    }
    os_ << "e\n";
    for (std::size_t i = 0; i < count_; ++i) {
        point_type const & point_ = points_[i];
        for (G const & coordinate_ : point_) {
            os_ << coordinate_ << ' ';
        }
        os_ << i << '\n';
    }
    os_ << "e\n";
    for (size_type i = 0; i < facets_count_; ++i) {
        auto const & vertices_ = facets_[i].vertices_;
        for (size_type const vertex_ : vertices_) {
            for (G const & coordinate_ : points_[vertex_]) {
                os_ << coordinate_ << ' ';
            }
            os_ << '\n';
        }
        point_type const & first_vertex_ = points_[vertices_.front()];
        for (G const & coordinate_ : first_vertex_) {
            os_ << coordinate_ << ' ';
        }
        os_ << "\ne\n";
    }
    os_ << std::endl;
    return EXIT_SUCCESS;
}
