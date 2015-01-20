#include "quickhull.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <deque>
#include <set>

#include <cstdlib>
#include <cstdio>
#include <cassert>

int
main(int argc, char * argv[])
{
    using size_type = std::size_t;

    std::ifstream ifs_;
    if (!(argc < 2)) {
        ifs_.open(argv[1]);
        if (!ifs_.is_open()) {
            std::cout << std::flush;
            std::cerr << "cannot open the file" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::istream & in_ = ifs_.is_open() ? ifs_ : std::cin;
    std::cout << "#read file: " << ((argc < 2) ? "stdin" : argv[1]) << '\n';
    std::string line_;
    if (!std::getline(in_, line_)) {
        std::cout << std::flush;
        std::cerr << "io: dimension line" << std::endl;
        return EXIT_FAILURE;
    }
    std::istringstream iss_;
    size_type dimension_ = 0;
    {
        iss_.str(line_);
        if (!(iss_ >> dimension_)) {
            std::cout << std::flush;
            std::cerr << "io: dimension" << std::endl;
            return EXIT_FAILURE;
        }
        {
            using char_type = typename std::string::value_type;
            std::cout << "#command line:";
            std::istreambuf_iterator< char_type > const ibeg(iss_), iend;
            std::copy(ibeg, iend, std::ostreambuf_iterator< char_type >(std::cout));
            std::cout << '\n';
        }
        iss_.clear();
    }
    if (!std::getline(in_, line_)) {
        std::cout << std::flush;
        std::cerr << "io: count line" << std::endl;
        return EXIT_FAILURE;
    }
    using value_type = float;
    using point = std::valarray< value_type >;
    using points = std::vector< point >;
    size_type count_ = 0;
    {
        iss_.str(line_);
        if (!(iss_ >> count_)) {
            std::cout << std::flush;
            std::cerr << "io: count" << std::endl;
            return EXIT_FAILURE;
        }
        iss_.clear();
    }
    if (!(dimension_ < count_)) {
        std::cout << std::flush;
        std::cerr << "io: points count is less than or equal than dimensionality" << std::endl;
        return EXIT_FAILURE;
    }
    points points_(count_);
    for (size_type i = 0; i < count_; ++i) {
        if (!std::getline(in_, line_)) {
            std::cout << std::flush;
            std::cerr << "io: line count error" << std::endl;
            return EXIT_FAILURE;
        }
        point & point_ = points_[i];
        point_.resize(dimension_);
        {
            iss_.str(line_);
            for (size_type j = 0; j < dimension_; ++j) {
                if (!(iss_ >> point_[j])) {
                    std::cout << std::flush;
                    std::cerr << "io: bad value at line " << j << " of data" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            iss_.clear();
        }
    }
    //std::cout.rdbuf()->pubsetbuf(nullptr, 0);
    std::cout << "#D = " << dimension_ << '\n';
    std::cout << "#N = " << count_ << '\n';
    using quick_hull_type = quick_hull< typename points::const_iterator >;
    quick_hull_type quick_hull_(dimension_);
    typename quick_hull_type::point_list initial_simplex_;
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        {
            steady_clock::time_point const start = steady_clock::now();
            initial_simplex_ = quick_hull_.create_simplex(std::cbegin(points_), std::cend(points_));
            size_type const basis_size_ = initial_simplex_.size();
            steady_clock::time_point const end = steady_clock::now();
            std::cout << "#simplex time = " << duration_cast< microseconds >(end - start).count() << "us\n";
            if (basis_size_ != dimension_ + 1) {
                std::cout << std::flush;
                std::cerr << "cannot create a simplex: size of basis: " << basis_size_ << std::endl;
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
    os_ << "set view equal xyz\n";
    switch (dimension_) {
    case 2 : {
        os_ << "plot";
        break;
    }
    case 3 : {
        os_ << "splot";
        break;
    }
    default : {
        std::cout << std::flush;
        std::cerr << "dimensionality value (" << dimension_ << ") is out of supported range: cannot generate output" << std::endl;
        return EXIT_FAILURE;
    }
    }
    os_ << " '-' with points notitle pointtype 7 pointsize 1.5 linetype 1, '-' with points notitle, '-' with labels offset character 0, character 1 notitle";
    for (size_type i = 0; i < facets_count_; ++i) {
        os_ << ", '-' with lines notitle";
    }
    os_ << ";\n";
    for (auto const p : initial_simplex_) {
        point const & point_ = *p;
        for (value_type const & coordinate_ : point_) {
            os_ << coordinate_ << ' ';
        }
        os_ << '\n';
    }
    os_ << "e\n";
    for (size_type i = 0; i < count_; ++i) {
        point const & point_ = points_[i];
        for (value_type const & coordinate_ : point_) {
            os_ << coordinate_ << ' ';
        }
        os_ << '\n';
    }
    os_ << "e\n";
    for (size_type i = 0; i < count_; ++i) {
        point const & point_ = points_[i];
        for (value_type const & coordinate_ : point_) {
            os_ << coordinate_ << ' ';
        }
        os_ << i << '\n';
    }
    os_ << "e\n";
    for (size_type i = 0; i < facets_count_; ++i) {
        auto const & vertices_ = facets_[i].vertices_;
        for (auto const vertex_ : vertices_) {
            for (value_type const & coordinate_ : *vertex_) {
                os_ << coordinate_ << ' ';
            }
            os_ << '\n';
        }
        point const & first_vertex_ = *vertices_.front();
        for (value_type const & coordinate_ : first_vertex_) {
            os_ << coordinate_ << ' ';
        }
        os_ << "\ne\n";
    }
    os_ << std::flush;
    return EXIT_SUCCESS;
}
