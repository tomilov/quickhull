#include "quickhull.hpp"

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <set>
#include <numeric>
#include <valarray>
#include <vector>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>

#ifdef __linux__
#define RED(str) __extension__ "\e[1;31m" str "\e[0m"
#else
#define RED(str) str
#endif

int
main(int argc, char * argv[])
{
    using size_type = std::size_t;

    std::ostream & out_ = std::cout;
    std::ostream & err_ = std::cerr;

    std::ifstream ifs_;
    if (argc == 2) {
        ifs_.open(argv[1]);
        if (!ifs_.is_open()) {
            out_ << std::flush;
            err_ << "cannot open file" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::istream & in_ = ifs_.is_open() ? ifs_ : std::cin;

    in_.sync_with_stdio(false);
    out_.sync_with_stdio(false);
    err_.sync_with_stdio(false);

    out_ << "#input file: " << ((argc < 2) ? "stdin" : argv[1]) << '\n';
    std::string line_;
    if (!std::getline(in_, line_)) {
        out_ << std::flush;
        err_ << "io: dimension line" << std::endl;
        return EXIT_FAILURE;
    }
    std::istringstream iss_;
    size_type dimension_ = 0;
    {
        iss_.str(line_);
        if (!(iss_ >> dimension_)) {
            out_ << std::flush;
            err_ << "io: dimension" << std::endl;
            return EXIT_FAILURE;
        }
        {
            using char_type = typename std::string::value_type;
            out_ << "#command line:";
            std::istreambuf_iterator< char_type > const ibeg(iss_), iend;
            std::copy(ibeg, iend, std::ostreambuf_iterator< char_type >(out_));
            out_ << '\n';
        }
        iss_.clear();
    }
    if (!std::getline(in_, line_)) {
        out_ << std::flush;
        err_ << "io: count line" << std::endl;
        return EXIT_FAILURE;
    }
    using value_type = double;
    using point = std::valarray< value_type >;
    using points = std::vector< point >;
    size_type count_ = 0;
    {
        iss_.str(line_);
        if (!(iss_ >> count_)) {
            out_ << std::flush;
            err_ << "io: count" << std::endl;
            return EXIT_FAILURE;
        }
        iss_.clear();
    }
    if (!(dimension_ < count_)) {
        out_ << std::flush;
        err_ << "io: points count is less than or equal than dimensionality" << std::endl;
        return EXIT_FAILURE;
    }
    points points_(count_);
    for (size_type i = 0; i < count_; ++i) {
        if (!std::getline(in_, line_)) {
            out_ << std::flush;
            err_ << "io: line count error" << std::endl;
            return EXIT_FAILURE;
        }
        point & point_ = points_[i];
        point_.resize(dimension_);
        {
            iss_.str(line_);
            for (size_type j = 0; j < dimension_; ++j) {
                if (!(iss_ >> point_[j])) {
                    out_ << std::flush;
                    err_ << "io: bad value at line " << j << " of data" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            iss_.clear();
        }
    }
    //out_.rdbuf()->pubsetbuf(nullptr, 0);
    out_ << "#D = " << dimension_ << '\n';
    out_ << "#N = " << count_ << '\n';
    using quick_hull_type = quick_hull< typename points::const_iterator >;
#if 0
    using std::sqrt;
    value_type const eps = sqrt(std::numeric_limits< value_type >::epsilon()); // use relaxed constraints for input, which supposedly should produce a plenty of complanar facets, like 'rbox D3 27 M3,4'
#else
    value_type const eps = std::numeric_limits< value_type >::epsilon();
#endif
    quick_hull_type quick_hull_(dimension_, eps);
    typename quick_hull_type::point_array initial_simplex_;
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        {
            steady_clock::time_point const start = steady_clock::now();
            initial_simplex_ = quick_hull_.create_initial_simplex(std::cbegin(points_), std::cend(points_));
            size_type const basis_size_ = initial_simplex_.size();
            steady_clock::time_point const end = steady_clock::now();
            out_ << "#simplex time = " << duration_cast< microseconds >(end - start).count() << "us\n";
            if (basis_size_ != dimension_ + 1) {
                out_ << std::flush;
                err_ << "cannot create a simplex: size of basis: " << basis_size_ << std::endl;
                return EXIT_FAILURE;
            }
        }
        {
            steady_clock::time_point const start = steady_clock::now();
            quick_hull_.create_convex_hull();
            steady_clock::time_point const end = steady_clock::now();
            out_ << "#quickhull time = " << duration_cast< microseconds >(end - start).count() << "us\n";
            if (!quick_hull_.check()) {
                out_ << std::flush;
                err_ << RED("resulting structure is not valid convex polytope") << std::endl;
                return EXIT_FAILURE;
            }
        }
    }
    auto const & facets_ = quick_hull_.facets_;
    size_type const facets_count_ = facets_.size();
    out_ << "#number of facets: " << facets_count_ << std::endl;
    std::ostream & gnuplot_ = out_;
    gnuplot_ << "clear\n";
    gnuplot_ << "set autoscale\n";
    gnuplot_ << "set view equal xyz\n";
    switch (dimension_) {
    case 2 : {
        gnuplot_ << "plot";
        break;
    }
    case 3 : {
        gnuplot_ << "splot";
        break;
    }
    default : {
        out_ << std::flush;
        err_ << "dimensionality value (" << dimension_ << ") is out of supported range: cannot generate output" << std::endl;
        return EXIT_FAILURE;
    }
    }
    gnuplot_ << " '-' with points notitle pointtype 4 pointsize 1.5 linetype 1"
                ", '-' with points notitle"
                ", '-' with labels offset character 0, character 1 notitle";
    for (size_type i = 0; i < facets_count_; ++i) {
        gnuplot_ << ", '-' with lines notitle"
                    ", '-' with points notitle pointtype 6 pointsize 1.5 linetype 4";
    }
    gnuplot_ << ";\n";
    for (auto const p : initial_simplex_) {
        point const & point_ = *p;
        for (value_type const & coordinate_ : point_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << '\n';
    }
    gnuplot_ << "e\n";
    for (size_type i = 0; i < count_; ++i) {
        point const & point_ = points_[i];
        for (value_type const & coordinate_ : point_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << '\n';
    }
    gnuplot_ << "e\n";
    for (size_type i = 0; i < count_; ++i) {
        point const & point_ = points_[i];
        for (value_type const & coordinate_ : point_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << i << '\n';
    }
    gnuplot_ << "e\n";
    for (size_type i = 0; i < facets_count_; ++i) {
        auto const & facet_ = facets_[i];
        auto const & vertices_ = facet_.vertices_;
        for (auto const vertex_ : vertices_) {
            for (value_type const & coordinate_ : *vertex_) {
                gnuplot_ << coordinate_ << ' ';
            }
            gnuplot_ << '\n';
        }
        std::cerr << std::endl;
        point const & first_vertex_ = *vertices_.front();
        for (value_type const & coordinate_ : first_vertex_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << "\n";
        gnuplot_ << "e\n";
        for (auto const p : facet_.coplanar_) {
            for (value_type const & coordinate_ : *p) {
                gnuplot_ << coordinate_ << ' ';
            }
            gnuplot_ << '\n';
        }
        gnuplot_ << "e\n";
    }
    gnuplot_ << std::flush;
    return EXIT_SUCCESS;
}
