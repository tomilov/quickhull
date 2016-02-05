#ifdef _DEBUG
#include <iostream>
#endif
#include <quickhull.hpp>

#include <limits>
#include <iterator>
#include <algorithm>
#include <valarray>
#include <vector>
#include <forward_list>
#include <string>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>

#include <cmath>
#include <cstdlib>

#ifdef __linux__
#define TERM_COLOR(code)   __extension__ "\e[1;" #code "m"
#define TERM_COLOR_RED     TERM_COLOR(31)
#define TERM_COLOR_GREEN   TERM_COLOR(32)
#define TERM_COLOR_BLUE    TERM_COLOR(34)
#define TERM_COLOR_DEFAULT __extension__ "\e[0m"
#else
#define TERM_COLOR_RED     ""
#define TERM_COLOR_GREEN   ""
#define TERM_COLOR_BLUE    ""
#define TERM_COLOR_DEFAULT ""
#endif

int
main(int argc, char * argv[]) // rbox D3 t 100 | bin/qh | gnuplot -p
{
    std::ostream & err_ = std::cerr;
    std::ostream & log_ = std::clog;

    // choose input source
    std::ifstream ifs_;
    if (argc == 2) {
        ifs_.open(argv[1]);
        if (!ifs_.is_open()) {
            err_ << "error: cannot open file '" << argv[1] << "'" << std::endl;
            return EXIT_FAILURE;
        }
    }
    std::istream & in_ = (ifs_.is_open() ? ifs_ : std::cin);
    using size_type = std::size_t;

    // select type
#if 0
    // RandomAccessIterator
    using value_type = float;
    using point = std::valarray< value_type >;
    using points = std::vector< point >;
#else
    // ForwardIterator
    using value_type = double;
    using point = std::forward_list< value_type >;
    using points = std::forward_list< point >;
#endif

    std::string line_;
    std::istringstream iss_;

    // fill container with points
    size_type dimension_ = 0;
    {
        if (!std::getline(in_, line_)) {
            err_ << "error: input: missing dimension line" << std::endl;
            return false;
        }
        iss_.str(line_);
        if (!(iss_ >> dimension_)) {
            err_ << "error: input: dimension format" << std::endl;
            return false;
        }
        log_ << "dimensionality of input is " << dimension_ << std::endl;
        if (!(1 < dimension_)) {
            err_ << "error: input: dimensionality value is not greater then one" << std::endl;
            return false;
        }
        {
            using char_type = typename std::istream::char_type;
            log_ << "rbox command line:";
            std::istreambuf_iterator< char_type > const ibeg(iss_), iend;
            std::copy(ibeg, iend, std::ostreambuf_iterator< char_type >(log_));
            log_ << std::endl;
        }
        iss_.clear();
    }
    size_type count_ = 0;
    {
        if (!std::getline(in_, line_)) {
            err_ << "error: input: missing count line" << std::endl;
            return false;
        }
        iss_.str(line_);
        if (!(iss_ >> count_)) {
            err_ << "error: input: format of count" << std::endl;
            return false;
        }
        iss_.clear();
        log_ << "input points count = " << count_ << std::endl;
        if (!(dimension_ < count_)) {
            err_ << "error: input: points count is not greater then dimensionality" << std::endl;
            return false;
        }
    }
    points points_(count_);
    for (point & point_ : points_) {
        if (!std::getline(in_, line_)) {
            err_ << "error: input: wrong line count" << std::endl;
            return false;
        }
        point_.resize(dimension_);
        {
            iss_.str(line_);
            auto c = std::begin(point_);
            for (size_type j = 0; j < dimension_; ++j) {
                if (!(iss_ >> *c)) {
                    err_ << "error: input: bad corodinate value at line " << j << " of data" << std::endl;
                    return false;
                }
                ++c;
            }
            iss_.clear();
        }
    }

    // set epsilon (can be zero)
    //value_type const zero = value_type(0);
    value_type const eps = std::numeric_limits< value_type >::epsilon();
    log_ << "epsilon = " << eps << std::endl;

    // define and setup QH class instance
    using quick_hull_type = quick_hull< typename points::const_iterator >;
    quick_hull_type quick_hull_(dimension_, eps); // (1)
    quick_hull_.add_points(std::cbegin(points_), std::cend(points_)); // (2)
    auto const initial_simplex_ = quick_hull_.get_affine_basis(); // (3)

    // run the algorithm
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
    using std::chrono::steady_clock;
    { // create initial simplex
        steady_clock::time_point const start = steady_clock::now();
        quick_hull_.create_initial_simplex(std::cbegin(initial_simplex_),
                                           std::prev(std::cend(initial_simplex_))); // (4)
        auto const delta = duration_cast< microseconds >(steady_clock::now() - start).count();
        log_ << "simplex time = " << delta << "us" << std::endl;
        size_type const basis_size_ = initial_simplex_.size();
        if (basis_size_ != quick_hull_.dimension_ + 1) { // (5)
            err_ << "error: algorithm: cannot construct a simplex. Degenerated input set. Size of basis: "
                      << basis_size_ << std::endl;
            return false;
        }
    }
    { // create convex hull
        steady_clock::time_point const start = steady_clock::now();
        quick_hull_.create_convex_hull(); // (6)
        auto const delta = duration_cast< microseconds >(steady_clock::now() - start).count();
        log_ << "quickhull time = "
                  << TERM_COLOR_GREEN << delta << "us"
                  << TERM_COLOR_DEFAULT << std::endl;
    }
    log_ << "number of (convex hull) polyhedron facets is "
              << TERM_COLOR_BLUE << quick_hull_.facets_.size()
              << TERM_COLOR_DEFAULT << std::endl;
    if (!quick_hull_.check()) {
        err_ << TERM_COLOR_RED << "error: algorithm: resulting structure is not valid convex polytope"
                  << TERM_COLOR_DEFAULT << std::endl;
        return false;
    }

    // output
    std::ostream & gnuplot_ = std::cout;
    if (3 < quick_hull_.dimension_) {
        log_ << "dimensionality value " << quick_hull_.dimension_
                  << " is out of supported range: cannot generate output for this" << std::endl;
        return EXIT_SUCCESS;
    }
    gnuplot_ << "set view equal xyz\n"
                "set autoscale\n"
                "set key left\n"
                "set xrange [] writeback\n"
                "set yrange [] writeback\n"
                "set zrange [] writeback\n";
    gnuplot_ << "set title \'Points count is " << count_ << "\'\n";
    if (quick_hull_.dimension_ == 2) {
        gnuplot_ << "plot";
    } else if (quick_hull_.dimension_ == 3) {
        gnuplot_ << "splot";
    }
    gnuplot_ << " '-' with points notitle pointtype 4 pointsize 1.5 linetype 1"
                ", '-' with points notitle"
                ", '-' with labels offset character 0, character 1 notitle";
    for (auto const & facet_ : quick_hull_.facets_) {
        gnuplot_ << ", '-' with lines notitle";
        if (!facet_.coplanar_.empty()) {
            gnuplot_ << ", '-' with points notitle pointtype 6 pointsize 1.5 linetype 4";
        }
    }
    gnuplot_ << ";\n";
    {
        for (auto const & v : initial_simplex_) {
            point const & point_ = *v;
            for (value_type const & coordinate_ : point_) {
                gnuplot_ << coordinate_ << ' ';
            }
            gnuplot_ << '\n';
        }
        gnuplot_ << "e\n";
        {
            for (auto const & point_ : points_) {
                for (value_type const & coordinate_ : point_) {
                    gnuplot_ << coordinate_ << ' ';
                }
                gnuplot_ << '\n';
            }
            gnuplot_ << "e\n";
            size_type i = 0;
            for (auto const & point_ : points_) {
                for (value_type const & coordinate_ : point_) {
                    gnuplot_ << coordinate_ << ' ';
                }
                gnuplot_ << i << '\n';
                ++i;
            }
        }
        gnuplot_ << "e\n";
    }
    for (auto const & facet_ : quick_hull_.facets_) {
        auto const & vertices_ = facet_.vertices_;
        for (auto const & vertex_ : vertices_) {
            for (value_type const & coordinate_ : *vertex_) {
                gnuplot_ << coordinate_ << ' ';
            }
            gnuplot_ << '\n';
        }
        for (value_type const & coordinate_ : *vertices_.front()) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << "\n"
                    "e\n";
        if (!facet_.coplanar_.empty()) {
            for (auto const & v : facet_.coplanar_) {
                for (value_type const & coordinate_ : *v) {
                    gnuplot_ << coordinate_ << ' ';
                }
                gnuplot_ << '\n';
            }
            gnuplot_ << "e\n";
        }
    }
    gnuplot_ << std::flush;
    return EXIT_SUCCESS;
}
