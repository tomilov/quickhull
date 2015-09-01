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
#include <forward_list>
#include <iterator>
#include <utility>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstddef>

#ifdef __linux__
#define RED(str) __extension__ "\e[1;31m" str "\e[0m"
#else
#define RED(str) str
#endif

namespace
{

template< typename iterator >
struct indexed_iterator
{

    static_assert(std::is_base_of< std::forward_iterator_tag, typename std::iterator_traits< iterator >::iterator_category >::value);
    using iterator_category = std::random_access_iterator_tag; // not really true, but below machinery is enough for quickhull algorithm

    using value_type        = typename std::iterator_traits< iterator >::value_type;
    using difference_type   = typename std::iterator_traits< iterator >::difference_type;
    using pointer           = typename std::iterator_traits< iterator >::pointer;
    using reference         = typename std::iterator_traits< iterator >::reference;

    iterator base;
    std::size_t index;

    indexed_iterator &
    operator ++ ()
    {
        ++base; ++index;
        return *this;
    }

    bool
    operator < (indexed_iterator const & r) const
    {
        if (operator == (r)) {
            return false;
        }
        return (index < r.index);
    }

    difference_type
    operator - (indexed_iterator const & r) const
    {
        if (index < r.index) {
            return -static_cast< difference_type >(r.index - index);
        } else {
            return +static_cast< difference_type >(index - r.index);
        }
    }

    reference
    operator * () const
    {
        return *base;
    }

    bool
    operator == (indexed_iterator const & r) const
    {
        return (base == r.base);
    }

    bool
    operator != (indexed_iterator const & r) const
    {
        return !operator == (r);
    }

};

}

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
#if 0
    using value_type = float;
    using point = std::vector< value_type >;
    using points = std::vector< point >;
#else
    using value_type = double;
    using point = std::forward_list< value_type >;
    using points = std::forward_list< point >;
#endif
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
    using iterator = typename points::iterator;
    using point_iterator = std::conditional_t< (std::is_base_of< std::random_access_iterator_tag, iterator >{}), iterator, indexed_iterator< iterator > >;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
    auto const pbeg = point_iterator{std::begin(points_)};
    auto const pend = point_iterator{std::end(points_)};
#pragma clang diagnostic pop
    auto p = pbeg;
    for (size_type i = 0; i < count_; ++i) {
        if (!std::getline(in_, line_)) {
            out_ << std::flush;
            err_ << "io: line count error" << std::endl;
            return EXIT_FAILURE;
        }
        point & point_ = *p;
        ++p;
        point_.resize(dimension_);
        {
            iss_.str(line_);
            auto c = std::begin(point_);
            for (size_type j = 0; j < dimension_; ++j) {
                if (!(iss_ >> *c)) {
                    out_ << std::flush;
                    err_ << "io: bad value at line " << j << " of data" << std::endl;
                    return EXIT_FAILURE;
                }
                ++c;
            }
            iss_.clear();
        }
    }
    assert(p == pend);
    //out_.rdbuf()->pubsetbuf(nullptr, 0);
    out_ << "#D = " << dimension_ << '\n';
    out_ << "#N = " << count_ << '\n';
    using quick_hull_type = quick_hull< point_iterator >;
#if 0
    using std::sqrt;
    value_type const eps = sqrt(std::numeric_limits< value_type >::epsilon()); // do use relaxed constraints for input, which supposedly should produce a plenty of complanar facets, like 'rbox D3 27 M3,4'
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
            initial_simplex_ = quick_hull_.create_initial_simplex(pbeg, pend);
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
    auto const facets_count_ = std::distance(std::cbegin(facets_), std::cend(facets_));
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
    for (auto i = facets_count_; 0 < i; --i) {
        gnuplot_ << ", '-' with lines notitle"
                    ", '-' with points notitle pointtype 6 pointsize 1.5 linetype 4";
    }
    gnuplot_ << ";\n";
    for (auto const v : initial_simplex_) {
        point const & point_ = *v;
        for (value_type const & coordinate_ : point_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << '\n';
    }
    gnuplot_ << "e\n";
    p = pbeg;
    for (size_type i = 0; i < count_; ++i) {
        point const & point_ = *p;
        ++p;
        for (value_type const & coordinate_ : point_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << '\n';
    }
    assert(p == pend);
    gnuplot_ << "e\n";
    p = pbeg;
    for (size_type i = 0; i < count_; ++i) {
        point const & point_ = *p;
        ++p;
        for (value_type const & coordinate_ : point_) {
            gnuplot_ << coordinate_ << ' ';
        }
        gnuplot_ << i << '\n';
    }
    assert(p == pend);
    gnuplot_ << "e\n";
    for (auto const & facet_ : facets_) {
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
        for (auto const v : facet_.coplanar_) {
            for (value_type const & coordinate_ : *v) {
                gnuplot_ << coordinate_ << ' ';
            }
            gnuplot_ << '\n';
        }
        gnuplot_ << "e\n";
    }
    gnuplot_ << std::flush;
    return EXIT_SUCCESS;
}
