#ifdef _DEBUG
#include <iostream>
#endif
#include <quickhull.hpp>

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
#include <algorithm>
#include <limits>
#include <memory>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <cstddef>

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

template< typename value_type = float,
          typename point = std::valarray< value_type >,
          typename points = std::vector< point > >
struct test_quickhull
{

    using size_type = std::size_t;

    std::ostream & err_;
    std::ostream & log_;

    test_quickhull(std::ostream & _err = std::cerr,
                   std::ostream & _log = std::clog)
        : err_(_err)
        , log_(_log)
    {
        err_.sync_with_stdio(false);
        log_.sync_with_stdio(false);
    }

    size_type dimension_ = 0;
    size_type count_ = 0;
    points points_;

    std::string line_;
    std::istringstream iss_;

    bool
    input(std::istream & _in)
    {
        {
            if (!std::getline(_in, line_)) {
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
        {
            if (!std::getline(_in, line_)) {
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
                err_ << "error: input: points count is not greater than to dimensionality" << std::endl;
                return false;
            }
        }
        points_ = points(count_);
        for (point & point_ : points_) {
            if (!std::getline(_in, line_)) {
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
        return true;
    }

    friend
    std::istream &
    operator >> (std::istream & _in, test_quickhull & _qh)
    {
        if (!_qh.input(_in)) {
            _in.setstate(std::ios::failbit);
        }
        return _in;
    }

    struct gnuplot
    {

        std::ostream & err_;
        std::ostream & log_;

#if 0
        using point_iterator_type = typename points::const_iterator;
#elif 0
        using point_iterator_type = typename points::pointer; // universe_.push_back(std::addressof(*p));
#else
        using point_iterator_type = typename points::iterator;
#endif

        using quick_hull_type = quick_hull< point_iterator_type >;


        value_type const zero = value_type(0);
#if 0
        value_type const eps = sqrt(std::numeric_limits< value_type >::epsilon());
#elif 0
        value_type const eps = std::numeric_limits< value_type >::epsilon();
#else
        value_type const eps = zero;
#endif

        quick_hull_type quick_hull_;

        typename quick_hull_type::point_list universe_;

        struct point_iterator_less
        {

            bool
            operator () (point_iterator_type const & _lhs,
                         point_iterator_type const & _rhs) const
            {
                return pointer_less_(std::addressof(*_lhs), std::addressof(*_rhs));
            }

        private :

            std::less< typename std::iterator_traits< point_iterator_type >::value_type const * > pointer_less_;

        };

        gnuplot(std::ostream & _err,
                std::ostream & _log,
                size_type const _dimension,
                points & _points)
            : err_(_err)
            , log_(_log)
            , quick_hull_(_dimension, eps)
        {
            auto p = std::begin(_points);
            auto const pend = std::end(_points);
            if (operator bool ()) {
                quick_hull_.add_points(p, pend);
            } else {
                while (p != pend) {
                    universe_.push_back(p);
                    ++p;
                }
                universe_.sort(point_iterator_less{});
            }
        }

        typename quick_hull_type::point_list initial_simplex_;

        bool
        operator () (bool const _use_simplex_heuristic = false)
        {
            assert(initial_simplex_.empty());
            assert(quick_hull_.facets_.empty());

            if (operator bool ()) {
                initial_simplex_ = quick_hull_.get_affine_basis();
            } else {
                log_ << "\ncount of points to process is " << universe_.size() << std::endl;
                auto u = std::cbegin(universe_);
                if (!_use_simplex_heuristic) {
                    for (size_type i = 0; i <= quick_hull_.dimension_; ++i) {
                        initial_simplex_.push_back(*u);
                        ++u;
                    }
                }
                quick_hull_.add_points(u, std::cend(universe_));
                if (_use_simplex_heuristic) {
                    initial_simplex_ = quick_hull_.get_affine_basis();
                }
            }

            log_ << "epsilon = " << eps << std::endl;

            using std::chrono::duration_cast;
            using std::chrono::microseconds;
            using std::chrono::steady_clock;
            {
                steady_clock::time_point const start = steady_clock::now();
                quick_hull_.create_initial_simplex(std::cbegin(initial_simplex_),
                                                   std::prev(std::cend(initial_simplex_)));
                auto const delta = duration_cast< microseconds >(steady_clock::now() - start).count();
                log_ << "simplex time = " << delta << "us" << std::endl;
                size_type const basis_size_ = initial_simplex_.size();
                if (basis_size_ != quick_hull_.dimension_ + 1) {
                    err_ << "error: algorithm: cannot construct a simplex. Degenerated input set. Size of basis: "
                         << basis_size_ << std::endl;
                    return false;
                }
            }
            {
                steady_clock::time_point const start = steady_clock::now();
                quick_hull_.create_convex_hull();
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
            return true;
        }

        mutable bool f = true;

        bool
        output(std::ostream & _out) const
        {
            if (operator bool ()) {
                log_ << "dimensionality value " << quick_hull_.dimension_
                     << " is out of supported range: cannot generate output" << std::endl;
                return false;
            }
            if (f) {
                f = false;
                _out << "set view equal xyz\n"
                        "set autoscale\n"
                        "set key left\n"
                        "set xrange [] writeback\n"
                        "set yrange [] writeback\n"
                        "set zrange [] writeback\n";
            } else {
                _out << "pause mouse\n\n"
                        "clear\n"
                        "set xrange restore\n"
                        "set yrange restore\n"
                        "set zrange restore\n";
            }
            _out << "set title \'Points count is " << universe_.size() << "\'\n";
            if (quick_hull_.dimension_ == 2) {
                _out << "plot";
            } else if (quick_hull_.dimension_ == 3) {
                _out << "splot";
            }
            _out << " '-' with points notitle pointtype 4 pointsize 1.5 linetype 1"
                    ", '-' with points notitle"
                    ", '-' with labels offset character 0, character 1 notitle";
            for (auto const & facet_ : quick_hull_.facets_) {
                _out << ", '-' with lines notitle";
                if (!facet_.coplanar_.empty()) {
                    _out << ", '-' with points notitle pointtype 6 pointsize 1.5 linetype 4";
                }
            }
            _out << ";\n";
            {
                for (auto const & v : initial_simplex_) {
                    point const & point_ = *v;
                    for (value_type const & coordinate_ : point_) {
                        _out << coordinate_ << ' ';
                    }
                    _out << '\n';
                }
                _out << "e\n";
                {
                    for (auto const & p : universe_) {
                        point const & point_ = *p;
                        for (value_type const & coordinate_ : point_) {
                            _out << coordinate_ << ' ';
                        }
                        _out << '\n';
                    }
                    _out << "e\n";
                    size_type i = 0;
                    for (auto const & p : universe_) {
                        point const & point_ = *p;
                        for (value_type const & coordinate_ : point_) {
                            _out << coordinate_ << ' ';
                        }
                        _out << i << '\n';
                        ++i;
                    }
                }
                _out << "e\n";
            }
            for (auto const & facet_ : quick_hull_.facets_) {
                auto const & vertices_ = facet_.vertices_;
                for (auto const & vertex_ : vertices_) {
                    for (value_type const & coordinate_ : *vertex_) {
                        _out << coordinate_ << ' ';
                    }
                    _out << '\n';
                }
                for (value_type const & coordinate_ : *vertices_.front()) {
                    _out << coordinate_ << ' ';
                }
                _out << "\n"
                        "e\n";
                if (!facet_.coplanar_.empty()) {
                    for (auto const & v : facet_.coplanar_) {
                        for (value_type const & coordinate_ : *v) {
                            _out << coordinate_ << ' ';
                        }
                        _out << '\n';
                    }
                    _out << "e\n";
                }
            }
            return true;
        }

        friend
        std::ostream &
        operator << (std::ostream & _out, gnuplot const & _gnuplot)
        {
            if (!_gnuplot.output(_out)) {
                throw std::runtime_error("can't generate output due to invalid conditions");
            }
            return _out;
        }

        bool
        slice_layer(bool const _coplanar = false)
        {
            assert(!quick_hull_.facets_.empty());
            std::set< point_iterator_type, point_iterator_less > surface_points_;
            for (auto const & facet_ : quick_hull_.facets_) {
                surface_points_.insert(std::cbegin(facet_.vertices_), std::cend(facet_.vertices_));
                if (_coplanar) {
                    surface_points_.insert(std::cbegin(facet_.coplanar_), std::cend(facet_.coplanar_));
                }
            }
            assert(!surface_points_.empty());
            log_ << "count of points at convex hull is " << surface_points_.size() << std::endl;
            {
                auto ibeg = std::cbegin(universe_);
                auto const iend = std::cend(universe_);
                auto sbeg = std::cbegin(surface_points_);
                point_iterator_less point_iterator_less_;
                assert(sbeg != std::cend(surface_points_));
                while (ibeg != iend) {
                    if (point_iterator_less_(*ibeg, *sbeg)) {
                        ++ibeg;
                    } else {
                        if (!point_iterator_less_(*sbeg, *ibeg)) {
                            universe_.erase(ibeg++);
                        }
                        ++sbeg;
                    }
                }
            }
            initial_simplex_.clear();
            quick_hull_.facets_.clear();
            return (quick_hull_.dimension_ < universe_.size());
        }

        explicit
        operator bool () const
        {
            return (3 < quick_hull_.dimension_);
        }

    };

    gnuplot
    operator () ()
    {
        return {err_, log_, dimension_, points_};
    }

};

int
main(int argc, char * argv[]) // rbox D3 t 100 | quickhull | gnuplot -p
{
    std::ostream & err_ = std::cerr;
    std::ostream & log_ = std::clog;
    err_.sync_with_stdio(false);
    log_.sync_with_stdio(false);

    log_ << "input file: " << ((argc < 2) ? "stdin" : argv[1]) << std::endl;

    std::ifstream ifs_;
    if (argc == 2) {
        ifs_.open(argv[1]);
        if (!ifs_.is_open()) {
            err_ << "error: cannot open file '" << argv[1] << "'" << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::istream & in_ = (ifs_.is_open() ? ifs_ : std::cin);
    in_.sync_with_stdio(false);

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
    test_quickhull< value_type, point, points > test_quickhull_(err_, log_);

    if (!(in_ >> test_quickhull_)) {
        return EXIT_FAILURE;
    }

    auto gnuplot_ = test_quickhull_();

    std::ostream & out_ = std::cout;
    out_.sync_with_stdio(false);
    out_.rdbuf()->pubsetbuf(nullptr, 0);
    if (!!gnuplot_) {
        gnuplot_(true); // rbox D10 t 30 | tee /tmp/q | bin/quickhull >/dev/null ; qconvex Qt Tv s TI /tmp/q
    } else {
        while (gnuplot_(true) && (out_ << gnuplot_) && gnuplot_.slice_layer(true)) { // try `rbox 1000 | bin/quickhull | gnuplot -p`
            continue;
        }
    }
    out_ << std::flush;
    return EXIT_SUCCESS;
}
