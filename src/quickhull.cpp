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
#define TERM_COLOR(code) __extension__ "\e[1;" #code "m"
#define TERM_COLOR_RED TERM_COLOR(31)
#define TERM_COLOR_GREEN TERM_COLOR(32)
#define TERM_COLOR_BLUE TERM_COLOR(34)
#define TERM_COLOR_DEFAULT __extension__ "\e[0m"
#else
#define TERM_COLOR_RED ""
#define TERM_COLOR_GREEN ""
#define TERM_COLOR_BLUE ""
#define TERM_COLOR_DEFAULT ""
#endif

template< typename value_type = float,
          typename point = std::valarray< value_type >,
          typename points = std::vector< point > >
struct qh
{

    using size_type = std::size_t;

    std::ostream & err_;
    std::ostream & log_;

    qh(std::ostream & _err = std::cerr,
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
    input(std::istream & _rbox)
    {
        {
            if (!std::getline(_rbox, line_)) {
                err_ << "error: io: missing dimension line" << std::endl;
                return false;
            }
            iss_.str(line_);
            if (!(iss_ >> dimension_)) {
                err_ << "error: io: dimension format" << std::endl;
                return false;
            }
            log_ << "dimensionality of input is " << dimension_ << std::endl;
            if (!(1 < dimension_)) {
                err_ << "error: io: dimensionality value is not greater then one" << std::endl;
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
            if (!std::getline(_rbox, line_)) {
                err_ << "error: io: missing count line" << std::endl;
                return false;
            }
            iss_.str(line_);
            if (!(iss_ >> count_)) {
                err_ << "error: io: format of count" << std::endl;
                return false;
            }
            iss_.clear();
            log_ << "input points count = " << count_ << std::endl;
            if (!(dimension_ < count_)) {
                err_ << "error: io: points count is less than or equal to dimensionality" << std::endl;
                return false;
            }
        }
        points_ = points(count_);
        for (point & point_ : points_) {
            if (!std::getline(_rbox, line_)) {
                err_ << "error: io: wrong line count" << std::endl;
                return false;
            }
            point_.resize(dimension_);
            {
                iss_.str(line_);
                auto c = std::begin(point_);
                for (size_type j = 0; j < dimension_; ++j) {
                    if (!(iss_ >> *c)) {
                        err_ << "error: io: bad corodinate value at line " << j << " of data" << std::endl;
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
    operator >> (std::istream & _rbox, qh & _qh)
    {
        _qh.input(_rbox);
        return _rbox;
    }

    value_type const eps = std::numeric_limits< value_type >::epsilon();

    using point_iterator_type = typename points::iterator;

    struct point_iterator_less
    {

        bool
        operator () (point_iterator_type const & _lhs, point_iterator_type const & _rhs) const
        {
            return (std::addressof(*_lhs) < std::addressof(*_rhs));
        }

    };

    using quick_hull_type = quick_hull< point_iterator_type >;

    using internal_set = typename quick_hull_type::point_list;
    using initial_simplex = typename quick_hull_type::point_array;
    using facets = typename quick_hull_type::facets;

    internal_set internal_set_;
    initial_simplex initial_simplex_;
    facets facets_;

    mutable bool f = false;

    void
    init()
    {
        auto p = std::begin(points_);
        auto const pend = std::end(points_);
        while (p != pend) {
            internal_set_.push_back(p);
            ++p;
        }
        internal_set_.sort(point_iterator_less{});
        f = false;
    }

    bool
    operator () (bool const _use_simplex_heuristic = false)
    {
        assert(dimension_ < internal_set_.size());
        assert(initial_simplex_.empty());
        assert(facets_.empty());

        log_ << "\ncount of points to process is " << internal_set_.size() << std::endl;

        quick_hull_type quick_hull_(dimension_, eps);
        {
            auto p = std::cbegin(internal_set_);
            if (!_use_simplex_heuristic) {
                std::copy_n(p, (dimension_ + 1), std::back_inserter(initial_simplex_));
            }
            std::copy(p, std::cend(internal_set_), std::back_inserter(quick_hull_.internal_set_));
        }

        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        {
            steady_clock::time_point const start = steady_clock::now();
            if (_use_simplex_heuristic) {
                initial_simplex_ = quick_hull_.create_initial_simplex();
            } else {
                quick_hull_.create_initial_simplex(std::cbegin(initial_simplex_),
                                                   std::prev(std::cend(initial_simplex_)));
            }
            log_ << "simplex time = " << duration_cast< microseconds >(steady_clock::now() - start).count()
                 << "us" << std::endl;
            size_type const basis_size_ = initial_simplex_.size();
            if (basis_size_ != dimension_ + 1) {
                err_ << "error: algorithm: cannot create a simplex. Degenerated input set. Size of basis: "
                     << basis_size_ << std::endl;
                return false;
            }
        }
        {
            steady_clock::time_point const start = steady_clock::now();
            quick_hull_.create_convex_hull();
            log_ << "quickhull time = "
                 << TERM_COLOR_GREEN << duration_cast< microseconds >(steady_clock::now() - start).count()
                 << "us" << TERM_COLOR_DEFAULT << std::endl;
            if (!quick_hull_.check()) {
                err_ << TERM_COLOR_RED << "error: algorithm: resulting structure is not valid convex polytope"
                     << TERM_COLOR_DEFAULT << std::endl;
                return false;
            }
        }
        log_ << "number of (convex hull) polyhedron facets is "
             << TERM_COLOR_BLUE << quick_hull_.facets_.size()
             << TERM_COLOR_DEFAULT << std::endl;
        facets_ = std::move(quick_hull_.facets_);
        return true;
    }

    bool
    output(std::ostream & _gnuplot) const
    {
        if (3 < dimension_) {
            log_ << "dimensionality value (" << dimension_
                 << ") is out of supported range: cannot generate output" << std::endl;
            return false;
        }
        _gnuplot << "clear\n";
        _gnuplot << "set view equal xyz\n";
        _gnuplot << "set autoscale\n";
        _gnuplot << "set key left\n";
        if (f) {
            _gnuplot << "set xrange restore; set yrange restore; set zrange restore\n";
        } else {
            f = true;
            _gnuplot << "set xrange [] writeback; set yrange [] writeback; set zrange [] writeback;\n";
        }
        _gnuplot << "set title \'Points count is " << internal_set_.size() << "\'\n";
        if (dimension_ == 2) {
            _gnuplot << "plot";
        } else if (dimension_ == 3) {
            _gnuplot << "splot";
        } else {
            assert(false);
        }
        _gnuplot << " '-' with points notitle pointtype 4 pointsize 1.5 linetype 1"
                    ", '-' with points notitle"
                    ", '-' with labels offset character 0, character 1 notitle";
        for (auto i = std::distance(std::cbegin(facets_), std::cend(facets_)); 0 < i; --i) {
            _gnuplot << ", '-' with lines notitle"
                        ", '-' with points notitle pointtype 6 pointsize 1.5 linetype 4";
        }
        _gnuplot << ";\n";
        {
            for (auto const v : initial_simplex_) {
                point const & point_ = *v;
                for (value_type const & coordinate_ : point_) {
                    _gnuplot << coordinate_ << ' ';
                }
                _gnuplot << '\n';
            }
            _gnuplot << "e\n";
            {
                for (auto const & p : internal_set_) {
                    point const & point_ = *p;
                    for (value_type const & coordinate_ : point_) {
                        _gnuplot << coordinate_ << ' ';
                    }
                    _gnuplot << '\n';
                }
                _gnuplot << "e\n";
                size_type i = 0;
                for (auto const & p : internal_set_) {
                    point const & point_ = *p;
                    for (value_type const & coordinate_ : point_) {
                        _gnuplot << coordinate_ << ' ';
                    }
                    _gnuplot << i << '\n';
                    ++i;
                }
            }
            _gnuplot << "e\n";
        }
        for (auto const & facet_ : facets_) {
            auto const & vertices_ = facet_.vertices_;
            for (auto const vertex_ : vertices_) {
                for (value_type const & coordinate_ : *vertex_) {
                    _gnuplot << coordinate_ << ' ';
                }
                _gnuplot << '\n';
            }
            point const & first_vertex_ = *vertices_.front();
            for (value_type const & coordinate_ : first_vertex_) {
                _gnuplot << coordinate_ << ' ';
            }
            _gnuplot << "\n";
            _gnuplot << "e\n";
            for (auto const v : facet_.coplanar_) {
                for (value_type const & coordinate_ : *v) {
                    _gnuplot << coordinate_ << ' ';
                }
                _gnuplot << '\n';
            }
            _gnuplot << "e\n";
        }
        _gnuplot << "pause mouse\n";
        return true;
    }

    friend
    std::ostream &
    operator << (std::ostream & _gnuplot, qh const & _qh)
    {
        if (!_qh.output(_gnuplot)) {
            // exception?
        }
        return _gnuplot;
    }

    bool
    slice_layer()
    {
        assert(!facets_.empty());
        std::set< point_iterator_type, point_iterator_less > surface_points_;
        for (auto const & facet_ : facets_) {
            surface_points_.insert(std::cbegin(facet_.vertices_), std::cend(facet_.vertices_));
        }
        assert(!surface_points_.empty());
        log_ << "count of points at convex hull is " << surface_points_.size() << std::endl;
        {
            auto ibeg = std::cbegin(internal_set_);
            auto const iend = std::cend(internal_set_);
            auto sbeg = std::cbegin(surface_points_);
            point_iterator_less point_iterator_less_;
            assert(sbeg != std::cend(surface_points_));
            while (ibeg != iend) {
                if (point_iterator_less_(*ibeg, *sbeg)) {
                    ++ibeg;
                } else {
                    if (!point_iterator_less_(*sbeg, *ibeg)) {
                        internal_set_.erase(ibeg++);
                    }
                    ++sbeg;
                }
            }
        }
        initial_simplex_.clear();
        facets_.clear();
        log_ << std::endl;
        return (dimension_ < internal_set_.size());
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

    std::istream & rbox_ = (ifs_.is_open() ? ifs_ : std::cin);
    rbox_.sync_with_stdio(false);

#if 0
    using value_type = float;
    using point = std::vector< value_type >;
    using points = std::vector< point >;
#else
    using value_type = double;
    using point = std::forward_list< value_type >;
    using points = std::forward_list< point >;
#endif
    qh< value_type, point, points > qh_(err_, log_);

    if (!(rbox_ >> qh_)) {
        return EXIT_FAILURE;
    }

    std::ostream & gnuplot_ = std::cout;
    gnuplot_.sync_with_stdio(false);
    //out_.rdbuf()->pubsetbuf(nullptr, 0);
    qh_.init();
    while (qh_(true)) {
        if (!qh_.output(gnuplot_)) {
            break;
        }
        if (!qh_.slice_layer()) {
            break;
        }
    }
    gnuplot_ << std::flush;
    // test: rbox D10 t 30 | tee /tmp/q | bin/quickhull >/dev/null ; qconvex Qt Tv s TI /tmp/q
    return EXIT_SUCCESS;
}
