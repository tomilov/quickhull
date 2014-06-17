#pragma once

#include <vector>
#include <deque>
#include <set>
#include <map>
#include <list>
#include <valarray>
#include <iterator>
#include <algorithm>
#include <utility>
#include <numeric>
#ifdef _DEBUG
#include <ostream>
#include <iostream>
#endif

#include <cassert>

template< typename points_type >
struct convex_hull
{

    using size_type = std::size_t;

    size_type const dimension_;

    using point_type = typename points_type::value_type;
    using G = typename point_type::value_type;

private : // math (simple functions, matrices, etc)

    G
    abs(G const & _x) const
    {
        return (_x < G(0)) ? -_x : _x;
    }

    using row_type = std::valarray< G >;
    using matrix_type = std::vector< row_type >;

    matrix_type matrix_;
    matrix_type shadow_matrix_; // storage for source matrix
    matrix_type minor_;

    void
    transpose()
    { // to cheaper filling columns with ones
        for (size_type row = 0; row < dimension_; ++row) {
            row_type & row_ = shadow_matrix_[row];
            for (size_type col = 1 + row; col < dimension_; ++col) {
                using std::swap;
                swap(shadow_matrix_[col][row], row_[col]);
            }
        }
    }

    void
    restore_matrix()
    {
        for (size_type row = 0; row < dimension_; ++row) {
            matrix_[row] = shadow_matrix_[row];
        }
    }

    void
    restore_matrix(size_type const _identity) // load matrix from storage with replacing of _identity column with ones
    {
        for (size_type row = 0; row < dimension_; ++row) {
            row_type & to_ = matrix_[row];
            if (row == _identity) {
                to_ = G(1);
            } else {
                to_ = shadow_matrix_[row];
            }
        }
    }

    void
    square_matrix(size_type const _size) // matrix_ = shadow_matrix_ * transposed shadow_matrix_
    {
        assert(_size < dimension_);
        for (size_type row = 0; row < _size; ++row) {
            row_type & lhs_ = matrix_[row];
            row_type const & row_ = shadow_matrix_[row];
            auto const rbeg = std::begin(row_);
            auto const rend = std::end(row_);
            for (size_type col = 0; col < _size; ++col) {
                lhs_[col] = std::inner_product(rbeg, rend, std::begin(shadow_matrix_[col]), G(0));
            }
        }
    }

    G
    det(size_type const _size) // based on LU factorization part of LUP decomposition algorithm
    {
        G det_(1);
        for (size_type i = 0; i < _size; ++i) {
            size_type p_ = i;
            G max_ = abs(matrix_[p_][i]);
            size_type pivot_ = p_;
            while (++p_ < _size) {
                G y_ = abs(matrix_[p_][i]);
                if (max_ < y_) {
                    max_ = y_;
                    pivot_ = p_;
                }
            }
            if (!(G(0) < max_)) { // regular?
                return G(0); // singular
            }
            row_type & ri_ = matrix_[i];
            if (pivot_ != i) {
                det_ = -det_; // each permutation flips sign of det
                ri_.swap(matrix_[pivot_]);
            }
            G & dia_ = ri_[i];
            G const inv_ = G(1) / dia_;
            det_ *= std::move(dia_); // det is multiple of diagonal elements
            for (size_type j = 1 + i; j < _size; ++j) {
                matrix_[j][i] *= inv_;
            }
            for (size_type a = 1 + i; a < _size; ++a) {
                row_type & a_ = minor_[a - 1];
                G const & ai_ = matrix_[a][i];
                for (size_type b = 1 + i; b < _size; ++ b) {
                    a_[b - 1] = ai_ * ri_[b];
                }
            }
            for (size_type a = 1 + i; a < _size; ++a) {
                row_type const & a_ = minor_[a - 1];
                row_type & ra_ = matrix_[a];
                for (size_type b = 1 + i; b < _size; ++ b) {
                    ra_[b] -= a_[b - 1];
                }
            }
        }
        return det_;
    }

public :

    using normal_type = std::vector< G >;
    using point_list = std::list< size_type >;
    using point_set = std::set< size_type >;
    using point_array = std::vector< size_type >;
    using point_deque = std::deque< size_type >;

    points_type const & points_;
    point_list internal_set_;

    convex_hull(size_type const _dimension, points_type const & _points)
        : dimension_(_dimension)
        , matrix_(_dimension)
        , shadow_matrix_(_dimension)
        , minor_(_dimension - 1)
        , points_(_points)
    {
        size_type const minor_size_ = dimension_ - 1;
        for (size_type row = 0; row < minor_size_; ++row) {
            matrix_[row].resize(dimension_);
            shadow_matrix_[row].resize(dimension_);
            minor_[row].resize(minor_size_);
        }
        matrix_.back().resize(dimension_);
        shadow_matrix_.back().resize(dimension_);
    }

    struct facet;

    using facet_set = std::set< size_type >;

    struct facet // (d - 1)-dimensional faces
    {

        point_array vertices_; // d points : oriented
        bool outward_;         // is top-oriented
        normal_type normal_;
        G D;
        facet_set neighbours_;
        point_deque outside_set_; // if not empty, then first point is furthest from this facet
        point_deque coplanar_;

        template< typename InputIterator >
        facet(InputIterator first, InputIterator mid, InputIterator last,
              bool const _outward)
            : vertices_(first, std::prev(mid))
            , outward_(_outward)
        {
            vertices_.insert(vertices_.cend(), mid, last);
            normal_.resize(vertices_.size());
        }

        facet(point_array && _vertices,
              bool const _outward,
              size_type const _neighbour)
            : vertices_(std::move(_vertices))
            , outward_(_outward)
            , normal_(vertices_.size())
            , neighbours_({_neighbour})
        { ; }

        G
        distance(point_type const & _point) const
        {
            return std::inner_product(normal_.cbegin(), normal_.cend(), std::begin(_point), D);
        }

        G
        cos_of_dihedral_angle(facet const & _other) const // for faces merging in the future
        {
            return std::inner_product(normal_.cbegin(), normal_.cend(), _other.normal_.cbegin(), G(0));
        }

        bool
        rank(G const & _nearer, G const & _further) const
        {
            if (outward_) {
                return (_nearer < _further);
            } else {
                return (_further < _nearer);
            }
        }

        bool
        above(point_type const & _point) const
        {
            G const distance_ = distance(_point);
            if (outward_) {
                return (G(0) < distance_);
            } else {
                return (distance_ < -G(0));
            }
        }

        bool
        above(G const & _distance) const
        {
            if (outward_) {
                return (G(0) < _distance);
            } else {
                return (_distance < -G(0));
            }
        }

#ifdef _DEBUG
        std::ostream &
        print_info(std::ostream & _out) const
        {
            _out << '[' << (outward_ ? "up" : "dn") << " p: ";
            for (size_type const v : vertices_) {
                _out << v << ';';
            }
            _out << " o: ";
            for (size_type const o : outside_set_) {
                _out << o << ';';
            }
            _out << " c: ";
            for (size_type const c : coplanar_) {
                _out << c << ';';
            }
            _out << " n: ";
            for (size_type const n : neighbours_) {
                _out << n << ';';
            }
            _out << " @: ";
            for (G const & x : normal_) {
                _out << x << ';';
            }
            return _out << ']';
        }
#endif

    };

    using facet_map = std::map< size_type, facet >;

    facet_map facets_;

private : // geometry and basic operation on geometric primitives

    using facet_iterator = typename facet_map::iterator;
    using facets_type = std::vector< size_type >;
    using vertices_sets_type = std::map< size_type, point_set >;

    vertices_sets_type ordered_; // ordered, but not oriented vertices of facets

    void
    set_hyperplane_equation(facet & _facet)
    {
        for (size_type row = 0; row < dimension_; ++row) {
            std::copy_n(std::begin(points_[_facet.vertices_[row]]), dimension_, std::begin(shadow_matrix_[row]));
        }
        transpose();
        G N(0);
        for (size_type i = 0; i < dimension_; ++i) {
            restore_matrix(i);
            G n = det(dimension_);
            N += n * n;
            _facet.normal_[i] = std::move(n);
        }
        using std::sqrt;
        N = -sqrt(N); // minus, else orientation and distance are contradirectional
        restore_matrix();
        _facet.D = -det(dimension_) / N;
        for (size_type i = 0; i < dimension_; ++i) {
            _facet.normal_[i] /= N;
        }
    }

    template< typename ...args >
    facet &
    add_facet(size_type const _facet_key, args &&... _args)
    {
        auto const f = facets_.emplace_hint(facets_.cend(), _facet_key, facet(std::forward< args >(_args)...));
        facet & facet_ = f->second;
        set_hyperplane_equation(facet_);
        ordered_[_facet_key].insert(facet_.vertices_.cbegin(), facet_.vertices_.cend());
        return facet_;
    }

    // http://math.stackexchange.com/questions/822741/
    template< typename vertices >
    G
    orientation(vertices const & _vertices, point_type const & _apex)
    { // rows_-dimensional oriented hypervolume of corresponding parallelotope, strictly positive value for subspaces
        size_type const rows_ = _vertices.size();
        assert(!(dimension_ < rows_));
        row_type & origin_ = shadow_matrix_.back();
        std::copy_n(std::begin(_apex), dimension_, std::begin(origin_));
        auto vertex = _vertices.cbegin();
        if (rows_ == dimension_) {
            for (size_type row = 0; row < dimension_; ++row) {
                assert(vertex != _vertices.cend());
                std::copy_n(std::begin(points_[*vertex]), dimension_, std::begin(matrix_[row]));
                ++vertex;
            }
            for (size_type row = 0; row < dimension_; ++row) { // vectorize
                matrix_[row] -= origin_;
            }
            return det(dimension_);
        } else {
            for (size_type row = 0; row < rows_; ++row) {
                assert(vertex != _vertices.cend());
                std::copy_n(std::begin(points_[*vertex]), dimension_, std::begin(shadow_matrix_[row]));
                ++vertex;
            }
            for (size_type row = 0; row < rows_; ++row) { // vectorize
                shadow_matrix_[row] -= origin_;
            }
            square_matrix(rows_);
            using std::sqrt;
            return sqrt(det(rows_));
        }
    }

    template< typename vertices >
    G
    orientation(vertices const & _vertices, size_type const _apex)
    {
        return orientation(_vertices, points_[_apex]);
    }

    G
    steal_best(point_list & _from, point_list & _to)
    {
        auto it = _from.begin();
        auto const end = _from.end();
        G orientation_ = orientation(_to, *it);
        auto furthest = it;
        while (++it != end) {
            G const o_ = orientation(_to, *it);
            if (abs(orientation_) < abs(o_)) {
                orientation_ = o_;
                furthest = it;
            }
        }
        _to.splice(_to.cend(), _from, furthest);
        return orientation_;
    }

    using ranking_type = std::multimap< G, size_type >;
    using ranking_meta_type = std::map< size_type, typename ranking_type::iterator >;

    ranking_type ranking_;
    ranking_meta_type ranking_meta_;

    void
    rank(G const _orientation, size_type const _facet)
    {
        if (G(0) < _orientation) {
            auto const r = ranking_.emplace(_orientation, _facet);
            ranking_meta_.emplace(_facet, r);
        }
    }

    void
    unrank(size_type const _facet)
    {
        auto const r = ranking_meta_.find(_facet);
        if (r != ranking_meta_.end()) {
            ranking_.erase(r->second);
            ranking_meta_.erase(r);
        }
    }

    size_type
    get_furthest(size_type const _bad_value) const
    {
        if (ranking_.empty()) {
            return _bad_value;
        } else {
            auto const r = std::prev(ranking_.cend());
            return r->second;
        }
    }

    G
    partition(facet & _facet, point_list & _points)
    {
        auto it = _points.begin();
        auto const end = _points.end();
        G distance_(0);
        while (it != end) {
            auto const next = std::next(it);
            size_type const p = *it;
            G const d_ = _facet.distance(points_[p]);
            if (_facet.above(d_)) {
                if (_facet.outside_set_.empty() || _facet.rank(distance_, d_)) {
                    distance_ = d_;
                    _facet.outside_set_.push_front(p);
                } else {
                    _facet.outside_set_.push_back(p);
                }
                _points.erase(it);
            } else if (!(G(0) < abs(d_))) { // coplanar
                _facet.coplanar_.push_back(p);
            }
            it = next;
        }
        return abs(distance_);
    }

    void
    adjacency(facets_type const & _newfacets)
    {
        auto const nend = _newfacets.cend();
        for (auto first = _newfacets.cbegin(); first != nend; ++first) {
            size_type const f = *first;
            point_set & first_ = ordered_.at(f);
            auto const lbeg = first_.cbegin();
            auto const lend = first_.cend();
            facet & first_facet_ = facets_.at(f);
            for (auto second = std::next(first); second != nend; ++second) {
                size_type const s = *second;
                point_set & second_ = ordered_.at(s);
                auto const rend = second_.cend();
                auto r = second_.cbegin();
                auto l = lbeg;
                bool lgood = false;
                bool rgood = false;
                while ((l != lend) && (r != rend)) {
                    size_type const left = *l;
                    size_type const right = *r;
                    if (left < right) {
                        if (lgood) {
                            lgood = false;
                            break;
                        } else {
                            lgood = true;
                        }
                        ++l;
                    } else {
                        if (right < left) {
                            if (rgood) {
                                lgood = false;
                                break;
                            } else {
                                rgood = true;
                            }
                        } else {
                            ++l;
                        }
                        ++r;
                    }
                }
                if (lgood != ((l != lend) && (++l == lend))) {
                    if (rgood != ((r != rend) && (++r == rend))) {
                        first_facet_.neighbours_.insert(s);
                        facets_.at(s).neighbours_.insert(f);
                    }
                }
            }
        }
    }

public : // largest possible simplex heuristic, convex hull algorithm

    bool
    create_simplex()
    {
        {
            size_type const size_ = points_.size();
            for (size_type i = 0; i < size_; ++i) {
                internal_set_.push_back(i);
            }
        }
        point_list vertices_;
        vertices_.splice(vertices_.cend(), internal_set_, internal_set_.begin());
        for (size_type i = 0; i < dimension_; ++i) {
            G const orientation_ = steal_best(internal_set_, vertices_);
            if (!(G(0) < abs(orientation_))) {
                return false; // can't find linearly independent point
            }
        }
        assert(vertices_.size() == dimension_ + 1); // (d + 1) vertices defining a d-simplex
        internal_set_.splice(internal_set_.cend(), vertices_, vertices_.begin());
        assert(vertices_.size() == dimension_); // d vertices defining a facet
        bool outward_ = !(G(0) < steal_best(internal_set_, vertices_)); // is top oriented?
        auto const vbeg = vertices_.cbegin();
        auto const vend = vertices_.cend();
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) { // creation of rest d facets of the simplex
            size_type const n = facets_.size();
            facet & facet_ = add_facet(n, vbeg, exclusive, vend, outward_);
            rank(partition(facet_, internal_set_), n);
            outward_ = !outward_;
        }
        assert(dimension_ + 1 == facets_.size()); // simplex
        { // adjacency
            auto const fbeg = facets_.begin();
            auto const fend = facets_.end();
            for (auto i = fbeg; i != fend; ++i) {
                facet_set & neighbours_ = i->second.neighbours_;
                for (auto j = fbeg; j != fend; ++j) {
                    if (j != i) {
                        neighbours_.insert(j->first);
                    }
                }
            }
        }
#ifdef _DEBUG
        show_scene(std::cerr);
        print_info(std::cerr);
        pause(std::cerr, "simplex");
#endif
        return true;
    }

    bool
    create_convex_hull()
    {/*
        if (!create_simplex()) {
            return false;
        }*/
        size_type facet_key = facets_.size(); // unique key for facets_
        assert(facet_key == dimension_ + 1);
        point_list outside_set_;
        facet_set visited_;
        facet_set viewable_;
        facet_set visible_facets_;
        facets_type newfacets_;
        for (size_type furthest = get_furthest(facet_key); furthest != facet_key; furthest = get_furthest(facet_key)) {
            facet & facet_ = facets_.at(furthest);
            size_type const apex = facet_.outside_set_.front();
            point_type const & apex_ = points_[apex];
            visible_facets_ = {furthest};
            { // find visible facets
                visited_ = {furthest};
                viewable_ = facet_.neighbours_;
                while (!viewable_.empty()) {
                    auto const first = viewable_.begin();
                    size_type const f = *first;
                    facet const & watchable_ = facets_.at(f);
                    if (watchable_.above(apex_)) {
                        visible_facets_.insert(f);
                        std::set_difference(watchable_.neighbours_.cbegin(), watchable_.neighbours_.cend(),
                                            visited_.cbegin(), visited_.cend(),
                                            std::inserter(viewable_, viewable_.end()));
                    }
                    visited_.insert(f);
                    viewable_.erase(first);
                }
            }
            // the boundary of visible facets is the set of horizon ridges
            // Each ridge signifies the adjacency of two facets.
            auto const vfend = visible_facets_.end();
            for (size_type const v : visible_facets_) {
                facet const & visible_facet_ = facets_.at(v);
                point_array const & vertices_ = visible_facet_.vertices_;
                for (size_type const n : visible_facet_.neighbours_) {
                    if (visible_facets_.find(n) == vfend) { // neighbour is not visible
                        point_set const & horizon_ = ordered_.at(n);
                        auto const hend = horizon_.end();
                        point_array ridge_; // horizon ridge + furthest point -> new facet
                        ridge_.reserve(dimension_);
                        for (size_type const p : vertices_) { // facets intersection with keeping of points order as in visible facet
                            auto const h = horizon_.find(p);
                            if (h == hend) {
                                ridge_.push_back(apex); // insert furthest point instead of inner point of visible facet
                            } else {
                                ridge_.push_back(p);
                            }
                        }
                        assert(ridge_.size() == dimension_); // ridge_ contains newfacet vertices (ridge + current furthest point)
                        { // replace visible facet became internal with newly created facet in neighbours set
                            facet & horizon_facet_ = facets_.at(n);
                            horizon_facet_.neighbours_.erase(v);
                            horizon_facet_.neighbours_.insert(horizon_facet_.neighbours_.cend(), facet_key);
                        }
                        newfacets_.push_back(facet_key);
                        add_facet(facet_key, std::move(ridge_), visible_facet_.outward_, n);
                        ++facet_key;
                    }
                }
            }
            adjacency(newfacets_);
#ifdef _DEBUG
            show_scene(std::cerr, newfacets_, visible_facets_, furthest);
            print_info(std::cerr);
            pause(std::cerr, "new facets, visible facets and furthest point");
#endif
            facet_.outside_set_.pop_front(); // already assorted
            for (size_type const v : visible_facets_) { // remove visible facets and gather outside points from them
                auto const visible_facet = facets_.find(v);
                assert(visible_facet != facets_.end());
                facet const & visible_facet_ = visible_facet->second;
                outside_set_.insert(outside_set_.cend(),
                                    visible_facet_.outside_set_.cbegin(),
                                    visible_facet_.outside_set_.cend());
                facets_.erase(visible_facet);
                ordered_.erase(v);
                unrank(v);
            }
            for (size_type const n : newfacets_) {
                rank(partition(facets_.at(n), outside_set_), n);
            }
            newfacets_.clear();
            internal_set_.splice(internal_set_.cend(), outside_set_);
#ifdef _DEBUG
            show_scene(std::cerr);
            print_info(std::cerr);
            pause(std::cerr, "quick hull step end");
#endif
        }
        ordered_.clear();
#ifdef _DEBUG
        std::cerr << "print 'convex hull'\n";
#endif
        return true;
    }

#ifdef _DEBUG
    std::ostream &
    show_scene(std::ostream & _out, facets_type const & _newfacets, facet_set const & _visible_facets, size_type const _furthest) const
    {
        _out << "clear\n";
        if (!facets_.empty()) {
            _out << "splot ";
            //size_type const size_ = facets_.size();
            auto const fbeg = facets_.cbegin();
            auto const fend = facets_.cend();
            auto const nend = _newfacets.cend();
            auto const vend = _visible_facets.cend();
            for (auto f_ = fbeg; f_ != fend; ++f_) {
                size_type const f = f_->first;
                _out << "'-' with lines notitle linetype 1 ";
                if (std::find(_newfacets.cbegin(), nend, f) != nend) {
                    _out << "linecolor rgb 'red'";
                } else if (_visible_facets.find(f) != vend) {
                    _out << "linecolor rgb 'blue'";
                } else if (_furthest == f) {
                    _out << "linecolor rgb 'green'";
                } else {
                    _out << "linecolor rgb 'black'";
                }
                //facet const & facet_ = f_->second;
                _out << ", ";
            }
            _out << "'-' with points notitle pointtype 6 pointsize 1, ";
            _out << "'-' with points notitle pointtype 1 linecolor rgb 'purple', ";
            _out << "'-' with labels notitle offset character 0, character 0.5 font 'Times,6'";
            _out << ";\n";
            for (auto const & f : facets_) {
                point_array const & vertices_ = f.second.vertices_;
                for (size_type const v : vertices_) {
                    point_type const & point_ = points_[v];
                    std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                    _out << '\n';
                }
                point_type const & point_ = points_[vertices_.front()];
                std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                _out << "\ne\n";
            }
            {
                point_type const & apex_ = points_[facets_.at(_furthest).outside_set_.front()];
                std::copy(std::begin(apex_), std::end(apex_), std::ostream_iterator< G >(_out, " "));
                _out << "\ne\n";
            }
            for (point_type const & point_ : points_) {
                std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                _out << '\n';
            }
            _out << "e\n";
            size_type const points_count = points_.size();
            for (size_type p = 0; p < points_count; ++p) {
                point_type const & point_ = points_[p];
                std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                _out << p << '\n';
            }
            _out << "e\n";
        }
        return _out;
    }

    std::ostream &
    show_scene(std::ostream & _out) const
    {
        _out << "clear\n";
        if (!facets_.empty()) {
            _out << "splot ";
            //size_type const size_ = facets_.size();
            auto const fbeg = facets_.cbegin();
            auto const fend = facets_.cend();
            for (auto f_ = fbeg; f_ != fend; ++f_) {
                _out << "'-' with lines notitle linetype 1 linecolor rgb 'black', ";
            }
            _out << "'-' with points notitle pointtype 1 linecolor rgb 'purple', ";
            _out << "'-' with labels notitle offset character 0, character 0.5 font 'Times,6', ";
            _out << "'-' with vectors head filled linetype 1 linecolor rgb 'green', ";
            _out << "'-' with points pointtype 7 linecolor rgb 'red'";
            _out << ";\n";
            for (auto const & f : facets_) {
                point_array const & vertices_ = f.second.vertices_;
                for (size_type const v : vertices_) {
                    point_type const & point_ = points_[v];
                    std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                    _out << '\n';
                }
                point_type const & point_ = points_[vertices_.front()];
                std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                _out << "\ne\n";
            }
            for (point_type const & point_ : points_) {
                std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                _out << '\n';
            }
            _out << "e\n";
            size_type const points_count = points_.size();
            for (size_type p = 0; p < points_count; ++p) {
                point_type const & point_ = points_[p];
                std::copy(std::begin(point_), std::end(point_), std::ostream_iterator< G >(_out, " "));
                _out << p << '\n';
            }
            _out << "e\n";
            std::vector< point_type > centres_(facets_.size());
            auto const cbeg = centres_.begin();
            auto const cend = centres_.end();
            auto cit = cbeg;
            for (auto const & f : facets_) {
                facet const & facet_ = f.second;
                assert(cit != cend);
                point_type & centre_ = *cit++;
                centre_.resize(dimension_, G(0));
                for (size_type const v : facet_.vertices_) {
                    point_type const & vertex_ = points_[v];
                    for (size_type i = 0; i < dimension_; ++i) {
                        centre_[i] += vertex_[i];
                    }
                }
                G const norm_ = G(1) / G(dimension_);
                for (G & coordinate_ : centre_) {
                    coordinate_ *= norm_;
                    _out << coordinate_ << ' ';
                }
                for (G const & component_ : facet_.normal_) {
                    _out << (facet_.outward_ ? component_ : -component_) << ' ';
                }
                _out << '\n';
            }
            _out << "e\n";
            for (point_type const & centre_ : centres_) {
                std::copy_n(std::begin(centre_), dimension_, std::ostream_iterator< G >(_out, " "));
                _out << '\n';
            }
            _out << "e\n";
        }
        return _out;
    }

    std::ostream &
    print_info(std::ostream & _out) const
    {
        _out << "print 'facets (" << facets_.size() << "): ";
        for (auto const & f : facets_) {
            _out << f.first << ';';
        }
        _out << "'\n";
        _out << "print '{'\n";
        for (auto const & f : facets_) {
            _out << "print '\tfacet #" << f.first << " = ";
            f.second.print_info(_out) << "'\n";
        }
        _out << "print '}'\n";
        return _out;
    }

    std::ostream &
    pause(std::ostream & _out, std::string const & _text = "") const
    {
        return _out << "pause -1 '" << _text << "'\n";
    }
#endif

};
