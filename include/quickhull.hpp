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

#include <cassert>

template< typename points_type >
struct convex_hull
{

    using size_type = std::size_t;

    using point_type = typename points_type::value_type;
    using G = typename point_type::value_type;

private : // math (simple functions, matrices, etc)

    G
    abs(G const & _x) const
    {
        return (_x < G(0)) ? -_x : _x;
    }

    G
    sqrtsign(G const & _x) const
    {
        using std::sqrt;
        return (_x < G(0)) ? -sqrt(-_x) : sqrt(_x);
    }

    using row_type = std::valarray< G >;
    using matrix_type = std::vector< row_type >;

    matrix_type matrix_;
    matrix_type shadow_matrix_; // storage for source matrix
    matrix_type minor_;

    void
    vectorize(size_type const _rows)
    {
        row_type const & origin_ = shadow_matrix_.back();
        for (size_type row = 0; row < _rows; ++row) {
            shadow_matrix_[row] -= origin_;
        }
    }

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
            std::copy_n(std::begin(shadow_matrix_[row]), dimension_, std::begin(matrix_[row]));
        }
    }

    void
    restore_matrix(size_type const _identity)
    {
        for (size_type row = 0; row < dimension_; ++row) {
            auto const to = std::begin(matrix_[row]);
            if (row == _identity) {
                std::fill_n(to, dimension_, G(1));
            } else {
                std::copy_n(std::begin(shadow_matrix_[row]), dimension_, to);
            }
        }
    }

    void
    square_matrix(size_type const _size)
    {
        assert(_size < dimension_);
        for (size_type row = 0; row < _size; ++row) {
            row_type & lhs_ = matrix_[row];
            row_type const & row_ = shadow_matrix_[row];
            auto const rbeg = std::begin(row_);
            auto const rend = std::end(row_);
            for (size_type col = 0; col < _size; ++col) {
                row_type const & col_ = shadow_matrix_[col];
                lhs_[col] = std::inner_product(rbeg, rend, std::begin(col_), G(0));
            }
        }
    }

    G
    det(size_type const _size) // based on LU factorization
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
                det_ = -det_;
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

    size_type const dimension_;
    points_type const & points_;
    point_list internal_set_;

    convex_hull(size_type const _dimension, points_type const & _points)
        : dimension_(_dimension)
        , points_(_points)
        , matrix_(_dimension + 1)
        , shadow_matrix_(_dimension + 1)
        , minor_(_dimension)
    {
        for (size_type row = 0; row < dimension_; ++row) {
            matrix_[row].resize(dimension_ + 1);
            shadow_matrix_[row].resize(dimension_ + 1);
            minor_[row].resize(dimension_);
        }
        matrix_.back().resize(dimension_ + 1);
        shadow_matrix_.back().resize(dimension_ + 1);
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
        cos_dihedral_angle(facet const & _other) const // for faces merging in the future
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
        above(G const & _distance) const
        {
            if (outward_) {
                return (G(0) < _distance);
            } else {
                return (_distance < -G(0));
            }
        }

    };

    using facets_map = std::map< size_type, facet >;

    facets_map facets_;

private : // geometry and basic operation on geometric primitives

    using facet_iterator = typename facets_map::iterator;
    using facets_type = std::deque< size_type >;
    using vertices_sets_type = std::map< size_type, point_set >;

    vertices_sets_type ordered_; // ordered, but not oriented vertices of facets

    void
    set_hyperplane_equation(facet & _facet)
    {
        for (size_type row = 0; row < dimension_; ++row) {
            point_type const & vertex_ = points_[_facet.vertices_[row]];
            std::copy_n(std::begin(vertex_), dimension_, std::begin(shadow_matrix_[row]));
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
        point_set & sorted_vertices_ = ordered_[_facet_key];
        sorted_vertices_.insert(facet_.vertices_.cbegin(), facet_.vertices_.cend());
        return facet_;
    }

    bool
    above(facet const & _facet, size_type const _point) const
    {
        return _facet.above(_facet.distance(points_[_point]));
    }

    // http://math.stackexchange.com/questions/822741/
    template< typename vertices >
    G
    orientation(vertices const & _vertices, point_type const & _apex)
    {
        size_type const rows_ = _vertices.size();
        assert(!(dimension_ < rows_));
        auto vertex = _vertices.cbegin();
        if (rows_ == dimension_) { // dimension_-dimensional oriented hypervolume
            for (size_type row = 0; row < dimension_; ++row) {
                assert(vertex != _vertices.cend());
                row_type & row_ = matrix_[row];
                std::copy_n(std::begin(points_[*vertex]), dimension_, std::begin(row_));
                row_[dimension_] = G(1);
                ++vertex;
            }
            auto const origin = std::begin(matrix_.back());
            std::copy_n(std::begin(_apex), dimension_, origin);
            origin[dimension_] = G(1);
            return det(dimension_ + 1);
        } else { // rows_-dimensional hypervolumes for rows_-dimensional subspaces
            for (size_type row = 0; row < rows_; ++row) {
                assert(vertex != _vertices.cend());
                row_type & row_ = shadow_matrix_[row];
                std::copy_n(std::begin(points_[*vertex]), dimension_, std::begin(row_));
                ++vertex;
            }
            auto const origin = std::begin(shadow_matrix_.back());
            std::copy_n(std::begin(_apex), dimension_, origin);
            vectorize(rows_);
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
        _to.splice(_to.end(), _from, furthest);
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
        auto const nend = _newfacets.end();
        for (auto first = _newfacets.begin(); first != nend; ++first) {
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

public : // quick hull algorithm

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
        vertices_.splice(vertices_.end(), internal_set_, internal_set_.begin());
        for (size_type i = 1; i < dimension_; ++i) {
            G const orientation_ = steal_best(internal_set_, vertices_);
            if (!(G(0) < orientation_)) {
                return false; // can't find linearly independent point
            }
        }
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
        for (size_type furthest = get_furthest(facet_key); furthest != facet_key; furthest = get_furthest(facet_key)) {
            facet & facet_ = facets_.at(furthest);
            size_type const apex = facet_.outside_set_.front();
            facet_.outside_set_.pop_front();
            facet_set visible_facets_{furthest};
            { // find visible facets
                facet_set pool_ = facet_.neighbours_;
                facet_set visited_{furthest};
                while (!pool_.empty()) {
                    auto const first = pool_.begin();
                    size_type const f = *first;
                    facet const & facet_ = facets_.at(f);
                    if (above(facet_, apex)) {
                        visible_facets_.insert(f);
                        std::set_difference(facet_.neighbours_.cbegin(), facet_.neighbours_.cend(),
                                            visited_.cbegin(), visited_.cend(),
                                            std::inserter(pool_, pool_.end()));
                    }
                    visited_.insert(f);
                    pool_.erase(first);
                }
            }
            // the boundary of visible facets is the set of horizon ridges
            // Each ridge signifies the adjacency of two facets.
            facets_type newfacets_;
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
                        { // replace visible facet became internal with newly created facet in adjacency
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
            internal_set_.splice(internal_set_.end(), outside_set_);
        }
        ordered_.clear();
        return true;
    }

};
