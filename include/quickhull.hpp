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
#include <stdexcept>
#include <functional>
#include <tuple>

#include <cassert>

#include <iostream>

struct bad_geometry
        : std::runtime_error
{

    explicit
    bad_geometry(char const * _what)
        : std::runtime_error(_what)
    { ; }

    explicit
    bad_geometry(std::string const & _what)
        : std::runtime_error(_what)
    { ; }

};

template< typename point_type >
struct convex_hull
{

    using size_type = std::size_t;

    using G = typename point_type::value_type;
    using point_refs_type = std::vector< std::reference_wrapper< point_type const > >;
    using point_list = std::list< size_type >;
    using point_set = std::set< size_type >;
    using points_type = std::deque< size_type >;

    size_type dimension_;
    point_refs_type points_;
    point_list internal_set_;

    template< typename ForwardIterator >
    convex_hull(size_type const _dimension,
                ForwardIterator _first, ForwardIterator _last)
        : dimension_(_dimension)
        , points_(_first, _last)
        , matrix_(dimension_ + 1)
        , minor_(dimension_)
    {
        assert(0 < dimension_);
#ifndef NDEBUG
        for (point_type const & point_ : points_) {
            assert(point_.size() == dimension_); // dimensionalities of input points does not match
        }
#endif
        for (size_type i = 0; i < dimension_; ++i) {
            matrix_[i].resize(dimension_ + 1);
            minor_[i].resize(dimension_);
        }
        matrix_[dimension_].resize(dimension_ + 1);
    }

    struct facet;

    using facet_set = std::set< size_type >;

    struct facet // (d - 1)-dimensional faces
    {

        points_type vertices_; // d points : oriented
        bool outward_;         // is top-oriented
        facet_set neighbours_;
        points_type outside_set_; // if not empty, then first point is furthest from this facet
        points_type coplanar_;

        template< typename ForwardIterator >
        facet(ForwardIterator first, ForwardIterator mid, ForwardIterator last,
              bool const _outward)
            : vertices_(first, std::prev(mid))
            , outward_(_outward)
        {
            vertices_.insert(vertices_.cend(), mid, last);
        }

        facet(points_type && _vertices,
              bool const _outward,
              size_type const _neighbour)
            : vertices_(std::move(_vertices))
            , outward_(_outward)
            , neighbours_({_neighbour})
        {
        }

        bool
        order(G const & _nearer, G const & _further) const
        {
            if (outward_) {
                return (_nearer < _further);
            } else {
                return (_further < _nearer);
            }
        }

    };

    using facets_map = std::map< size_type, facet >;
    using facet_iterator = typename facets_map::iterator;
    using facets_type = std::deque< size_type >;
    using vertices_sets_type = std::map< size_type, point_set >;

    facets_map facets_;
    vertices_sets_type ordered_; // ordered, but not oriented vertices of facets

    template< typename ...args >
    facet &
    add_facet(size_type const _facet_key, args &&... _args)
    {
        auto const f = facets_.emplace_hint(facets_.cend(), _facet_key, facet(std::forward< args >(_args)...));
        facet & facet_ = f->second;
        point_set & points_ = ordered_[_facet_key];
        points_.insert(facet_.vertices_.cbegin(), facet_.vertices_.cend());
        return facet_;
    }

    G
    abs(G const & _x) const
    {
        return (_x < G(0.0L)) ? -_x : _x;
    }

    G
    abs(G && _x) const
    {
        return (_x < G(0.0L)) ? -std::move(_x) : std::move(_x);
    }

    bool
    below(facet const & _facet, G const & _orientation) const
    {
        if (_facet.outward_) {
            return (G(0.0L) < _orientation);
        } else {
            return (_orientation < -G(0.0L));
        }
    }

    using row_type = std::valarray< G >;
    using matrix_type = std::vector< row_type >;

    matrix_type matrix_;
    matrix_type minor_;

    G
    det(size_type const _size) // based on LU factorization
    {
        G det_(1.0L);
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
            if (!(G(0.0L) < max_)) { // regular?
                return G(0.0L); // singular
            }
            row_type & ri_ = matrix_[i];
            if (pivot_ != i) {
                det_ = -det_;
                ri_.swap(matrix_[pivot_]);
            }
            G & dia_ = ri_[i];
            G const inv_ = G(1.0L) / dia_;
            det_ *= std::move(dia_);
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

    G
    det()
    {
        return det(dimension_ + 1);
    }

    // http://math.stackexchange.com/questions/822741/
    template< typename vertices >
    G
    orientation(vertices const & _vertices, point_type const & _apex)
    {
        size_type const size_ = _vertices.size(); // dimensionality of the subspace of interest
        assert(!(_apex.size() < size_));
        auto v_ = _vertices.cbegin();
        for (size_type i = 0; i < size_; ++i) {
            assert(v_ != _vertices.cend());
            point_type const & vertex_ = points_.at(*v_);
            ++v_;
            assert(!(vertex_.size() < size_));
            for (size_type j = 0; j < size_; ++j) {
                matrix_[i][j] = vertex_[j];
            }
            matrix_[i][size_] = G(1.0L);
        }
        for (size_type j = 0; j < size_; ++j) {
            matrix_[size_][j] = _apex[j];
        }
        matrix_[size_][size_] = G(1.0L);
        return det(size_ + 1);
    }

    template< typename vertices >
    G
    orientation(vertices const & _vertices, size_type const _apex)
    {
        return orientation(_vertices, points_.at(_apex));
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
        if (!(G(0.0L) < abs(orientation_))) {
            throw bad_geometry("can't find linearly independent point");
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
        if (G(0.0L) < _orientation) {
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
        G orientation_(0.0L);
        while (it != end) {
            auto const next = std::next(it);
            size_type const p = *it;
            G const o_ = orientation(_facet.vertices_, p);
            if (below(_facet, o_)) {
                if (_facet.outside_set_.empty() || _facet.order(orientation_, o_)) {
                    orientation_ = o_;
                    _facet.outside_set_.push_front(p);
                } else {
                    _facet.outside_set_.push_back(p);
                }
                _points.erase(it);
            } else if (!(G(0.0L) < abs(o_))) { // coplanar
                _facet.coplanar_.push_back(p);
            }
            it = next;
        }
        return abs(orientation_);
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
            for (auto second = std::next(first); second != nend; ++second) {
                size_type const s = *second;
                point_set & second_ = ordered_.at(s);
                auto const rend = second_.cend();
                auto r = second_.cbegin();
                auto l = lbeg;
                bool lgood = false;
                bool rgood = false;
                while (l != lend) {
                    if (r == rend) {
                        lgood = (lgood != ((l != lend) && (++l == lend)));
                        break;
                    }
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
                if (lgood) {
                    if (rgood != ((r != rend) && (++r == rend))) {
                        facets_.at(f).neighbours_.insert(s);
                        facets_.at(s).neighbours_.insert(f);
                    }
                }
            }
        }
    }

    void
    create_simplex()
    {
        {
            size_type const size_ = points_.size();
            assert(dimension_ < size_);
            for (size_type i = 0; i < size_; ++i) {
                internal_set_.push_back(i);
            }
        }
        point_list vertices_;
        vertices_.splice(vertices_.end(), internal_set_, internal_set_.begin());
        for (size_type i = 0; i < dimension_; ++i) {
            steal_best(internal_set_, vertices_);
        }
        assert(vertices_.size() == 1 + dimension_); // (d + 1) vertices defining a simplex
        internal_set_.splice(internal_set_.end(), vertices_, vertices_.begin());
        assert(vertices_.size() == dimension_); // d vertices defining a facet
        bool outward_ = !(G(0.0L) < steal_best(internal_set_, vertices_)); // is top oriented?
        auto const vbeg = vertices_.cbegin();
        auto const vend = vertices_.cend();
#ifndef NDEBUG
        point_type inner_point_ = points_.at(*vbeg);
        {
            auto it = vbeg;
            while (++it != vend) {
                inner_point_ += points_.at(*it);
            }
            inner_point_ /= G(1 + dimension_);
        }
#endif
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) { // creation of rest d facets of the simplex
            size_type const n = facets_.size();
            facet & facet_ = add_facet(n, vbeg, exclusive, vend, outward_);
            rank(partition(facet_, internal_set_), n);
            assert(outward_ == !(G(0.0L) < orientation(facet_.vertices_, inner_point_)));
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
    }

    void
    create_convex_hull()
    {
        create_simplex();
        size_type facet_key = facets_.size(); // unique key for facets_
        assert(facet_key == dimension_ + 1);
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
                    if (below(facet_, orientation(facet_.vertices_, apex))) {
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
                points_type const & vertices_ = visible_facet_.vertices_;
                for (size_type const n : visible_facet_.neighbours_) {
                    if (visible_facets_.find(n) == vfend) { // neighbour is not visible
                        point_set const & horizon_ = ordered_.at(n);
                        auto const hend = horizon_.end();
                        points_type ridge_; // horizon ridge + furthest point -> new facet
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
            point_list outside_set_;
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
    }

};
