#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <deque>
#include <set>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <utility>
#include <numeric>
#include <exception>

#include <cassert>

struct bad_geometry
        : std::exception
{

    virtual
    ~bad_geometry() noexcept = default;

    bad_geometry() = default;

    explicit
    bad_geometry(const char * const _what)
        : what_(_what)
    { ; }

    explicit
    bad_geometry(std::string const & _what)
        : what_(_what)
    { ; }

    virtual
    const char *
    what() const noexcept override
    {
        return what_.c_str();
    }

private :

    std::string const what_ = "bad_geometry";

};

template< typename F >
F
determinant(boost::numeric::ublas::matrix< F > _m)
{
    using size_type = typename boost::numeric::ublas::matrix< F >::size_type;
    size_type const size_ = _m.size1();
    assert(_m.size2() == size_);
    boost::numeric::ublas::permutation_matrix< size_type > pm_(size_);
    if (0 < boost::numeric::ublas::lu_factorize(_m, pm_)) {
        return F(0.0L); // singular matrix
    } else {
        F determinant_(1.0L);
        for (size_type i = 0; i < size_; ++i) {
            if (i == pm_(i)) {
                determinant_ *= +_m(i, i);
            } else {
                determinant_ *= -_m(i, i);
            }
        }
        return determinant_;
    }
}

template< typename point_type >
struct convex_hull
{

    using size_type = std::size_t;

    using G = typename point_type::value_type;
    using point_refs_type = std::deque< std::reference_wrapper< point_type const > >;
    using point_list = std::list< size_type >;
    using point_set = std::set< size_type >;
    using points_type = std::deque< size_type >;

    template< typename ForwardIterator >
    convex_hull(ForwardIterator _first, ForwardIterator _last)
        : dimension_(_first->size())
        , points_(_first, _last)
    {
        assert(0 < dimension_);
        for (point_type const & point_ : points_) {
            if (point_.size() != dimension_) {
                throw bad_geometry("dimensionalities does not match");
            }
        }
    }

    size_type dimension_;
    point_refs_type points_;
    point_list internal_set_;

    struct facet;

    using facet_set = std::set< size_type >;

    struct facet // (d - 1)-dimensional faces
    {

        template< typename ForwardIterator >
        facet(ForwardIterator first, ForwardIterator mid, ForwardIterator last,
              bool const _outward)
            : vertices_(first, std::prev(mid))
            , points_(vertices_.cbegin(), vertices_.cend())
            , outward_(_outward)
        {
            auto const rest = vertices_.insert(vertices_.cend(), mid, last);
            points_.insert(points_.cend(), rest, vertices_.end());
            std::sort(points_.begin(), points_.end());
        }

        facet(points_type && _vertices,
              bool const _outward,
              size_type const _neighbour)
            : vertices_(std::move(_vertices))
            , points_(vertices_.cbegin(), vertices_.cend())
            , outward_(_outward)
            , neighbours_({_neighbour})
        {
            std::sort(points_.begin(), points_.end());
        }

        points_type vertices_; // oriented
        points_type points_;   // non-oriented
        bool outward_;
        facet_set neighbours_;
        points_type outside_set_; // if not empty, then first point is furthest for this facet

        bool
        further(G const & _nearer, G const & _further) const
        {
            if (outward_) {
                return (_nearer < _further);
            } else {
                return (_further < _nearer);
            }
        }

    };

    using facets_map = std::map< size_type, facet >;
    using facets_type = std::deque< size_type >;

    facets_map facets_;

    bool
    below(facet const & _facet, G const & _orientation) const
    {
        if (_facet.outward_) {
            return (G(0.0L) < _orientation);
        } else {
            return (_orientation < -G(0.0L));
        }
    }

    // http://math.stackexchange.com/questions/822741/
    template< typename vertices >
    G
    orientation(vertices const & _vertices, point_type const & _apex) const
    {
        size_type const size_ = _vertices.size(); // dimensionality of the subspace of interest
        assert(!(_apex.size() < size_));
        boost::numeric::ublas::matrix< G > m_(size_ + 1, size_ + 1);
        auto v_ = _vertices.cbegin();
        for (size_type i = 0; i < size_; ++i) {
            assert(v_ != _vertices.cend());
            point_type const & vertex_ = points_.at(*v_);
            ++v_;
            assert(!(vertex_.size() < size_));
            for (size_type j = 0; j < size_; ++j) {
                m_(i, j) = vertex_[j];
            }
            m_(i, size_) = G(1.0L);
        }
        for (size_type j = 0; j < size_; ++j) {
            m_(size_, j) = _apex[j];
        }
        m_(size_, size_) = G(1.0L);
        return determinant(std::move(m_));
    }

    template< typename vertices >
    G
    orientation(vertices const & _vertices, size_type const _apex) const
    {
        return orientation(_vertices, points_.at(_apex));
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

    G
    steal_best(point_list & _from, point_list & _to) const
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
        points_type & outside_set_ = _facet.outside_set_;
        G orientation_(0.0L);
        while (it != end) {
            auto const next = std::next(it);
            G const o_ = orientation(_facet.vertices_, *it);
            if (below(_facet, o_)) {
                if (outside_set_.empty() || _facet.further(orientation_, o_)) {
                    orientation_ = o_;
                    outside_set_.push_front(*it);
                } else {
                    outside_set_.push_back(*it);
                }
                _points.erase(it);
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
            facet & first_ = facets_.at(f);
            auto const lbeg = first_.points_.cbegin();
            auto const lend = first_.points_.cend();
            for (auto second = std::next(first); second != nend; ++second) {
                size_type const s = *second;
                facet & second_ = facets_.at(s);
                auto const rend = second_.points_.cend();
                auto r = second_.points_.cbegin();
                auto l = lbeg;
                bool lgood = false;
                bool rgood = false;
                while (l != lend) {
                    if (r == rend) {
                        lgood = (lgood != (++l == lend)); // xor
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
                rgood = (rgood != (++r == rend)); // xor
                if (lgood && rgood) {
                    first_.neighbours_.insert(s);
                    second_.neighbours_.insert(f);
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
        assert(vertices_.size() == 1 + dimension_); // (N + 1) vertices defining a simplex
        internal_set_.splice(internal_set_.end(), vertices_, vertices_.begin());
        assert(vertices_.size() == dimension_); // N vertices defining a facet
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
        auto const fend = facets_.end();
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) {
            size_type const n = facets_.size();
            auto const f = facets_.emplace_hint(fend, n, facet(vbeg, exclusive, vend, outward_));
            facet & facet_ = f->second;
            rank(partition(facet_, internal_set_), n);
            assert(outward_ == !(G(0.0L) < orientation(facet_.vertices_, inner_point_)));
            outward_ = !outward_;
        }
        assert(dimension_ + 1 == facets_.size()); // simplex
        {
            auto const fbeg = facets_.begin();
            for (auto i = fbeg; i != fend; ++i) {
                facet_set & neighbours_ = i->second.neighbours_;
                for (auto j = fbeg; j != fend; ++j) {
                    if (j != i) {
                        neighbours_.emplace_hint(neighbours_.end(), j->first);
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
        auto const fend = facets_.end();
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
                        facet & horizon_facet_ = facets_.at(n);
                        point_set horizon_(horizon_facet_.points_.cbegin(),
                                           horizon_facet_.points_.cend()); // `linear in N if the range is already sorted'
                        auto const hend = horizon_.end();
                        points_type ridge_; // horizon ridge + furthest point -> new facet
                        for (size_type const p : vertices_) { // facets intersection with keeping of points order as in visible facet
                            auto const h = horizon_.find(p);
                            if (h == hend) {
                                ridge_.push_back(apex);
                            } else {
                                ridge_.push_back(p);
                                horizon_.erase(h);
                            }
                        }
                        assert(horizon_.size() == 1); // horizon_ contains the only invisible point beyond the horizon
                        assert(ridge_.size() == dimension_); // ridge_ contains newfacet vertices (ridge + current furthest point)
                        { // replace visible facet became internal with newly created facet in adjacency
                            horizon_facet_.neighbours_.erase(v);
                            horizon_facet_.neighbours_.insert(horizon_facet_.neighbours_.cend(), facet_key);
                        }
                        newfacets_.push_back(facet_key);
                        facets_.emplace_hint(fend, facet_key, facet(std::move(ridge_), visible_facet_.outward_, n));
                        ++facet_key;
                    }
                }
            }
            adjacency(newfacets_);
            point_list outside_set_;
            for (size_type const v : visible_facets_) {
                auto const visible_facet = facets_.find(v);
                assert(visible_facet != fend);
                facet const & visible_facet_ = visible_facet->second;
                outside_set_.insert(outside_set_.cend(),
                                    visible_facet_.outside_set_.cbegin(),
                                    visible_facet_.outside_set_.cend());
                facets_.erase(visible_facet);
                unrank(v);
            }
            for (size_type const n : newfacets_) {
                rank(partition(facets_.at(n), outside_set_), n);
            }
            internal_set_.splice(internal_set_.cend(), outside_set_);
        }
    }

};
