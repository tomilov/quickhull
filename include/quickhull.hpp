#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <valarray>
#include <deque>
#include <set>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <utility>
#include <numeric>

#include <cassert>

struct bad_geometry
        : std::exception
{

    ~bad_geometry() noexcept = default;

    bad_geometry() = default;

    bad_geometry(const char * const _what)
        : what_(_what)
    { ; }

    virtual
    const char *
    what() const noexcept
    {
        return what_;
    }

private :

    const char * const what_ = "bad_get: failed value get using get()";

};

using boolean_type = bool;
using size_type = std::size_t;

template< typename G >
G
determinant(boost::numeric::ublas::matrix< G > _m)
{
    size_type const size_ = _m.size1();
    assert(_m.size2() == size_);
    boost::numeric::ublas::permutation_matrix< size_type > pm_(size_);
    if (0 < boost::numeric::ublas::lu_factorize(_m, pm_)) {
        return G(0.0L); // singular matrix
    } else {
        G determinant_(1.0L);
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

template< typename G >
struct convex_hull
{

    using self = convex_hull< G >;
    using point_type = std::valarray< G >;
    using points_type = std::deque< std::reference_wrapper< point_type const > >;
    using point_list = std::list< size_type >;
    using point_set = std::set< size_type >;

    convex_hull() = default;

    template< typename ForwardIterator >
    convex_hull(ForwardIterator _first, ForwardIterator _last)
        : dimension_(_first->size())
        , points_(_first, _last)
    {
        assert(0 < dimension_);
    }

    void
    append(point_type const & _point)
    {
        points_.push_back(_point);
    }

    void
    append(point_type && _point)
    {
        points_.push_back(std::move(_point));
    }

    size_type
    dimension() const
    {
        return dimension_;
    }

    boolean_type
    quickhull()
    {
        return {};
    }

    //private :

    size_type dimension_;
    points_type points_;

    struct hyperplane // oriented hyperplane
    {

        point_type unit_normal_;
        point_type offset_; // offset from the origin

    };

    struct facet;

    using facet_set = std::set< size_type >;

    struct facet // (d - 1)-dimensional faces
    {

        facet(point_list && _vertices,
              boolean_type const _outward)
            : vertices_(std::move(_vertices))
            , outward_(_outward)
        { ; }

        template< typename ForwardIterator >
        facet(ForwardIterator first, ForwardIterator last,
              boolean_type const _outward)
            : vertices_(first, last)
            , outward_(_outward)
        { ; }

        boolean_type
        further(G const & _nearer, G const & _further) const
        {
            if (outward_) {
                return (_nearer < _further);
            } else {
                return (_further < _nearer);
            }
        }

        boolean_type
        above(G const & _volume) const
        {
            if (outward_) {
                return (G(0.0L) < _volume);
            } else {
                return (_volume < -G(0.0L));
            }
        }

        point_list vertices_;
        boolean_type outward_;
        facet_set neighbours_;
        point_list outside_set_; // if not empty, then first point is best for this facet

    };

    using facets_type = std::map< size_type, facet >;
    using facet_iterator = typename facets_type::iterator;

    facets_type facets_;

    struct ridge // The (d - 2)-dimensional faces
    {

        points_type vertices_;

    };

    hyperplane
    get_oriented_hyperplane(points_type const & _points) const
    {
        return {};
    }

    G
    signed_distance_to_hyperplane(hyperplane const & _hyperplane, point_type const & _point) const
    {
        return (_hyperplane.unit_normal_ * (_point - _hyperplane.offset_)).sum();
    }

    // http://math.stackexchange.com/questions/822741/
    template< typename vertices >
    G
    volume(vertices const & _vertices, point_type const & _apex) const
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
    volume(vertices const & _vertices, size_type const _apex) const
    {
        return volume(_vertices, points_.at(_apex));
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
        G volume_ = volume(_to, *it);
        auto best = it;
        while (++it != end) {
            G const v_ = volume(_to, *it);
            if (abs(volume_) < abs(v_)) {
                volume_ = v_;
                best = it;
            }
        }
        if (!(G(0.0L) < abs(volume_))) {
            throw bad_geometry("can't find linearly independent point");
        }
        _to.splice(_to.end(), _from, best);
        return volume_;
    }

    points_type
    pricked_set(point_list const & _vertices, size_type const _nth) const
    {
        assert(_vertices.size() == dimension_ + 1);
        points_type pricked_set_;
        auto vertex_ = _vertices.cbegin();
        for (size_type i = 0; i <= dimension_; ++i) {
            if (i != _nth) {
                pricked_set_.emplace_back(*vertex_);
            }
            ++vertex_;
        }
        return pricked_set_;
    }

    G
    partition(facet & _facet, point_list & _points)
    {
        auto it = _points.begin();
        auto const end = _points.end();
        point_list & outside_set_ = _facet.outside_set_;
        auto const oend = outside_set_.end(); // remains valid for std::list
        G volume_(0.0L);
        while (it != end) {
            auto const next = std::next(it);
            G const v_ = volume(_facet.vertices_, *it);
            if (_facet.above(v_)) {
                if (outside_set_.empty() || _facet.further(volume_, v_)) {
                    volume_ = v_;
                    outside_set_.splice(outside_set_.begin(), _points, it);
                } else {
                    outside_set_.splice(oend, _points, it);
                }
            }
            it = next;
        }
        return volume_;
    }

    template< typename ForwardIterator >
    facet
    make_facet(ForwardIterator beg, ForwardIterator mid, ForwardIterator end,
               boolean_type const _outward)
    {
        facet facet_(beg, std::prev(mid), _outward);
        point_list & vertices_ = facet_.vertices_;
        vertices_.insert(vertices_.end(), mid, end);
        return facet_;
    }

    facet_iterator
    create_simplex()
    {
        point_list point_list_;
        {
            size_type const size_ = points_.size();
            assert(dimension_ < size_);
            for (size_type i = 0; i < size_; ++i) {
                point_list_.push_back(i);
            }
        }
        point_list vertices_;
        vertices_.splice(vertices_.end(), point_list_, point_list_.begin());
        for (size_type i = 0; i < dimension_; ++i) {
            steal_best(point_list_, vertices_);
        }
        assert(vertices_.size() == 1 + dimension_); // (N + 1) vertices defines a simplex
        point_list_.splice(point_list_.end(), vertices_, vertices_.begin());
        assert(vertices_.size() == dimension_); // N vertices defines a facet
        boolean_type outward_ = !(G(0.0L) < steal_best(point_list_, vertices_)); // top oriented?
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
        auto furthest = fend;
        G volume_(0.0L);
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) {
            auto const f = facets_.emplace_hint(fend, facets_.size(), make_facet(vbeg, exclusive, vend, outward_));
            facet & facet_ = f->second;
            G const v_ = abs(partition(facet_, point_list_));
            if (volume_ < v_) {
                volume_ = v_;
                furthest = f;
            }
            assert(outward_ == !(G(0.0L) < volume(facet_.vertices_, inner_point_)));
            outward_ = !outward_;
        }
        //point_list_.clear();
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
        if (furthest == fend) {
            std::cout << "convex hull is simplex" << std::endl; // rbox D3 5 n t5 > points.txt
        }
        return furthest;
    }

    void
    create_convex_hull()
    {
        size_type facet_key_ = facets_.size(); // unique key for facets_
        auto const fend = facets_.end();
        auto furthest = create_simplex();
        while (furthest != fend) {
            facet const & facet_ = furthest->second;
            facet_set visible_facets_{furthest->first};
            facet_set neighbours_ = facet_.neighbours_;
            point_list const & outside_set_ = facet_.outside_set_;
            size_type const apex = outside_set_.front();
            while (!neighbours_.empty()) {
                auto const first = neighbours_.begin();
                size_type const f = *first;
                auto const candidate = facets_.find(f);
                assert(candidate != facets_.end());
                facet const & candidate_ = candidate->second;
                if (G(0.0L) < volume(candidate_.vertices_, apex)) { // if point is above the neighbour, then add they to visible set
                    visible_facets_.insert(f);
                    neighbours_.insert(candidate_.neighbours_.cbegin(),
                                       candidate_.neighbours_.cend());
                }
                neighbours_.erase(first);
            }
            // the boundary of visible_facets_ is the set of horizon ridges
            // Each ridge signifies the adjacency of two facets.
            // ... qh_makenew_simplicial ...
            auto const vfend = visible_facets_.end();
            for (size_type const v : visible_facets_) {
                assert(facets_.find(v) != facets_.end());
                facet const & visible_facet_ = facets_.at(v);
                point_list const & vertices_ = visible_facet_.vertices_;
                for (size_type const n : visible_facet_.neighbours_) {
                    if (visible_facets_.find(n) == vfend) { // facets intersection with keeping of order
                        facet & horizon_facet_ = facets_.at(n);
                        point_set horizon_(horizon_facet_.vertices_.cbegin(),
                                           horizon_facet_.vertices_.cend()); // n * log(n) +
                        auto const hend = horizon_.end();
                        point_list ridge_;
                        for (size_type const p : vertices_) { // n *
                            auto const h = horizon_.find(p); // (log(n) +
                            if (h == hend) {
                                ridge_.push_back(apex);
                            } else {
                                ridge_.push_back(*h);
                                horizon_.erase(h); // const)
                            }
                        }
                        assert(horizon_.size() == 1);
                        assert(ridge_.size() == vertices_.size());
                        // ... qh_makenewfacet ...
                        { // replace visible facet became internal with newly created facet
                            horizon_facet_.neighbours_.erase(v);
                            horizon_facet_.neighbours_.insert(facet_key_);
                        }
                        auto const r = facets_.emplace_hint(fend, facet_key_, facet(std::move(ridge_), visible_facet_.outward_));
                        ++facet_key_;
                    }
                }
            }
        }
    }

};
