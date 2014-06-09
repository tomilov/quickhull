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

using boolean_type = bool;
using size_type = std::size_t;

template< typename G >
G
determinant(boost::numeric::ublas::matrix< G > _m)
{
    size_type const size_ = _m.size1();
    assert(_m.size2() == size_);
    boost::numeric::ublas::permutation_matrix< size_type > pm_(size_);
    if (boost::numeric::ublas::lu_factorize(_m, pm_) == 0) {
        G determinant_(1.0L);
        for (size_type i = 0; i < size_; ++i) {
            if (i == pm_(i)) {
                determinant_ *= +_m(i, i);
            } else {
                determinant_ *= -_m(i, i);
            }
        }
        return determinant_;
    } else {
        return G(0.0L); // singular matrix
    }
}

template< typename G >
struct convex_hull
{

    using self = convex_hull< G >;
    using point_type = std::valarray< G >;
    using points_type = std::deque< std::reference_wrapper< point_type const > >;
    using point_list = std::list< std::reference_wrapper< point_type const > >;

    convex_hull() = default;

    template< typename ForwardIterator >
    convex_hull(ForwardIterator _first, ForwardIterator _last)
        : dimension_(_first->size())
        , points_(_first, _last)
    { ; }

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

    using facets_type = std::map< size_type, facet >;
    size_type facet_id_ = 0;

    using facet_set_type = std::set< size_type >;

    struct facet // (d - 1)-dimensional faces
    {

        facet() = default;

        points_type vertices_;
        boolean_type upward_;
        facet_set_type neighboring_facets_;
        point_list outside_set_; // if not empty, then first point is furthest from this facet

    };

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

    template< typename vertices > // vertices is point_list or points_type
    G
    orientation(vertices const & _vertices, point_type const & _point) const
    {
        /*
            $\displaystyle
            \Delta_{H}(x)=
            \left\vert
            \begin{array}{ccccc}
             p_{1,1} & p_{1,2} & \cdots & p_{1,n} & 1 \\
             p_{2,1} & p_{2,2} & \cdots & p_{2,n} & 1 \\
             \vdots & \vdots & \vdots & \vdots & \vdots \\
             p_{n,1} & p_{n,2} & \cdots & p_{n,n} & 1 \\
             x_{1} & x_{2} & \cdots & x_{n} & 1
            \end{array}
            \right\vert$
         */
        size_type const size_ = _vertices.size(); // dimensionality of the subspace of interest
        assert(!(_point.size() < size_));
        boost::numeric::ublas::matrix< G > m_(size_ + 1, size_ + 1);
        auto v_ = _vertices.cbegin();
        for (size_type i = 0; i < size_; ++i) {
            assert(v_ != _vertices.cend());
            point_type const & vertex_ = *v_;
            ++v_;
            assert(!(vertex_.size() < size_));
            for (size_type j = 0; j < size_; ++j) {
                m_(i, j) = vertex_[j];
            }
            m_(i, size_) = G(1.0L);
        }
        for (size_type j = 0; j < size_; ++j) {
            m_(size_, j) = _point[j];
        }
        m_(size_, size_) = G(1.0L);
        return determinant(std::move(m_));
        // if orientation is less than (dimension_-th-root of epsilon), then we have a coplanar set of points
    }

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

    G
    abs(G const & _x) const
    {
        return (_x < G(0.0L)) ? -_x : +_x;
    }

    G
    abs(G && _x) const
    {
        return (_x < G(0.0L)) ? -std::move(_x) : std::move(_x);
    }

    G
    steal_furthest(point_list & _from, point_list & _to) const // move from from_ to to_ furthest (in sense of _to.size()-dimensional subspace distance) point
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

    points_type
    hole_set(point_list const & _vertices, size_type const _nth) const
    {
        assert(_vertices.size() == dimension_ + 1);
        points_type hole_set_;
        auto vertex_ = _vertices.cbegin();
        for (size_type i = 0; i <= dimension_; ++i) {
            if (i != _nth) {
                hole_set_.emplace_back(*vertex_);
            }
            ++vertex_;
        }
        return hole_set_;
    }

    void
    create_simplex()
    {
        assert(dimension_ < points_.size());
        point_list vertices_{points_.front()}; // not the optimal choise, but would then be rejudged
        point_list source_points_(std::next(points_.cbegin()), points_.cend());
        {
            for (size_type i = 0; i < dimension_; ++i) {
                G const orientation_ = steal_furthest(source_points_, vertices_);
                if (!(G(0.0L) < abs(orientation_))) {
                    throw bad_geometry("can't select a (dim + 1) set of noncoplanar points");
                }
            }
            for (size_type i = 0; i < dimension_; ++i) { // rejudge feasibility of all the points again
                source_points_.splice(source_points_.end(), vertices_, vertices_.begin());
                G const orientation_ = steal_furthest(source_points_, vertices_);
                if (!(G(0.0L) < abs(orientation_))) {
                    throw bad_geometry("can't select a (dim + 1) set of noncoplanar points");
                }
            }
            { // last
                source_points_.splice(source_points_.end(), vertices_, vertices_.begin());
                points_type first_vertices_(vertices_.cbegin(), vertices_.cend());
                G const o_ = steal_furthest(source_points_, vertices_);
                if (!(G(0.0L) < abs(o_))) {
                    throw bad_geometry("can't select a (dim + 1) set of noncoplanar points");
                }
#ifndef NDEBUG
                point_type inner_point_;
                {
                    auto it = vertices_.cbegin();
                    inner_point_ = *it;
                    auto const end = vertices_.cend();
                    while (++it != end) {
                        inner_point_ += *it;
                    }
                    inner_point_ /= G(1 + dimension_);
                }
#endif
                vertices_.splice(vertices_.begin(), vertices_, std::prev(vertices_.end()));
                boolean_type upward_ = (G(0.0L) < o_);
                facets_.emplace_hint(facets_.end(), 0, facet{std::move(first_vertices_), upward_});
                for (size_type i = 1; i <= dimension_; ++i) {
                    upward_ = !upward_;
                    facets_.emplace_hint(facets_.end(), i, facet{hole_set(vertices_, i), upward_});
                    assert(upward_ == (G(0.0L) < orientation(std::prev(facets_.end())->second.vertices_, inner_point_)));
                }
            }
        }
#if 0
        {
            auto const beg = facets_.begin();
            auto const end = facets_.end();
            for (auto i = beg; i != end; ++i) {
                facet_set_type & neighboring_facets_ = i->second.neighboring_facets_;
                for (auto j = beg; j != end; ++j) {
                    if (j != i) {
                        neighboring_facets_.emplace_hint(neighboring_facets_.end(), j->first);
                    }
                }
            }
            // set hyperplane equation for each facet here (does we need in?)
            for (auto & facet_ : facets_) {
                auto it = source_points_.begin();
                auto const end = source_points_.end();
                point_list & outside_set_ = facet_.second.outside_set_;
                auto const oend = outside_set_.end(); // remains valid for std::list
                G orientation_(0.0L);
                while (it != end) {
                    auto const next = std::next(it);
                    G const o_ = orientation(facet_.second.vertices_, *it);
                    if (G(0.0L) < o_) {
                        if (outside_set_.empty() || (orientation_ < o_)) {
                            orientation_ = o_;
                            outside_set_.splice(outside_set_.begin(), source_points_, it);
                        } else {
                            outside_set_.splice(oend, source_points_, it);
                        }
                    }
                    it = next;
                }
                if (source_points_.empty()) {
                    break;
                }
            }
        }
        for (;;) {
            auto current = facets_.begin();
            auto const end = facets_.end();
            while (current != end) {
                if (!current->second.outside_set_.empty()) {
                    facet const & facet_ = current->second;
                    point_list const & outside_set_ = facet_.outside_set_;
                    if (!outside_set_.empty()) {
                        point_type const & furthest_point_ = outside_set_.front();
                        facet_set_type visible_facets_{current->first};
                        facet_set_type neighboring_facets_ = facet_.neighboring_facets_;
                        while (!neighboring_facets_.empty()) {
                            auto const first = neighboring_facets_.begin();
                            size_type const f_ = *first;
                            auto const candidate = facets_.find(f_);
                            assert(candidate != facets_.end());
                            facet const & candidate_facet_ = candidate->second;
                            if (G(0.0L) < orientation(candidate_facet_.vertices_, furthest_point_)) { // if point is above the neighbor, then add they to visible set
                                visible_facets_.insert(f_);
                                neighboring_facets_.insert(candidate_facet_.neighboring_facets_.cbegin(),
                                                           candidate_facet_.neighboring_facets_.cend());
                            }
                            neighboring_facets_.erase(first);
                        }
                    }
                }
            }
            // ?
        }
#endif
    }

};
