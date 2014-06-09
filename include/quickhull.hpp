#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <valarray>
#include <deque>
#include <set>
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
        : dimension_(*_first.size())
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

    struct facet // (d - 1)-dimensional faces
    {

        points_type vertices_;
        std::deque< std::reference_wrapper< facet const > > neighboring_facets_;
        hyperplane hyperplane_; // hyperplane equation

    };

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

    G
    orientation(point_list const & _vertices, point_type const & _point) const
    {
        size_type const size_ = _vertices.size(); // dimension of the subspace of interest
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
    }

    struct bad_geometry
    {

    };

    size_type
    random(size_type const) const
    {
        return 0;
    }

    void
    steal_farthest(point_list & _from, point_list & _to) const
    {
        auto it = _from.begin();
        auto const end = _from.end();
        G orientation_ = orientation(_to, *it);
        auto farthest = it;
        while (++it != end) {
            if (orientation_ < orientation(_to, *it)) {
                farthest = it;
            }
        }
        _to.splice(_to.end(), _from, farthest);
        return (orientation_ != G(0.0L));
    }

    point_list
    create_simplex() const
    {
        assert(dimension_ < points_.size());
        point_list points_list_(points_.cbegin() + 1, points_.cend());
        point_list simplex_{points_.front()};
        for (size_type i = 0; i < dimension_; ++i) {
            if (!steal_farthest(points_list_, simplex_)) {
                throw bad_geometry(); // linear dependent points
            }
        }
        return simplex_;
    }

};
