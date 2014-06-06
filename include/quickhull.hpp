#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <valarray>
#include <deque>
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
        bool sign_ = false;
        for (size_type i = 0; i < size_; ++i) {
            if (i != pm_(i)) {
                sign_ = !sign_;
            }
            determinant_ *= _m(i, i);
        }
        if (sign_) {
            return -determinant_;
        } else {
            return +determinant_;
        }
    } else {
        return G(0.0L); // singular matrix
    }
}

template< typename G >
using point = std::valarray< G >;

template< typename G >
struct convex_hull
{

    using self = convex_hull< G >;
    using point_type = point< G >;
    using points_type = std::deque< point_type >;

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

    }

private :

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

    points_type
    create_simplex() const
    {
        points_type simplex_;
        std::copy_n(points_, dimension_, std::back_inserter(simplex_));

        return simplex_;
    }

};
