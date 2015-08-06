/* Quickhull algorithm implementation
 *
 * Copyright (c) 2014-2015, Anatoliy V. Tomilov
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following condition is met:
 * Redistributions of source code must retain the above copyright notice, this condition and the following disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once

#include <valarray>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <utility>
#include <numeric>
#include <limits>
#include <functional>
#ifdef _DEBUG
#include <iostream>
#endif

#include <cstdint>
#include <cmath>
#include <cassert>

template< typename point_iterator >
struct quick_hull
{

    static_assert(std::is_base_of< std::random_access_iterator_tag, typename std::iterator_traits< point_iterator >::iterator_category >::value);

    using size_type = std::size_t;
    using difference_type = std::intptr_t;

    using point = typename std::iterator_traits< point_iterator >::value_type;
    using value_type = std::decay_t< decltype(*std::cbegin(*std::declval< point_iterator >())) >;

    size_type const dimension_;
    value_type const eps;

    value_type const zero = value_type(0);
    value_type const one = value_type(1);

    using vector = std::valarray< value_type >;

private :

    using matrix = std::vector< vector >;

    matrix matrix_;
    matrix shadow_matrix_;
    matrix minor_;

public :

    quick_hull(size_type const _dimension,
               value_type _eps = std::numeric_limits< value_type >::epsilon())
        : dimension_(_dimension)
        , eps(std::move(_eps))
        , matrix_(dimension_)
        , shadow_matrix_(dimension_)
        , minor_(dimension_)
    {
        assert(1 < dimension_);
        assert(!(eps < zero));
        for (size_type r = 0; r < dimension_; ++r) {
            matrix_[r].resize(dimension_);
            shadow_matrix_[r].resize(dimension_);
        }
        size_type const minor_size = dimension_ - 1;
        for (size_type r = 1; r < minor_size; ++r) {
            minor_[r].resize(minor_size);
        }
        minor_.front().resize(dimension_);
        minor_.back().resize(dimension_);
    }

    using point_array = std::vector< point_iterator >;
    using point_deque = std::deque< point_iterator >;
    using point_list = std::list< point_iterator >;
    using facet_array = std::vector< size_type >;

    struct facet // (d - 1)-dimensional face
    {

        // each neighbouring facet lies against corresponding vertex and vice versa
        point_array vertices_; // dimension_ points (oriented)
        facet_array neighbours_; // dimension_ neighbouring facets

        point_list outside_; // if empty, then is convex hull's facet, else the first point (i.e. outside_.front()) is the furthest point from this facet
        point_deque coplanar_; // coplanar points, for resulting convex hull it is guaranted that they lies within the facet or on a facet's ridge (in later case these points can be non-unique)

        // equation of hyperplane supported the facet
        vector normal_; // components of normalized normal vector
        value_type D; // distance from the origin to the hyperplane

        void
        init(size_type const _dimension,
             point_array const & _vertices,
             size_type const _against,
             point_iterator const & _apex,
             size_type const _neighbour)
        {
            static_cast< void >(_dimension);
            assert(vertices_.size() == _dimension);
            vertices_ = _vertices;
            vertices_[_against] = _apex;
            assert(neighbours_.size() == _dimension);
            neighbours_[_against] = _neighbour;
            assert(normal_.size() == _dimension);
        }

        facet(size_type const _dimension,
              point_array const & _vertices,
              size_type const _against,
              point_iterator const & _apex,
              size_type const _neighbour)
            : vertices_(_vertices)
            , neighbours_(_dimension)
            , normal_(_dimension)
        {
            assert(vertices_.size() == _dimension);
            vertices_[_against] = _apex;
            neighbours_[_against] = _neighbour;
        }

        facet(size_type const _dimension,
              point_array const & _simplex,
              size_type const _vertex)
            : normal_(_dimension)
        {
            assert(_simplex.size() == _dimension + 1);
            vertices_.reserve(_dimension);
            neighbours_.reserve(_dimension);
            for (size_type v = 0; v <= _dimension; ++v) {
                if (v != _vertex) {
                    vertices_.push_back(_simplex[v]);
                    neighbours_.push_back(v);
                }
            }
            if (((_dimension - _vertex) % 2) == 0) {
                using std::swap;
                swap(vertices_.front(), vertices_.back());
                swap(neighbours_.front(), neighbours_.back());
            }
        }

        value_type
        distance(point const & _point) const
        {
            return std::inner_product(std::cbegin(normal_), std::cend(normal_), std::cbegin(_point), D);
        }

    };

    using facets = std::deque< facet >;

    facets facets_;

    value_type
    cos_of_dihedral_angle(facet const & _this, facet const & _other) const
    {
        return std::inner_product(std::cbegin(_this.normal_), std::cend(_this.normal_), std::cbegin(_other.normal_), zero);
    }

private :

    void
    transpose() // transpose shadow matrix to cheaper filling columns with 1-s
    {
        for (size_type r = 0; r < dimension_; ++r) {
            vector & row_ = shadow_matrix_[r];
            for (size_type c = 1 + r; c < dimension_; ++c) {
                using std::swap;
                swap(shadow_matrix_[c][r], row_[c]);
            }
        }
    }

    void
    restore_matrix(size_type const _identity) // load matrix from storage, but replace _identity column with 1-s
    {
        for (size_type c = 0; c < dimension_; ++c) {
            vector & col_ = matrix_[c];
            if (c == _identity) {
                col_ = one;
            } else {
                col_ = shadow_matrix_[c];
            }
        }
    }

    void
    square_matrix(size_type const _size) // matrix_ = shadow_matrix_ * transposed shadow_matrix_
    {
        assert(_size < dimension_);
        for (size_type r = 0; r < _size; ++r) {
            vector & lhs_ = shadow_matrix_[r];
            vector const & row_ = matrix_[r];
            for (size_type c = 0; c < _size; ++c) {
                lhs_[c] = std::inner_product(std::cbegin(row_), std::cend(row_), std::cbegin(matrix_[c]), zero);
            }
        }
    }

    // based on LUP decomposition (complexity is (n^3 / 3 + n^2 / 2 - 5 * n / 6) vs (2 * n^3 / 3 + n^2 + n / 3 - 2) for QR decomposition via Householder reflections) http://math.stackexchange.com/a/93508/54348
    value_type
    det(matrix & _matrix, size_type const _dimension) // hottest function
    { // produces lower unit triangular matrix and upper triangular
        assert(0 < _dimension);
        value_type det_ = one;
        for (size_type i = 0; i < _dimension; ++i) {
            vector & ri_ = _matrix[i];
            using std::abs;
            value_type max_ = abs(ri_[i]);
            size_type pivot = i;
            {
                size_type p = i;
                while (++p < _dimension) {
                    value_type y_ = abs(_matrix[p][i]);
                    if (max_ < y_) {
                        max_ = std::move(y_);
                        pivot = p;
                    }
                }
            }
            if (!(eps < max_)) { // regular?
                return zero; // singular
            }
            if (pivot != i) {
                det_ = -det_; // each permutation flips sign of det
                ri_.swap(_matrix[pivot]);
            }
            value_type & dia_ = ri_[i];
            det_ *= dia_; // det is multiple of diagonal elements
            for (size_type j = 1 + i; j < _dimension; ++j) {
                _matrix[j][i] /= dia_;
            }
            for (size_type a = 1 + i; a < _dimension; ++a) {
                vector & a_ = minor_[a - 1];
                value_type const & ai_ = _matrix[a][i];
                for (size_type b = 1 + i; b < _dimension; ++b) {
                    a_[b - 1] = ai_ * ri_[b];
                }
            }
            for (size_type a = 1 + i; a < _dimension; ++a) {
                vector const & a_ = minor_[a - 1];
                vector & ra_ = _matrix[a];
                for (size_type b = 1 + i; b < _dimension; ++b) {
                    ra_[b] -= a_[b - 1];
                }
            }
        }
        return det_;
    }

    value_type
    det()
    {
        return det(matrix_, dimension_);
    }

    void
    set_hyperplane_equation(facet & _facet)
    {
        for (size_type r = 0; r < dimension_; ++r) {
            std::copy_n(std::cbegin(*_facet.vertices_[r]), dimension_, std::begin(shadow_matrix_[r]));
        }
        transpose();
        matrix_ = shadow_matrix_;
        _facet.D = -det();
        value_type N = zero;
        for (size_type i = 0; i < dimension_; ++i) {
            restore_matrix(i);
            value_type & n = _facet.normal_[i];
            n = det();
            N += n * n;
        }
        using std::sqrt;
        N = sqrt(std::move(N));
        _facet.normal_ /= N;
        _facet.D /= std::move(N);
    }

    bool
    orthonormalize(point_array const & _affine_space, size_type const _rank, vector const & _origin)
    {
        assert(!(dimension_ < _rank));
        assert(!(_affine_space.size() < _rank));
        auto vertex = std::begin(_affine_space);
        for (size_type r = 0; r < _rank; ++r) { // affine space -> vector space
            vector & row_ = shadow_matrix_[r];
            std::copy_n(std::cbegin(**vertex), dimension_, std::begin(row_));
            row_ -= _origin;
            ++vertex;
        }
        for (size_type i = 0; i < _rank; ++i) { // Householder transformation
            value_type norm_ = zero;
            vector & qri_ = shadow_matrix_[i]; // shadow_matrix_ is packed QR after all
            for (size_type j = i; j < dimension_; ++j) {
                value_type const & qrij_ = qri_[j];
                norm_ += qrij_ * qrij_;
            }
            using std::sqrt;
            norm_ = sqrt(norm_);
            if (!(eps < norm_)) {
                return false;
            }
            value_type & qrii_ = qri_[i];
            bool const sign_ = (zero < qrii_);
            value_type factor_ = norm_ * (norm_ + (sign_ ? qrii_ : -qrii_));
            if (!(eps < factor_)) {
                return false;
            }
            factor_ = one / sqrt(std::move(factor_));
            if (sign_) {
                qrii_ += norm_;
            } else {
                qrii_ -= norm_;
            }
            for (size_type k = i; k < dimension_; ++k) {
                qri_[k] *= factor_;
            }
            for (size_type j = i + 1; j < _rank; ++j) {
                vector & qrj_ = shadow_matrix_[j];
                value_type s_ = zero;
                for (size_type k = i; k < dimension_; ++k) {
                    s_ += qri_[k] * qrj_[k];
                }
                for (size_type k = i; k < dimension_; ++k) {
                    qrj_[k] -= qri_[k] * s_;
                }
            }
        }
        return true;
    }

    void
    forward_transformation(size_type const _rank) // calculation of Q
    {
        assert(!(dimension_ < _rank));
        for (size_type i = 0; i < _rank; ++i) {
            vector & qi_ = matrix_[i]; // matrix_ is Q after all
            qi_ = zero;
            qi_[i] = one;
            size_type j = _rank;
            while (0 < j) {
                --j;
                vector & qrj_ = shadow_matrix_[j]; // containing packed QR
                value_type s_ = zero;
                for (size_type k = j; k < dimension_; ++k) {
                    s_ += qrj_[k] * qi_[k];
                }
                for (size_type k = j; k < dimension_; ++k) {
                    qi_[k] -= qrj_[k] * s_;
                }
            }
        }
    }

    bool
    steal_best(point_list & _from, point_array & _to)
    {
        assert(!_to.empty());
        size_type const rank_ = _to.size() - 1;
        assert(rank_ < dimension_);
        vector & origin_ = matrix_[rank_];
        std::copy_n(std::cbegin(*_to.back()), dimension_, std::begin(origin_));
        if (!orthonormalize(_to, rank_, origin_)) {
            return false;
        }
        forward_transformation(rank_);
        vector & projection_ = minor_.back();
        vector & apex_ = minor_.front();
        value_type distance_ = zero;
        auto furthest = std::cend(_from);
        for (auto it = std::cbegin(_from); it != std::cend(_from); ++it) {
            std::copy_n(std::cbegin(**it), dimension_, std::begin(apex_));
            apex_ -= origin_; // turn translated space into vector space
            projection_ = apex_; // project onto orthogonal subspace
            for (size_type i = 0; i < rank_; ++i) {
                vector const & qi_ = matrix_[i];
                projection_ -= std::inner_product(std::cbegin(apex_), std::cend(apex_), std::cbegin(qi_), zero) * qi_;
            }
            projection_ *= projection_;
            using std::sqrt;
            value_type d_ = sqrt(projection_.sum()); // distance to the subspace
            if (distance_ < d_) {
                distance_ = std::move(d_);
                furthest = it;
            }
        }
        if (furthest == std::cend(_from)) {
            return false;
        }
        _to.push_back(std::move(*furthest));
        _from.erase(furthest);
        return true;
    }

    std::set< size_type, std::greater< size_type > > removed_facets_;

    size_type
    add_facet(point_array const & _vertices, size_type const _against, point_iterator const & _apex, size_type const _neighbour)
    {
        if (removed_facets_.empty()) {
            size_type const f = facets_.size();
            facets_.emplace_back(dimension_, _vertices, _against, _apex, _neighbour);
            facet & facet_ = facets_.back();
            set_hyperplane_equation(facet_);
            return f;
        } else {
            auto const rend = std::prev(std::cend(removed_facets_));
            size_type const f = *rend;
            removed_facets_.erase(rend);
            facet & facet_ = facets_[f];
            facet_.init(dimension_, _vertices, _against, _apex, _neighbour);
            set_hyperplane_equation(facet_);
            return f;
        }
    }

    // selecting of the best facet:

    using ranking = std::multimap< value_type, size_type >;
    using ranking_meta = std::unordered_map< size_type, typename ranking::iterator >;

    ranking ranking_;
    ranking_meta ranking_meta_;

    void
    rank(value_type && _orientation, size_type const _facet)
    {
        if (eps < _orientation) {
            ranking_meta_.emplace(_facet, ranking_.emplace(std::move(_orientation), _facet));
        }
    }

    void
    unrank(size_type const _facet)
    {
        auto const r = ranking_meta_.find(_facet);
        if (r != std::end(ranking_meta_)) {
            ranking_.erase(r->second);
            ranking_meta_.erase(r);
        }
        removed_facets_.insert(_facet);
    }

    value_type
    partition(facet & _facet, point_list & _points)
    {
        auto it = std::cbegin(_points);
        value_type distance_ = zero;
        while (it != std::cend(_points)) {
            auto const next = std::next(it);
            value_type d_ = _facet.distance(**it);
            if (eps < d_) {
                if (distance_ < d_) {
                    distance_ = std::move(d_);
                    _facet.outside_.splice(std::cbegin(_facet.outside_), _points, it);
                } else {
                    _facet.outside_.splice(std::cend(_facet.outside_), _points, it);
                }
            } else if (!(d_ < -eps)) {
                _facet.coplanar_.push_back(*it);
            }
            it = next;
        }
        return distance_;
    }

    size_type
    get_best_facet() const // select a facet with furthest (between all facets with non-empty outsides_ set) point
    {
        assert(ranking_meta_.size() == ranking_.size());
        return std::prev(std::cend(ranking_))->second;
    }

    // visibility from apex:

    using facet_unordered_set = std::unordered_set< size_type >;

    facet_unordered_set visited_;
    facet_unordered_set visible_;

    bool
    process_visibles(point_list & _outside, facet_array & _newfacets, size_type const _facet, point_iterator const & _apex) // traverse the graph of visible facets
    {
        if (!visited_.insert(_facet).second) {
            return (visible_.count(_facet) != 0);
        }
        facet & facet_ = facets_[_facet];
        if (!(eps < facet_.distance(*_apex))) {
            return false;
        }
        visible_.insert(_facet);
        _outside.splice(std::cend(_outside), std::move(facet_.outside_));
        facet_.coplanar_.clear();
        for (size_type vertex = 0; vertex < dimension_; ++vertex) {
            size_type const neighbour = facet_.neighbours_[vertex];
            if (!process_visibles(_outside, _newfacets, neighbour, _apex)) {
                size_type const newfacet = add_facet(facet_.vertices_, vertex, _apex, neighbour);
                _newfacets.push_back(newfacet);
                replace_neighbour(neighbour, _facet, newfacet);
                find_adjacent_facets(newfacet, vertex, _apex);
            }
        }
        // to track internal_set_: move all (not yet contained in internal_set_) facet_.vertices_ to internal_set_ if all the neighbours are visible
        unrank(_facet);
        return true;
    }

    void
    replace_neighbour(size_type const _facet, size_type const _from, size_type const _to)
    {
        if (_from != _to) {
            for (size_type & neighbour : facets_[_facet].neighbours_) {
                if (neighbour == _from) {
                    neighbour = _to;
                    return;
                }
            }
        }
    }

    // adjacency of new facets via its common ridges:

    struct ridge
    {

        facet & facet_;
        size_type const f_;
        size_type const against_;
        size_type const hash_;

        bool
        operator == (ridge const & _rhs) const noexcept
        {
            point_iterator const & lskip = facet_.vertices_[against_];
            point_iterator const & rskip = _rhs.facet_.vertices_[_rhs.against_];
            for (point_iterator const & l : facet_.vertices_) {
                if (l != lskip) {
                    bool found_ = false;
                    for (point_iterator const & r : _rhs.facet_.vertices_) {
                        if (r != rskip) {
                            if (l == r) {
                                found_ = true; // O(d^2) expensive
                                break;
                            }
                        }
                    }
                    if (!found_) {
                        return false;
                    }
                }
            }
            return true;
        }

    };

    struct ridge_hash
    {

        size_type
        operator () (ridge const & _ridge) const noexcept
        {
            return _ridge.hash_;
        }

    };

    std::unordered_set< ridge, ridge_hash > unique_ridges_;

    void
    find_adjacent_facets(size_type const _f, size_type const _skip, point_iterator const & _apex)
    {
        facet & facet_ = facets_[_f];
        std::hash< difference_type > point_hash_;
        size_type ridge_hash_ = 0;
        for (size_type v = 0; v < dimension_; ++v) {
            if (v != _skip) {
                ridge_hash_ ^= point_hash_(facet_.vertices_[v] - _apex);
            }
        }
        for (size_type against_ = 0; against_ < dimension_; ++against_) {
            if (against_ != _skip) { // neighbouring facet against _apex is known atm
                auto position = unique_ridges_.insert({facet_, _f, against_, (ridge_hash_ ^ point_hash_(facet_.vertices_[against_] - _apex))});
                if (!position.second) {
                    ridge const & ridge_ = *position.first;
                    ridge_.facet_.neighbours_[ridge_.against_] = _f;
                    facet_.neighbours_[against_] = ridge_.f_;
                    unique_ridges_.erase(position.first);
                }
            }
        }
    }

public : // hypervolume of simplex, largest possible simplex heuristic, convex hull algorithm:

    // http://math.stackexchange.com/questions/822741/
    value_type
    hypervolume(point_array const & _vertices) // hypervolume of parallelotope spanned on vectors from one of _vertices to all the rest
    {
        assert(!_vertices.empty());
        size_type const rank_ = _vertices.size() - 1;
        assert(!(dimension_ < rank_));
        vector & origin_ = minor_.back();
        std::copy_n(std::cbegin(*_vertices.back()), dimension_, std::begin(origin_));
        auto vertex = std::cbegin(_vertices);
        for (size_type r = 0; r < rank_; ++r) { // affine space -> vector space
            vector & row_ = matrix_[r];
            std::copy_n(std::cbegin(**vertex), dimension_, std::begin(row_));
            row_ -= origin_;
            ++vertex;
        }
        if (rank_ == dimension_) { // oriented hypervolume
            return det();
        } else { // non-oriented _rank-dimensional measure
            square_matrix(rank_);
            using std::sqrt;
            return sqrt(det(shadow_matrix_, rank_));
        }
    }

    point_array
    create_initial_simplex(point_iterator const _beg, point_iterator const _end)
    {
        // construction of reference set of universe of points
        point_array basis_;
        if (_beg == _end) {
            return basis_;
        }
        basis_.reserve(dimension_ + 1);
        basis_.push_back(_beg);
        point_list internal_set_;
        {
            auto it = _beg;
            while (++it != _end) {
                internal_set_.push_back(it);
            }
        }
        // selection of (dimension_ + 1) affinely independent points
        if (!steal_best(internal_set_, basis_)) {
            return basis_; // can't find affinely independent second point
        }
        { // reject 0-indexed point to rejudge it
            point_iterator & first_ = basis_.front();
            internal_set_.push_back(first_);
            first_ = std::move(basis_.back());
            basis_.pop_back();
        }
        for (size_type i = 0; i < dimension_; ++i) {
            if (!steal_best(internal_set_, basis_)) {
                return basis_; // can't find (i + 2) affinely independent points
            }
        }
        assert(basis_.size() == dimension_ + 1); // enough for simplex construction
        // simplex construction
        if (zero < hypervolume(basis_)) {
            using std::swap;
            swap(basis_.front(), basis_.back());
        }
        for (size_type newfacet = 0; newfacet <= dimension_; ++newfacet) {
            facets_.emplace_back(dimension_, basis_, newfacet);
            facet & newfacet_ = facets_.back();
            set_hyperplane_equation(newfacet_);
            rank(partition(newfacet_, internal_set_), newfacet);
        }
        return basis_;
    }

    // Barber, C. B., D.P. Dobkin, and H.T. Huhdanpaa, 1995. "The Quickhull Algorithm for Convex Hulls", ACM Transactions on Mathematical Software.
    void
    create_convex_hull()
    {
        assert(facets_.size() == dimension_ + 1);
        assert(removed_facets_.empty());
        point_list outside_;
        facet_array newfacets_;
        while (!ranking_.empty()) {
            size_type const best_facet = get_best_facet();
            point_list & o_ = facets_[best_facet].outside_;
            assert(!o_.empty());
            point_iterator const apex = std::move(o_.front());
            o_.pop_front();
            if (!process_visibles(outside_, newfacets_, best_facet, apex)) {
                assert(false);
            }
            visited_.clear();
            visible_.clear();
            assert(unique_ridges_.empty());
            for (size_type const newfacet : newfacets_) {
                rank(partition(facets_[newfacet], outside_), newfacet);
            }
            newfacets_.clear();
            outside_.clear(); // to track internal_set_: splice points from outside_ to internal_set_
        }
        assert(ranking_meta_.empty());
        { // compactify
            size_type source = facets_.size();
            for (size_type const destination : removed_facets_) {
                if (destination != --source) {
                    facet & facet_ = facets_[destination];
                    facet_ = std::move(facets_.back());
                    for (size_type const neighbour : facet_.neighbours_) {
                        replace_neighbour(neighbour, source, destination);
                    }
                }
                facets_.pop_back();
            }
            facets_.shrink_to_fit();
            removed_facets_.clear();
        }
    }

    // Kurt Mehlhorn, Stefan Näher, Thomas Schilz, Stefan Schirra, Michael Seel, Raimund Seidel, and Christian Uhrig.
    // Checking geometric programs or verification of geometric structures. In Proc. 12th Annu. ACM Sympos. Comput. Geom., pages 159–165, 1996.
    bool
    check(value_type const & _eps) const
    {
        std::set< point_iterator > surface_points_;
        size_type const facets_count_ = facets_.size();
        for (size_type f = 0; f < facets_count_; ++f) { // check whether the inner point is inside relative to each hull facet or not
            facet const & facet_ = facets_[f];
            surface_points_.insert(std::cbegin(facet_.vertices_), std::cend(facet_.vertices_));
            for (size_type const neighbour : facet_.neighbours_) {
                facet const & neighbour_ = facets_[neighbour];
                if (cos_of_dihedral_angle(facet_, neighbour_) < one) { // avoiding some of roundoff errors
                    for (size_type v = 0; v < dimension_; ++v) {
                        if (neighbour_.neighbours_[v] == f) { // opposite to facet_ vertex into neighbour_ facet
                            if (_eps < facet_.distance(*neighbour_.vertices_[v])) {
                                return false; // facet is not locally convex at the ridge, common for facet_ and neighbour_ facets
                            } else {
                                break;
                            }
                        }
                    }
                }
            }
        }
        assert(!surface_points_.empty());
        vector inner_point_(zero, dimension_);
        for (point_iterator const & point_ : surface_points_) {
            auto x = std::cbegin(*point_);
            for (size_type i = 0; i < dimension_; ++i) {
                inner_point_[i] += *x;
                ++x;
            }
        }
        inner_point_ /= value_type(surface_points_.size());
        facet const & first_ = facets_.front();
        if (!(first_.distance(inner_point_) < zero)) {
            return false; // inner point is not on negative side of the first facet, therefore structure is not convex
        }
        vector ray_(zero, dimension_);
        for (point_iterator const & vertex : first_.vertices_) {
            auto x = std::cbegin(*vertex);
            for (size_type i = 0; i < dimension_; ++i) {
                ray_[i] += *x;
                ++x;
            }
        }
        ray_ /= value_type(dimension_);
        ray_ -= inner_point_;
        if (!(zero < std::inner_product(std::cbegin(ray_), std::cend(ray_), std::cbegin(first_.normal_), zero))) {
            return false;
        }
        matrix g_(dimension_); // storage (d * (d + 1)) for Gaussian elimination with partial pivoting
        for (vector & row_ : g_) {
            row_.resize(dimension_ + 1);
        }
        vector intersection_point_(zero, dimension_);
        vector centroid_(zero, dimension_);
        for (size_type f = 1; f < facets_count_; ++f) {
            using std::abs;
            facet const & facet_ = facets_[f];
            value_type const numerator_ = facet_.distance(inner_point_);
            if (!(numerator_ < zero)) {
                return false; // inner point is not on negative side of all the facets, i.e. structure is not convex
            }
            value_type const denominator_ = std::inner_product(std::cbegin(ray_), std::cend(ray_), std::cbegin(facet_.normal_), zero);
            if (!(zero < denominator_)) { // ray is parallel to the plane or directed away from the plane
                continue;
            }
            intersection_point_ = ray_;
            intersection_point_ *= -(numerator_ / denominator_);
            intersection_point_ += inner_point_;
            for (size_type v = 0; v < dimension_; ++v) {
                auto beg = std::cbegin(*facet_.vertices_[v]);
                for (size_type r = 0; r < dimension_; ++r) {
                    g_[r][v] = *beg;
                    ++beg;
                }
            }
            for (size_type r = 0; r < dimension_; ++r) {
                vector & gr_ = g_[r];
                gr_[dimension_] = intersection_point_[r];
                centroid_[r] = gr_.sum();
            }
            centroid_ /= value_type(dimension_ + 1);
            for (size_type r = 0; r < dimension_; ++r) {
                vector & gr_ = g_[r];
                gr_ -= centroid_[r];
                centroid_[r] = gr_.max() - gr_.min();
            }
            centroid_ *= centroid_;
            using std::sqrt;
            centroid_ = sqrt(centroid_.sum()) / (one + one);
            centroid_ *= facet_.normal_;
            for (size_type r = 0; r < dimension_; ++r) { // shift by half of bounding box main diagonal length in perpendicular to facet direction
                g_[r] += centroid_[r];
            }
            for (size_type i = 0; i < dimension_; ++i) { // Gaussian elimination
                vector & gi_ = g_[i];
                value_type max_ = abs(gi_[i]);
                size_type pivot = i;
                {
                    size_type p = i;
                    while (++p < dimension_) {
                        value_type y_ = abs(g_[p][i]);
                        if (max_ < y_) {
                            max_ = std::move(y_);
                            pivot = p;
                        }
                    }
                }
                assert(_eps < max_); // vertex cannot match the origin after above transformations
                if (pivot != i) {
                    gi_.swap(g_[pivot]);
                }
                value_type & gii_ = gi_[i];
                for (size_type j = i + 1; j < dimension_; ++j) {
                    vector & gj_ = g_[j];
                    value_type & gji_ = gj_[i];
                    gji_ /= gii_;
                    for (size_type k = i + 1; k <= dimension_; ++k) {
                        gj_[k] -= gji_ * gi_[k];
                    }
                    gji_ = zero;
                }
            } // g_ is upper triangular now
            bool in_range_ = true;
            {
                size_type i = dimension_;
                while (0 < i) {
                    --i;
                    vector & gi_ = g_[i];
                    value_type & xi_ = gi_[dimension_];
                    for (size_type j = i + 1; j < dimension_; ++j) {
                        xi_ -= gi_[j] * g_[j][dimension_];
                    }
                    value_type const & gii_ = gi_[i];
                    assert(_eps < abs(gii_)); // vertex cannot match the origin
                    xi_ /= gii_;
                    if ((xi_ < zero) || (one < xi_)) {
                        in_range_ = false; // barycentric coordinate does not lie in [0;1] interval => miss
                        break;
                    }
                }
            }
            if (in_range_) {
                return false; // hit
            }
        }
        return true;
    }

    bool
    check() const
    {
        return check(eps);
    }

};
