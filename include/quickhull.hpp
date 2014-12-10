#pragma once

#include <vector>
#include <deque>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <list>
#include <valarray>
#include <iterator>
#include <algorithm>
#include <utility>
#include <numeric>

#include <cmath>
#include <cassert>

template< typename points, typename point = typename points::value_type, typename value_type = typename point::value_type >
struct quick_hull
{

    using size_type = std::size_t;

    size_type const dimension_;
    points const & points_;
    value_type const eps;

    quick_hull(size_type const _dimension,
               points const & _points,
               value_type const & _eps = std::numeric_limits< value_type >::epsilon())
        : dimension_(_dimension)
        , points_(_points)
        , eps(_eps)
        , matrix_(_dimension)
        , shadow_matrix_(_dimension)
        , minor_()
    {
        assert(!(_eps < zero));
        size_type const minor_size = dimension_ - 1;
        minor_.resize(minor_size);
        for (size_type r = 0; r < minor_size; ++r) {
            matrix_[r].resize(dimension_);
            shadow_matrix_[r].resize(dimension_);
            minor_[r].resize(minor_size);
        }
        matrix_.back().resize(dimension_);
        shadow_matrix_.back().resize(dimension_);
    }

    using point_list = std::list< size_type >;
    using facet_set = std::set< size_type >;
    using facet_unordered_set = std::unordered_set< size_type >;
    using facet_deque = std::deque< size_type >;
    using point_array = std::vector< size_type >;
    using point_deque = std::deque< size_type >;

    struct facet // (d - 1)-dimensional face
    {

        using normal = std::valarray< value_type >;

        point_array vertices_; // d points (oriented)
        facet_set neighbours_;
        point_list outside_; // if empty, then is convex hull's facet, else the first point (i.e. outside_.front()) is the furthest point from this facet
        point_deque coplanar_; // coplanar points

        // hyperplane equation
        normal normal_; // components of normalized normal vector
        value_type D; // distance fromt the origin to the hyperplane

        facet(typename point_list::const_iterator _first,
              typename point_list::const_iterator _middle,
              typename point_list::const_iterator _last)
            : vertices_(_first, std::prev(_middle))
        {
            vertices_.insert(std::cend(vertices_), _middle, _last);
            normal_.resize(vertices_.size());
        }

        facet(point_array && _vertices,
              size_type const _neighbour)
            : vertices_(std::move(_vertices))
            , neighbours_()
            , normal_(vertices_.size())
        {
            neighbours_.insert(_neighbour);
        }

        value_type
        distance(point const & _point) const
        {
            return std::inner_product(std::cbegin(normal_), std::cend(normal_), std::cbegin(_point), D);
        }

    };

    using facets_storage = std::deque< facet >;

    facets_storage facets_;
    point_list internal_set_;

    value_type
    cos_of_dihedral_angle(facet const & _this, facet const & _other) const // for faces merging in the future
    {
        return std::inner_product(std::cbegin(_this.normal_), std::cend(_this.normal_), std::cbegin(_other.normal_), zero);
    }

private :

    // math (simple functions, matrices, etc):

    value_type const zero = value_type(0);
    value_type const one = value_type(1);

    using row = std::valarray< value_type >;
    using matrix = std::vector< row >;

    matrix matrix_;
    matrix shadow_matrix_;
    matrix minor_;

    void
    transpose()
    { // to cheaper filling columns with ones
        for (size_type r = 0; r < dimension_; ++r) {
            row & row_ = shadow_matrix_[r];
            for (size_type c = 1 + r; c < dimension_; ++c) {
                using std::swap;
                swap(shadow_matrix_[c][r], row_[c]);
            }
        }
    }

    void
    restore_matrix()
    {
        for (size_type r = 0; r < dimension_; ++r) {
            matrix_[r] = shadow_matrix_[r];
        }
    }

    void
    restore_matrix(size_type const _identity) // load matrix from storage with replacing of _identity column with ones
    {
        for (size_type r = 0; r < dimension_; ++r) { // col
            row & row_ = matrix_[r];
            if (r == _identity) {
                row_ = one;
            } else {
                row_ = shadow_matrix_[r];
            }
        }
    }

    void
    square_matrix(size_type const _size) // matrix_ = shadow_matrix_ * transposed shadow_matrix_
    {
        assert(_size < dimension_);
        for (size_type r = 0; r < _size; ++r) {
            row & lhs_ = shadow_matrix_[r];
            row const & row_ = matrix_[r];
            auto const rbeg = std::cbegin(row_);
            auto const rend = std::cend(row_);
            for (size_type c = 0; c < _size; ++c) {
                lhs_[c] = std::inner_product(rbeg, rend, std::cbegin(matrix_[c]), zero);
            }
        }
    }

    value_type
    det(matrix & _matrix, size_type const _size) // based on LU factorization
    {
        value_type det_ = one;
        for (size_type i = 0; i < _size; ++i) {
            size_type p = i;
            using std::abs;
            value_type max_ = abs(_matrix[p][i]);
            size_type pivot = p;
            while (++p < _size) {
                value_type y_ = abs(_matrix[p][i]);
                if (max_ < y_) {
                    max_ = std::move(y_);
                    pivot = p;
                }
            }
            if (!(eps < max_)) { // regular?
                return zero; // singular
            }
            row & ri_ = _matrix[i];
            if (pivot != i) {
                det_ = -det_; // each permutation flips sign of det
                ri_.swap(_matrix[pivot]);
            }
            value_type & dia_ = ri_[i];
            value_type const inv_ = one / dia_;
            det_ *= std::move(dia_); // det is multiple of diagonal elements
            for (size_type j = 1 + i; j < _size; ++j) {
                _matrix[j][i] *= inv_;
            }
            for (size_type a = 1 + i; a < _size; ++a) {
                row & a_ = minor_[a - 1];
                value_type const & ai_ = _matrix[a][i];
                for (size_type b = 1 + i; b < _size; ++ b) {
                    a_[b - 1] = ai_ * ri_[b];
                }
            }
            for (size_type a = 1 + i; a < _size; ++a) {
                row const & a_ = minor_[a - 1];
                row & ra_ = _matrix[a];
                for (size_type b = 1 + i; b < _size; ++ b) {
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

    // geometry and basic operations on geometric primitives:

    using point_set = std::set< size_type >;
    using ordered_vertices = std::deque< point_set >;
    using facet_set_desc = std::set< size_type, std::greater< size_type > >;

    ordered_vertices ordered_; // ordered, but not oriented vertices of facets
    facet_set_desc removed_facets_;

    void
    set_hyperplane_equation(facet & _facet)
    {
        for (size_type r = 0; r < dimension_; ++r) {
            std::copy_n(std::cbegin(points_[_facet.vertices_[r]]), dimension_, std::begin(shadow_matrix_[r]));
        }
        transpose();
        value_type N = zero;
        for (size_type i = 0; i < dimension_; ++i) {
            restore_matrix(i);
            value_type n = det();
            N += n * n;
            _facet.normal_[i] = std::move(n);
        }
        using std::sqrt;
        N = sqrt(N);
        restore_matrix();
        _facet.normal_ /= N;
        _facet.D = -det() / std::move(N);
    }

    size_type
    add_facet(point_array && _vertices, size_type const _neighbour)
    {
        assert(ordered_.size() == facets_.size());
        if (removed_facets_.empty()) {
            size_type const f = facets_.size();
            facets_.emplace_back(std::move(_vertices), _neighbour);
            facet & facet_ = facets_.back();
            set_hyperplane_equation(facet_);
            ordered_.emplace_back(std::cbegin(facet_.vertices_),
                                  std::cend(facet_.vertices_));
            return f;
        } else {
            auto const rend = std::prev(std::cend(removed_facets_));
            size_type const f = *rend;
            removed_facets_.erase(rend);
            facet & facet_ = facets_[f];
            facet_.vertices_ = std::move(_vertices);
            assert(facet_.neighbours_.empty());
            facet_.neighbours_.insert(_neighbour);
            set_hyperplane_equation(facet_);
            ordered_[f].insert(std::cbegin(facet_.vertices_),
                               std::cend(facet_.vertices_));
            return f;
        }
    }

    // http://math.stackexchange.com/questions/822741/
    value_type
    hypervolume(point_list const & _vertices, point const & _apex)
    {
        size_type const rows_count = _vertices.size();
        assert(!(dimension_ < rows_count));
        row & origin_ = shadow_matrix_.back();
        std::copy_n(std::cbegin(_apex), dimension_, std::begin(origin_));
        auto vertex = std::cbegin(_vertices);
        for (size_type r = 0; r < rows_count; ++r) {
            std::copy_n(std::cbegin(points_[*vertex]), dimension_, std::begin(matrix_[r]));
            ++vertex;
        }
        for (size_type r = 0; r < rows_count; ++r) { // vectorize
            matrix_[r] -= origin_;
        }
        if (rows_count == dimension_) { // oriented hypervolume
            return det();
        } else { // non-oriented rows_count-dimensional measure
            square_matrix(rows_count);
            using std::sqrt;
            return sqrt(det(shadow_matrix_, rows_count));
        }
    }

    value_type
    hypervolume(point_list const & _vertices, size_type const _apex)
    {
        return hypervolume(_vertices, points_[_apex]);
    }

    value_type
    steal_best(point_list & _from, point_list & _to)
    {
        using std::abs;
        value_type hypervolume_ = zero;
        auto const end = std::cend(_from);
        auto furthest = end;
        for (auto it = std::cbegin(_from); it != end; ++it) {
            value_type v_ = hypervolume(_to, *it);
            if (abs(hypervolume_) < abs(v_)) {
                hypervolume_ = std::move(v_);
                furthest = it;
            }
        }
        if (furthest != end) {
            _to.splice(std::cend(_to), _from, furthest);
        }
        return hypervolume_;
    }

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
    }

    size_type
    get_best_facet() const
    {
        if (ranking_.empty()) {
            return facets_.size(); // special value
        } else {
            return std::prev(std::cend(ranking_))->second;
        }
    }

    value_type
    partition(facet & _facet, point_list & _points)
    {
        auto it = std::cbegin(_points);
        auto const end = std::cend(_points);
        value_type distance_ = zero;
        while (it != end) {
            auto const next = std::next(it);
            size_type const p = *it;
            value_type d_ = _facet.distance(points_[p]);
            if (eps < d_) {
                if ((distance_ < d_) || _facet.outside_.empty()) {
                    distance_ = std::move(d_);
                    _facet.outside_.splice(std::cbegin(_facet.outside_), _points, it);
                } else {
                    _facet.outside_.splice(std::cend(_facet.outside_), _points, it);
                }
            } else if (!(d_ < -eps)) { // coplanar
                _facet.coplanar_.push_back(p);
            }
            it = next;
        }
        return distance_;
    }

    void
    adjacency(facet_deque const & _newfacets)
    {
        auto const nend = std::cend(_newfacets);
        for (auto first = std::cbegin(_newfacets); first != nend; ++first) {
            size_type const f = *first;
            facet & first_facet_ = facets_[f];
            size_type neighbours_count = first_facet_.neighbours_.size();
            if (neighbours_count < dimension_) {
                point_set & first_ = ordered_[f];
                auto const lbeg = std::cbegin(first_);
                auto const lend = std::cend(first_);
                for (auto second = std::next(first); second != nend; ++second) {
                    size_type const s = *second;
                    point_set & second_ = ordered_[s];
                    auto const rend = std::cend(second_);
                    auto r = std::cbegin(second_);
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
                            facets_[s].neighbours_.insert(f);
                            if (!(++neighbours_count < dimension_)) {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

public : // largest possible simplex heuristic, convex hull algorithm

    point_list
    create_simplex()
    {
        assert(1 < dimension_);
        size_type const input_size = points_.size();
        assert(0 < input_size);
        for (size_type i = 0; i < input_size; ++i) {
            internal_set_.push_back(i);
        }
        point_list basis_;
        basis_.splice(std::cend(basis_), internal_set_, std::cbegin(internal_set_));
        steal_best(internal_set_, basis_);
        if (basis_.size() != 2) {
            return basis_; // can't find linearly independent second point
        }
        internal_set_.splice(std::cend(internal_set_), basis_, std::cbegin(basis_)); // rejudge 0-indexed point
        for (size_type i = 1; i < dimension_; ++i) {
            steal_best(internal_set_, basis_);
            if (basis_.size() != i + 1) {
                return basis_; // can't find (i + 1) linearly independent point
            }
        }
        value_type hypervolume_ = steal_best(internal_set_, basis_);
        if (basis_.size() != dimension_ + 1) {
            return basis_; // can't find linearly independent (d + 1) point
        }
        // simplex construction
        bool inward_ = (zero < hypervolume_); // is top oriented?
        auto const vbeg = std::cbegin(basis_);
        auto const vend = std::cend(basis_);
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) {
            size_type const newfacet = facets_.size();
            facets_.emplace_back(vbeg, exclusive, vend);
            facet & facet_ = facets_.back();
            inward_ = !inward_;
            if (inward_) {
                std::swap(facet_.vertices_.front(),
                          facet_.vertices_.back());
            }
            set_hyperplane_equation(facet_);
            ordered_.emplace_back(std::cbegin(facet_.vertices_), std::cend(facet_.vertices_));
            rank(partition(facet_, internal_set_), newfacet);
        }
        { // adjacency
            for (size_type i = 0; i < dimension_; ++i) {
                facet_set & neighbours_ = facets_[i].neighbours_;
                for (size_type j = 1 + i; j < dimension_ + 1; ++j) {
                    neighbours_.insert(j);
                    facets_[j].neighbours_.insert(i);
                }
            }
        }
        return basis_;
    }

    void
    create_convex_hull()
    {
        assert(facets_.size() == dimension_ + 1);
        assert(ordered_.size() == dimension_ + 1);
        assert(removed_facets_.empty());
        point_list outside_;
        facet_set visited_; // makes finding of visible facets in pravo output-sensitive
        auto const vend = std::end(visited_);
        facet_deque viewable_;
        facet_unordered_set visible_facets_;
        auto const vfend = std::end(visible_facets_); // not invalidated till the end
        facet_set neighbours_;
        facet_deque newfacets_;
        point_array vertices_;
        point_array ridge_; // horizon ridge + furthest point = new facet
        for (size_type best_facet = get_best_facet(); best_facet != facets_.size(); best_facet = get_best_facet()) {
            facet & best_facet_ = facets_[best_facet];
            assert(!best_facet_.outside_.empty());
            size_type const apex = best_facet_.outside_.front();
            best_facet_.outside_.pop_front();
            point const & apex_ = points_[apex];
            assert(visible_facets_.empty());
            visible_facets_.insert(best_facet);
            { // find visible facets
                assert(visited_.empty());
                assert(viewable_.empty());
                viewable_.assign(std::cbegin(best_facet_.neighbours_), std::cend(best_facet_.neighbours_));
                while (!viewable_.empty()) {
                    size_type const inner = viewable_.front();
                    viewable_.pop_front();
                    facet const & watchable_ = facets_[inner];
                    if (eps < watchable_.distance(apex_)) {
                        visible_facets_.insert(inner);
                        for (size_type const neighbour : watchable_.neighbours_) {
                            if (visited_.find(neighbour) != vend) {
                                viewable_.push_back(neighbour);
                            }
                        }
                    }
                    visited_.insert(inner);
                }
                visited_.clear();
            }
            assert(newfacets_.empty());
            for (size_type const visible_facet : visible_facets_) {
                facet & visible_facet_ = facets_[visible_facet];
                vertices_ = std::move(visible_facet_.vertices_);
                outside_.splice(outside_.cend(), std::move(visible_facet_.outside_));
                neighbours_ = std::move(visible_facet_.neighbours_);
                unrank(visible_facet);
                { // remove facet
                    removed_facets_.insert(visible_facet);
                    visible_facet_.coplanar_.clear();
                    ordered_[visible_facet].clear();
                }
                for (size_type const neighbour : neighbours_) {
                    if (visible_facets_.find(neighbour) == vfend) {
                        facet & horizon_facet_ = facets_[neighbour];
                        point_set const & horizon_ = ordered_[neighbour];
                        {
                            assert(ridge_.empty());
                            ridge_.reserve(dimension_);
                            auto const hend = std::cend(horizon_);
                            for (size_type const vertex : vertices_) { // facets intersection with keeping of points order as it is in visible facet
                                if (horizon_.find(vertex) == hend) {
                                    ridge_.push_back(apex);
                                } else {
                                    ridge_.push_back(vertex);
                                }
                            }
                            assert(ridge_.size() == dimension_); // facet
                        }
                        size_type const newfacet = add_facet(std::move(ridge_), neighbour);
                        newfacets_.push_back(newfacet);
                        facet_set & horizon_neighbours_ = horizon_facet_.neighbours_;
                        horizon_neighbours_.erase(visible_facet);
                        horizon_neighbours_.insert(newfacet);
                    }
                }
            }
            visible_facets_.clear();
            adjacency(newfacets_);
            for (size_type const newfacet : newfacets_) {
                rank(partition(facets_[newfacet], outside_), newfacet);
            }
            newfacets_.clear();
            internal_set_.splice(std::cend(internal_set_), std::move(outside_));
        }
        assert(ranking_.empty());
        assert(ranking_meta_.empty());
        { // compactify
            size_type source = facets_.size();
            for (size_type const destination : removed_facets_) {
                if (destination != --source) {
                    facet & facet_ = facets_[destination];
                    facet_ = std::move(facets_.back());
                    for (size_type const neighbour : facet_.neighbours_) {
                        facet_set & neighbours_ = facets_[neighbour].neighbours_;
                        neighbours_.erase(source);
                        neighbours_.insert(destination);
                    }
                }
                facets_.pop_back();
            }
            facets_.shrink_to_fit();
            removed_facets_.clear();
        }
        ordered_.clear();
        ordered_.shrink_to_fit();
    }

};
