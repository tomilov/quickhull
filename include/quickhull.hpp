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

#include <cassert>

template< typename points_type >
struct quick_hull
{

    using size_type = std::size_t;

    size_type const dimension_;

    using point_type = typename points_type::value_type;
    using G = typename point_type::value_type;

private : // math (simple functions, matrices, etc)

    G const zero = G(0);
    G const one = G(1);

    G
    abs(G const & _x) const
    {
        return (_x < zero) ? -_x : _x;
    }

    G
    abs(G && _x) const
    {
        return (_x < zero) ? -std::move(_x) : std::move(_x);
    }

    using row_type = std::valarray< G >;
    using matrix_type = std::vector< row_type >;

    matrix_type matrix_;
    matrix_type shadow_matrix_;
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
        for (size_type row = 0; row < dimension_; ++row) { // col
            row_type & row_ = matrix_[row];
            if (row == _identity) {
                row_ = one;
            } else {
                row_ = shadow_matrix_[row];
            }
        }
    }

    void
    square_matrix(size_type const _size) // matrix_ = shadow_matrix_ * transposed shadow_matrix_
    {
        assert(_size < dimension_);
        for (size_type row = 0; row < _size; ++row) {
            row_type & lhs_ = shadow_matrix_[row];
            row_type const & row_ = matrix_[row];
            auto const rbeg = std::cbegin(row_);
            auto const rend = std::cend(row_);
            for (size_type col = 0; col < _size; ++col) {
                lhs_[col] = std::inner_product(rbeg, rend, std::cbegin(matrix_[col]), zero);
            }
        }
    }

    G
    det(matrix_type & _matrix, size_type const _size) // based on LU factorization
    {
        G det_ = one;
        for (size_type i = 0; i < _size; ++i) {
            size_type p_ = i;
            G max_ = abs(_matrix[p_][i]);
            size_type pivot_ = p_;
            while (++p_ < _size) {
                G y_ = abs(_matrix[p_][i]);
                if (max_ < y_) {
                    max_ = std::move(y_);
                    pivot_ = p_;
                }
            }
            if (!(eps < max_)) { // regular?
                return zero; // singular
            }
            row_type & ri_ = _matrix[i];
            if (pivot_ != i) {
                det_ = -det_; // each permutation flips sign of det
                ri_.swap(_matrix[pivot_]);
            }
            G & dia_ = ri_[i];
            G const inv_ = one / dia_;
            det_ *= std::move(dia_); // det is multiple of diagonal elements
            for (size_type j = 1 + i; j < _size; ++j) {
                _matrix[j][i] *= inv_;
            }
            for (size_type a = 1 + i; a < _size; ++a) {
                row_type & a_ = minor_[a - 1];
                G const & ai_ = _matrix[a][i];
                for (size_type b = 1 + i; b < _size; ++ b) {
                    a_[b - 1] = ai_ * ri_[b];
                }
            }
            for (size_type a = 1 + i; a < _size; ++a) {
                row_type const & a_ = minor_[a - 1];
                row_type & ra_ = _matrix[a];
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
        return det(matrix_, dimension_);
    }

public :

    using point_array = std::vector< size_type >;
    using point_deque = std::deque< size_type >;
    using point_list = std::list< size_type >;
    using point_set = std::set< size_type >;

    points_type const & points_;
    point_list internal_set_;
    G const eps;

    quick_hull(size_type const _dimension, points_type const & _points,
               G const & _eps = std::numeric_limits< G >::epsilon())
        : dimension_(_dimension)
        , matrix_(_dimension)
        , shadow_matrix_(_dimension)
        , minor_(_dimension - 1)
        , points_(_points)
        , eps(_eps)
    {
        assert(!(_eps < zero));
        size_type const minor_size_ = dimension_ - 1;
        for (size_type row = 0; row < minor_size_; ++row) {
            matrix_[row].resize(dimension_);
            shadow_matrix_[row].resize(dimension_);
            minor_[row].resize(minor_size_);
        }
        matrix_.back().resize(dimension_);
        shadow_matrix_.back().resize(dimension_);
    }

    using facet_set = std::set< size_type >;
    using facet_unordered_set = std::unordered_set< size_type >;

    struct facet // (d - 1)-dimensional face
    {

        using normal_type = std::valarray< G >;

        point_array vertices_; // d points : oriented
        facet_set neighbours_;
        point_deque outside_; // if empty, then is convex hull's facet, else the first point (i.e. outside_.front()) is the furthest point from this facet
        point_deque coplanar_; // coplanar points

        // hyperplane equation
        normal_type normal_; // components of normalized normal vector
        G D; // distance fromt the origin to the hyperplane

        facet(typename point_list::const_iterator first,
              typename point_list::const_iterator middle,
              typename point_list::const_iterator last)
            : vertices_(first, std::prev(middle))
        {
            vertices_.insert(std::cend(vertices_), middle, last);
            normal_.resize(vertices_.size());
        }

        facet(point_array && _vertices,
              size_type const _neighbour)
            : vertices_(std::move(_vertices))
            , neighbours_({_neighbour})
            , normal_(vertices_.size())
        { ; }

        G
        distance(point_type const & _point) const
        {
            return std::inner_product(std::cbegin(normal_), std::cend(normal_), std::cbegin(_point), D);
        }

    };

    G
    cos_of_dihedral_angle(facet const & _this, facet const & _other) const // for faces merging in the future
    {
        return std::inner_product(std::cbegin(_this.normal_), std::cend(_this.normal_), std::cbegin(_other.normal_), zero);
    }

    using facets_storage_type = std::deque< facet >;

    facets_storage_type facets_;

private : // geometry and basic operations on geometric primitives

    using vertices_sets_type = std::deque< point_set >;
    using facets_type = std::deque< size_type >;
    using facet_set_desc = std::set< size_type, std::greater< size_type > >;

    vertices_sets_type ordered_; // ordered, but not oriented vertices of facets
    facet_set_desc removed_facets_;

    void
    set_hyperplane_equation(facet & _facet)
    {
        for (size_type row = 0; row < dimension_; ++row) {
            std::copy_n(std::cbegin(points_[_facet.vertices_[row]]), dimension_, std::begin(shadow_matrix_[row]));
        }
        transpose();
        G N = zero;
        for (size_type i = 0; i < dimension_; ++i) {
            restore_matrix(i);
            G n = det();
            N += n * n;
            _facet.normal_[i] = std::move(n);
        }
        using std::sqrt;
        N = one / sqrt(N);
        restore_matrix();
        _facet.D = -det() * N;
        _facet.normal_ *= std::move(N);
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
            facet_.neighbours_ = {_neighbour};
            set_hyperplane_equation(facet_);
            ordered_[f].insert(std::cbegin(facet_.vertices_),
                               std::cend(facet_.vertices_));
            return f;
        }
    }

    // http://math.stackexchange.com/questions/822741/
    G
    hypervolume(point_list const & _vertices, point_type const & _apex)
    {
        size_type const rows_ = _vertices.size();
        assert(!(dimension_ < rows_));
        row_type & origin_ = shadow_matrix_.back();
        std::copy_n(std::cbegin(_apex), dimension_, std::begin(origin_));
        auto vertex = std::cbegin(_vertices);
        for (size_type row = 0; row < rows_; ++row) {
            assert(vertex != std::cend(_vertices));
            std::copy_n(std::cbegin(points_[*vertex]), dimension_, std::begin(matrix_[row]));
            ++vertex;
        }
        for (size_type row = 0; row < rows_; ++row) { // vectorize
            matrix_[row] -= origin_;
        }
        if (rows_ == dimension_) { // oriented hypervolume
            return det();
        } else { // non-oriented rows_-dimensional measure
            square_matrix(rows_);
            using std::sqrt;
            return sqrt(det(shadow_matrix_, rows_));
        }
    }

    G
    hypervolume(point_list const & _vertices, size_type const _apex)
    {
        return hypervolume(_vertices, points_[_apex]);
    }

    G
    steal_best(point_list & _from, point_list & _to)
    {
        auto it = std::cbegin(_from);
        auto const end = std::cend(_from);
        G hypervolume_ = hypervolume(_to, *it);
        auto furthest = it;
        while (++it != end) {
            G v_ = hypervolume(_to, *it);
            if (abs(hypervolume_) < abs(v_)) {
                hypervolume_ = std::move(v_);
                furthest = it;
            }
        }
        if (eps < abs(hypervolume_)) {
            _to.splice(std::cend(_to), _from, furthest);
        }
        return hypervolume_;
    }

    using ranking_type = std::multimap< G, size_type >;
    using ranking_meta_type = std::unordered_map< size_type, typename ranking_type::iterator >;

    ranking_type ranking_;
    ranking_meta_type ranking_meta_;

    void
    rank(G && _orientation, size_type const _facet)
    {
        if (eps < _orientation) {
            auto r = ranking_.emplace(std::move(_orientation), _facet);
            ranking_meta_.emplace(_facet, std::move(r));
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
            auto const r = std::prev(std::cend(ranking_));
            return r->second;
        }
    }

    G
    partition(facet & _facet, point_list & _points)
    {
        auto it = std::cbegin(_points);
        auto const end = std::cend(_points);
        G distance_ = zero;
        while (it != end) {
            auto const next = std::next(it);
            size_type const p = *it;
            G d_ = _facet.distance(points_[p]);
            if (eps < d_) {
                if ((distance_ < d_) || _facet.outside_.empty()) {
                    distance_ = std::move(d_);
                    _facet.outside_.push_front(p);
                } else {
                    _facet.outside_.push_back(p);
                }
                _points.erase(it);
            } else if (!(d_ < -eps)) { // coplanar
                _facet.coplanar_.push_back(p);
            }
            it = next;
        }
        return abs(std::move(distance_));
    }

    void
    adjacency(facets_type const & _newfacets)
    {
        auto const nend = std::cend(_newfacets);
        for (auto first = std::cbegin(_newfacets); first != nend; ++first) {
            size_type const f = *first;
            point_set & first_ = ordered_[f];
            auto const lbeg = std::cbegin(first_);
            auto const lend = std::cend(first_);
            facet & first_facet_ = facets_[f];
            size_type neighbours_count_ = first_facet_.neighbours_.size();
            if (neighbours_count_ < dimension_) {
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
                            if (!(++neighbours_count_ < dimension_)) {
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
        size_type const size_ = points_.size();
        assert(dimension_ < size_);
        internal_set_.resize(size_);
        std::iota(std::begin(internal_set_), std::end(internal_set_), size_type(0));
        point_list basis_;
        basis_.splice(std::cend(basis_), internal_set_, std::cbegin(internal_set_));
        if (!(eps < abs(steal_best(internal_set_, basis_)))) {
            return basis_; // can't find linearly independent point
        }
        internal_set_.splice(std::cend(internal_set_), basis_, std::cbegin(basis_)); // rejudge 0-indexed point
        for (size_type i = 1; i < dimension_; ++i) {
            if (!(eps < abs(steal_best(internal_set_, basis_)))) {
                return basis_; // can't find linearly independent point
            }
        }
        assert(basis_.size() == dimension_); // facet
        G hypervolume_ = steal_best(internal_set_, basis_);
        if (!(eps < abs(hypervolume_))) {
            return basis_; // can't find linearly independent point
        }
        assert(basis_.size() == dimension_ + 1); // simplex
        bool inward_ = (zero < hypervolume_); // is top oriented?
        auto const vbeg = std::cbegin(basis_);
        auto const vend = std::cend(basis_);
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) {
            size_type const newfacet = facets_.size();
            facets_.emplace_back(vbeg, exclusive, vend);
            facet & facet_ = facets_.back();
            inward_ = !inward_;
            if (inward_) {
                using std::swap;
                swap(facet_.vertices_.front(),
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
        point_list outside_;
        facet_set visited_;
        facet_unordered_set viewable_;
        facet_set visible_facets_;
        auto const vfend = std::end(visible_facets_); // not invalidated till the end
        facet_set neighbours_;
        facets_type newfacets_;
        point_array vertices_;
        point_array ridge_; // horizon ridge + furthest point = new facet
        for (size_type best_facet = get_best_facet(); best_facet != facets_.size(); best_facet = get_best_facet()) {
            facet & best_facet_ = facets_[best_facet];
            size_type const apex = best_facet_.outside_.front();
            best_facet_.outside_.pop_front();
            point_type const & apex_ = points_[apex];
            visible_facets_ = {best_facet};
            { // find visible facets
                assert(visited_.empty());
                assert(viewable_.empty());
                viewable_.insert(std::cbegin(best_facet_.neighbours_),
                                 std::cend(best_facet_.neighbours_));
                while (!viewable_.empty()) {
                    auto const first = std::cbegin(viewable_);
                    size_type const f = *first;
                    facet const & watchable_ = facets_[f];
                    if (eps < watchable_.distance(apex_)) {
                        visible_facets_.insert(f);
                        std::set_difference(std::cbegin(watchable_.neighbours_),
                                            std::cend(watchable_.neighbours_),
                                            std::cbegin(visited_), std::cend(visited_),
                                            std::inserter(viewable_, std::end(viewable_)));
                    }
                    visited_.insert(f);
                    viewable_.erase(first);
                }
                visited_.clear();
            }
            assert(newfacets_.empty());
            for (size_type const visible_facet : visible_facets_) {
                facet & visible_facet_ = facets_[visible_facet];
                vertices_ = std::move(visible_facet_.vertices_);
                outside_.insert(std::cend(outside_), std::cbegin(visible_facet_.outside_), std::cend(visible_facet_.outside_));
                visible_facet_.outside_.clear();
                neighbours_ = std::move(visible_facet_.neighbours_);
                unrank(visible_facet);
                removed_facets_.insert(visible_facet);
                visible_facet_.coplanar_.clear();
                ordered_[visible_facet].clear();
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
                        horizon_neighbours_.insert(std::cend(horizon_neighbours_), newfacet);
                    }
                }
            }
            adjacency(newfacets_);
            for (size_type const newfacet : newfacets_) {
                rank(partition(facets_[newfacet], outside_), newfacet);
            }
            newfacets_.clear();
            internal_set_.splice(std::cend(internal_set_), std::move(outside_));
        }
        ordered_.clear();
        ordered_.shrink_to_fit();
        {
            size_type source_ = facets_.size();
            for (size_type const destination_ : removed_facets_) {
                if (destination_ != --source_) {
                    facet & facet_ = facets_[destination_];
                    facet_ = std::move(facets_.back());
                    for (size_type const neighbour : facet_.neighbours_) {
                        facet_set & neighbours_ = facets_[neighbour].neighbours_;
                        neighbours_.erase(source_);
                        neighbours_.insert(destination_);
                    }
                }
                facets_.pop_back();
            }
            facets_.shrink_to_fit();
        }
    }

};
