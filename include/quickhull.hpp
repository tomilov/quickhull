#pragma once

#include <valarray>
#include <vector>
#include <list>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <utility>
#include <numeric>

#include <cmath>
#include <cassert>

template< typename points,
          typename point = typename points::value_type,
          typename value_type = typename point::value_type >
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
        , origin_(zero, _dimension)
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

    using point_array = std::vector< size_type >;
    using point_list = std::list< size_type >;
    using facet_array = std::vector< size_type >;

    struct facet // (d - 1)-dimensional face
    {

        using normal = std::valarray< value_type >;

        point_array vertices_; // d points (oriented)
        point_list outside_; // if empty, then is convex hull's facet, else the first point (i.e. outside_.front()) is the furthest point from this facet
        facet_array neighbours_; // neighbouring facets

        // hyperplane equation
        normal normal_; // components of normalized normal vector
        value_type D; // distance fromt the origin to the hyperplane

        void
        init(point_array && _vertices,
             size_type const _neighbour)
        {
            vertices_ = std::move(_vertices);
            size_type const dimension_ = vertices_.size();
            neighbours_.reserve(dimension_);
            neighbours_.push_back(_neighbour);
            normal_.resize(dimension_);
        }

        facet(point_array && _vertices,
              size_type const _neighbour)
        {
            init(std::move(_vertices), _neighbour);
        }

        facet(typename point_list::const_iterator _first,
              typename point_list::const_iterator _middle,
              typename point_list::const_iterator _last)
            : vertices_(_first, std::prev(_middle))
        {
            vertices_.insert(std::cend(vertices_), _middle, _last);
            size_type const dimension_ = vertices_.size();
            neighbours_.reserve(dimension_);
            normal_.resize(dimension_);
        }

        value_type
        distance(point const & _point) const
        {
            return std::inner_product(std::cbegin(normal_), std::cend(normal_), std::cbegin(_point), D);
        }

    };

    using facets_storage = std::vector< facet >;

    facets_storage facets_;

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
    row origin_;

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
    det(matrix & _matrix, size_type const _size) // based on LUP decomposition (complexity is 2 * n^3 / 3 + O(n^2) vs 4 * n^3 / 3 + O(n^2) for QR decomposition via Householder reflections)
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

    std::vector< point_array > ordered_; // ordered, but not oriented vertices of facets
    std::set< size_type, std::greater< size_type > > removed_facets_;

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
            point_array & ordered_vertices_ = ordered_.back();
            std::sort(std::begin(ordered_vertices_), std::end(ordered_vertices_));
            return f;
        } else {
            auto const rend = std::prev(std::cend(removed_facets_));
            size_type const f = *rend;
            removed_facets_.erase(rend);
            facet & facet_ = facets_[f];
            facet_.init(std::move(_vertices), _neighbour);
            set_hyperplane_equation(facet_);
            point_array & ordered_vertices_ = ordered_[f];
            assert(ordered_vertices_.empty());
            ordered_vertices_.assign(std::cbegin(facet_.vertices_),
                                     std::cend(facet_.vertices_));
            std::sort(std::begin(ordered_vertices_), std::end(ordered_vertices_));
            return f;
        }
    }

    // http://math.stackexchange.com/questions/822741/
    value_type
    hypervolume(point_list const & _vertices, point const & _apex)
    {
        size_type const rows_count = _vertices.size();
        assert(!(dimension_ < rows_count));
        std::copy_n(std::cbegin(_apex), dimension_, std::begin(origin_));
        auto vertex = std::cbegin(_vertices);
        for (size_type r = 0; r < rows_count; ++r) { // affine space -> vector space
            row & row_ = matrix_[r];
            std::copy_n(std::cbegin(points_[*vertex]), dimension_, std::begin(row_));
            row_ -= origin_;
            ++vertex;
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

    void
    orthogonalize(point_list const & _affine_space, size_type const _rank)
    {
        std::copy_n(std::cbegin(points_[_affine_space.back()]), dimension_, std::begin(origin_));
        auto vertex = std::begin(_affine_space);
        for (size_type r = 0; r < _rank; ++r) { // affine space -> vector space
            row & row_ = shadow_matrix_[r];
            std::copy_n(std::cbegin(points_[*vertex]), dimension_, std::begin(row_));
            row_ -= origin_;
            ++vertex;
        }
        for (size_type i = 0; i < _rank; ++i) { // Householder transformation // shadow_matrix_ is packed QR
            value_type norm_ = zero;
            row & qri_ = shadow_matrix_[i];
            for (size_type j = i; j < dimension_; ++j) {
                value_type const & qrij_ = qri_[j];
                norm_ += qrij_ * qrij_;
            }
            using std::sqrt;
            norm_ = sqrt(norm_);
            assert(eps < norm_);
            value_type & qrii_ = qri_[i];
            bool const sign_ = (zero < qrii_);
            value_type factor_ = norm_ * (norm_ + (sign_ ? qrii_ : -qrii_));
            assert(eps < factor_);
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
                row & qrj_ = shadow_matrix_[j];
                value_type s_ = zero;
                for (size_type k = i; k < dimension_; ++k) {
                    s_ += qri_[k] * qrj_[k];
                }
                for (size_type k = i; k < dimension_; ++k) {
                    qrj_[k] -= qri_[k] * s_;
                }
            }
        }
        for (size_type i = 0; i < _rank; ++i) { // matrix_ is Q
            row & qi_ = matrix_[i];
            qi_ = zero;
            qi_[i] = one;
            size_type j = _rank;
            while (0 < j) {
                --j;
                row & qrj_ = shadow_matrix_[j];
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

    value_type
    distance(size_type const _apex, size_type const _rank)
    {
        row & apex_ = shadow_matrix_.front();
        auto const abeg = std::begin(apex_);
        auto const aend = std::end(apex_);
        std::copy_n(std::cbegin(points_[_apex]), dimension_, abeg);
        apex_ -= origin_;
        row & projection_ = shadow_matrix_.back();
        projection_ = apex_;
        for (size_type i = 0; i < _rank; ++i) {
            row const & qi_ = matrix_[i];
            projection_ -= std::inner_product(abeg, aend, std::cbegin(qi_), zero) * qi_;
        }
        projection_ *= projection_;
        using std::sqrt;
        return sqrt(projection_.sum());
    }

    value_type
    steal_best(point_list & _from, point_list & _to)
    {
        assert(!_from.empty());
        size_type const rank_ = _to.size() - 1;
        orthogonalize(_to, rank_);
        value_type distance_ = zero;
        auto const end = std::cend(_from);
        auto furthest = end;
        for (auto it = std::cbegin(_from); it != end; ++it) {
            value_type d_ = distance(*it, rank_);
            if (distance_ < d_) {
                distance_ = std::move(d_);
                furthest = it;
            }
        }
        value_type hypervolume_ = hypervolume(_to, *furthest);
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
    unrank_and_remove(size_type const _facet)
    {
        auto const r = ranking_meta_.find(_facet);
        if (r != std::end(ranking_meta_)) {
            ranking_.erase(r->second);
            ranking_meta_.erase(r);
        }
        removed_facets_.insert(_facet);
        ordered_[_facet].clear();
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
            }
            it = next;
        }
        return distance_;
    }

    void
    adjacency(facet_array const & _newfacets)
    {
        size_type const newfacets_count = _newfacets.size();
        for (size_type i = 0; i < newfacets_count; ++i) {
            size_type const f = _newfacets[i];
            point_array const & first_ = ordered_[f];
            facet & first_facet_ = facets_[f];
            size_type neighbours_count = first_facet_.neighbours_.size();
            if (neighbours_count < dimension_) {
                for (size_type j = i + 1; j < newfacets_count; ++j) {
                    size_type const s = _newfacets[j];
                    point_array const & second_ = ordered_[s];
                    size_type r = 0;
                    size_type l = 0;
                    bool lgood = false;
                    bool rgood = false;
                    while ((l != dimension_) && (r != dimension_)) {
                        size_type const left = first_[l];
                        size_type const right = second_[r];
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
                    if (lgood != ((l != dimension_) && (++l == dimension_))) {
                        if (rgood != ((r != dimension_) && (++r == dimension_))) {
                            first_facet_.neighbours_.push_back(s);
                            assert(facets_[s].neighbours_.size() < dimension_);
                            facets_[s].neighbours_.push_back(f);
                            if (!(++neighbours_count < dimension_)) {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    using facet_unordered_set = std::unordered_set< size_type >;

    facet_unordered_set visited_;
    facet_unordered_set bth_facets_;
    facet_unordered_set not_bth_facets_;

    bool
    is_invisible(size_type const _facet) const
    {
        return (0 == not_bth_facets_.count(_facet)) && (0 == bth_facets_.count(_facet));
    }

    bool
    process_visibles(size_type const _facet, point const & _apex)
    {
        visited_.insert(_facet);
        facet const & facet_ = facets_[_facet];
        if (eps < facet_.distance(_apex)) {
            bool bth_ = false;
            for (size_type const neighbour : facet_.neighbours_) {
                if (visited_.count(neighbour) == 0) {
                    if (process_visibles(neighbour, _apex)) {
                        bth_ = true;
                    }
                } else if (!bth_) {
                    if (is_invisible(neighbour)){
                        bth_ = true;
                    }
                }
            }
            if (bth_) {
                bth_facets_.insert(_facet);
            } else {
                not_bth_facets_.insert(_facet);
            }
            return false;
        } else {
            return true;
        }
    }

    void
    clear_bth()
    {
        visited_.clear();
        bth_facets_.clear();
        not_bth_facets_.clear();
    }

    void
    replace_neighbour(size_type const _facet, size_type const _from, size_type const _to)
    {
        if (_from == _to) {
            return;
        }
        facet_array & neighbours_ = facets_[_facet].neighbours_;
        for (size_type & neighbour_ : neighbours_) {
            if (neighbour_ == _from) {
                neighbour_ = _to;
                return;
            }
        }
    }

public : // largest possible simplex heuristic, convex hull algorithm

    point_list
    create_simplex()
    {
        assert(1 < dimension_);
        size_type const input_size = points_.size();
        point_list basis_;
        if (!(0 < input_size)) {
            return basis_;
        }
        point_list internal_set_;
        for (size_type i = 1; i < input_size; ++i) {
            internal_set_.push_back(i);
        }
        basis_.push_back(0);
        steal_best(internal_set_, basis_);
        if (basis_.size() != 2) {
            return basis_; // can't find affine independent second point
        }
        internal_set_.splice(std::cend(internal_set_), basis_, std::cbegin(basis_)); // rejudge 0-indexed point
        for (size_type i = 1; i < dimension_; ++i) {
            steal_best(internal_set_, basis_);
            if (basis_.size() != i + 1) {
                return basis_; // can't find (i + 1) affine independent point
            }
        }
        value_type hypervolume_ = steal_best(internal_set_, basis_);
        if (basis_.size() != dimension_ + 1) {
            return basis_; // can't find affine independent (d + 1) point
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
            point_array & ordered_vertices_ = ordered_.back();
            std::sort(std::begin(ordered_vertices_), std::end(ordered_vertices_));
            rank(partition(facet_, internal_set_), newfacet);
        }
        { // adjacency
            for (size_type i = 0; i <= dimension_; ++i) {
                facet_array & neighbours_ = facets_[i].neighbours_;
                for (size_type j = 0; j <= dimension_; ++j) {
                    if (i != j) {
                        neighbours_.push_back(j);
                    }
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
        point_array vertices_;
        facet_array neighbours_;
        point_array ridge_; // horizon ridge + furthest point = new facet
        facet_array newfacets_;
        for (size_type best_facet = get_best_facet(); best_facet != facets_.size(); best_facet = get_best_facet()) {
            point_list & best_facet_outsides_ = facets_[best_facet].outside_;
            assert(!best_facet_outsides_.empty());
            size_type const apex = best_facet_outsides_.front();
            best_facet_outsides_.pop_front();
            process_visibles(best_facet, points_[apex]);
            assert(outside_.empty());
            for (size_type const not_bth_facet : not_bth_facets_) {
                facet & facet_ = facets_[not_bth_facet];
                outside_.splice(outside_.cend(), std::move(facet_.outside_));
                facet_.vertices_.clear();
                facet_.neighbours_.clear();
                unrank_and_remove(not_bth_facet);
            }
            assert(newfacets_.empty());
            for (size_type const bth_facet : bth_facets_) {
                facet & facet_ = facets_[bth_facet];
                outside_.splice(outside_.cend(), std::move(facet_.outside_));
                vertices_ = std::move(facet_.vertices_);
                neighbours_ = std::move(facet_.neighbours_);
                unrank_and_remove(bth_facet);
                for (size_type const neighbour : neighbours_) {
                    if (is_invisible(neighbour)) {
                        point_array const & horizon_ = ordered_[neighbour];
                        {
                            assert(ridge_.empty());
                            ridge_.reserve(dimension_);
                            for (size_type const vertex : vertices_) { // facets intersection with keeping of points order as it is in visible facet
                                if (std::binary_search(std::cbegin(horizon_), std::cend(horizon_), vertex)) {
                                    ridge_.push_back(vertex);
                                } else {
                                    ridge_.push_back(apex);
                                }
                            }
                            assert(ridge_.size() == dimension_); // facet
                        }
                        size_type const newfacet = add_facet(std::move(ridge_), neighbour);
                        newfacets_.push_back(newfacet);
                        replace_neighbour(neighbour, bth_facet, newfacet);
                    }
                }
            }
            clear_bth();
            adjacency(newfacets_);
            for (size_type const newfacet : newfacets_) {
                rank(partition(facets_[newfacet], outside_), newfacet);
            }
            newfacets_.clear();
            outside_.clear();
        }
        assert(outside_.empty());
        assert(ranking_.empty());
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
        ordered_.clear();
        ordered_.shrink_to_fit();
    }

};
