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
    }

    size_type dimension_;
    point_refs_type points_;

    struct hyperplane // oriented hyperplane
    {

        point_type unit_normal_;
        point_type offset_; // offset from the origin

    };

    struct facet;

    using facet_set = std::set< size_type >;

    struct facet // (d - 1)-dimensional faces
    {

        template< typename ForwardIterator >
        facet(ForwardIterator first, ForwardIterator mid, ForwardIterator last,
              boolean_type const _outward)
            : vertices_(first, std::prev(mid))
            , points_(vertices_.cbegin(), vertices_.cend())
            , outward_(_outward)
        {
            auto const rest = vertices_.insert(vertices_.cend(), mid, last);
            points_.insert(points_.cend(), rest, vertices_.end());
            std::sort(points_.begin(), points_.end());
        }

        facet(points_type && _vertices,
              boolean_type const _outward,
              size_type const _neighbour)
            : vertices_(std::move(_vertices))
            , points_(vertices_.cbegin(), vertices_.cend())
            , outward_(_outward)
            , neighbours_({_neighbour})
        {
            std::sort(points_.begin(), points_.end());
        }

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
        below(G const & _volume) const
        {
            if (outward_) {
                return (G(0.0L) < _volume);
            } else {
                return (_volume < -G(0.0L));
            }
        }

        points_type vertices_; // oriented
        points_type points_;   // non-oriented
        boolean_type outward_;
        facet_set neighbours_;
        points_type outside_set_; // if not empty, then first point is best for this facet

        friend
        std::ostream &
        operator << (std::ostream & _out, facet const & _facet)
        {
            _out << "d " << std::boolalpha << _facet.outward_ << " : v ";
            std::copy(_facet.vertices_.cbegin(), _facet.vertices_.cend(), std::ostream_iterator< size_type >(_out, ";"));
            _out << " : n ";
            std::copy(_facet.neighbours_.cbegin(), _facet.neighbours_.cend(), std::ostream_iterator< size_type >(_out, ";"));
            _out << " : o ";
            std::copy(_facet.outside_set_.cbegin(), _facet.outside_set_.cend(), std::ostream_iterator< size_type >(_out, ";"));
            return _out;
        }

    };

    using facets_map = std::map< size_type, facet >;
    using facet_iterator = typename facets_map::iterator;
    using facets_type = std::deque< size_type >;

    facets_map facets_;

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
        auto furthest = it;
        while (++it != end) {
            G const v_ = volume(_to, *it);
            if (abs(volume_) < abs(v_)) {
                volume_ = v_;
                furthest = it;
            }
        }
        if (!(G(0.0L) < abs(volume_))) {
            throw bad_geometry("can't find linearly independent point");
        }
        _to.splice(_to.end(), _from, furthest);
        return volume_;
    }

    using ranking_type = std::multimap< G, size_type >;
    using ranking_meta_type = std::map< size_type, typename ranking_type::iterator >;
    ranking_type ranking_;
    ranking_meta_type ranking_meta_;

    void
    rank(G const _volume, size_type const _facet)
    {
        if (G(0.0L) < _volume) {
            auto const r = ranking_.emplace(_volume, _facet);
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
        G volume_(0.0L);
        while (it != end) {
            auto const next = std::next(it);
            G const v_ = volume(_facet.vertices_, *it);
            if (_facet.below(v_)) {
                if (outside_set_.empty() || _facet.further(volume_, v_)) {
                    volume_ = v_;
                    outside_set_.push_front(*it);
                } else {
                    outside_set_.push_back(*it);
                }
                _points.erase(it);
            }
            it = next;
        }
        return abs(volume_);
    }

    struct counter
            : std::iterator< std::output_iterator_tag, void, void, void, void >
    {

        counter(size_type & _counter)
            : counter_(_counter)
        { ; }

        counter &
        operator ++ ()
        {
            ++counter_;
            return *this;
        }

        counter
        operator ++ (int)
        {
            return {counter_++};
        }

        counter &
        operator * ()
        {
            return *this;
        }

        template< typename T >
        counter &
        operator = (T &&)
        {
            return *this;
        }

        operator size_type () const
        {
            return counter_;
        }

    private :

        size_type & counter_;

    };

    void
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
        for (auto exclusive = vend; exclusive != vbeg; --exclusive) {
            size_type const n = facets_.size();
            auto const f = facets_.emplace_hint(fend, n, facet(vbeg, exclusive, vend, outward_));
            facet & facet_ = f->second;
            G const v_ = partition(facet_, point_list_);
            rank(v_, n);
            assert(outward_ == !(G(0.0L) < volume(facet_.vertices_, inner_point_)));
            outward_ = !outward_;
        }
        assert(dimension_ + 1 == facets_.size());
        std::cout << "inner points: ";
        std::copy(point_list_.cbegin(), point_list_.cend(), std::ostream_iterator< size_type >(std::cout, ";"));
        std::cout << "inner points completely removed" << std::endl;
        point_list_.clear();
        std::cout << std::endl;
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
        std::cout << "initial simplex:" << std::endl;
        for (auto const & f : facets_) {
            std::cout << "initial facet #" << f.first << " = " << f.second << std::endl;
        }
        std::cout << std::endl;
    }

    void
    create_convex_hull()
    {
        create_simplex();
        size_type facet_key = facets_.size(); // unique key for facets_
        assert(facet_key == dimension_ + 1);
        auto const fend = facets_.end();
        for (size_type furthest = get_furthest(facet_key); furthest != facet_key; furthest = get_furthest(facet_key)) {
            std::cout << furthest << " ? " << facet_key << std::endl;
            facet & facet_ = facets_.at(furthest);
            std::cout << "best facet #" << furthest << " = " << facet_ << std::endl;
            facet_set visible_facets_{furthest};
            size_type const apex = facet_.outside_set_.front();
            facet_.outside_set_.pop_front();
            std::cout << "furthest point is p#" << apex << std::endl;
            {
                facet_set pool_ = facet_.neighbours_;
                facet_set visited_{furthest};
                while (!pool_.empty()) {
                    auto const first = pool_.begin();
                    size_type const f = *first;
                    facet const & candidate_ = facets_.at(f);
                    if (candidate_.below(volume(candidate_.vertices_, apex))) { // if point is above the neighbour, then add they to visible set
                        visible_facets_.insert(f);
                        std::cout << "is visible facet #" << f << std::endl;
                        std::set_difference(candidate_.neighbours_.cbegin(), candidate_.neighbours_.cend(),
                                            visited_.cbegin(), visited_.cend(),
                                            std::inserter(pool_, pool_.end()));
                    } else {
                        std::cout << "is invisible facet #" << f << std::endl;
                    }
                    visited_.insert(f);
                    pool_.erase(first);
                    std::cout << "pool = ";
                    std::copy(pool_.cbegin(), pool_.cend(), std::ostream_iterator< size_type >(std::cout, ";"));
                    std::cout << std::endl;
                }
            }
            std::cout << "visible facets: ";
            std::copy(visible_facets_.cbegin(), visible_facets_.cend(), std::ostream_iterator< size_type >(std::cout, ";"));
            std::cout << std::endl;
            std::cout << "for each visible facet (boundary finding):" << std::endl;
            // the boundary of visible_facets_ is the set of horizon ridges
            // Each ridge signifies the adjacency of two facets.
            facets_type newfacets_;
            auto const vfend = visible_facets_.end();
            for (size_type const v : visible_facets_) {
                facet const & visible_facet_ = facets_.at(v);
                std::cout << " visible facet #" << v << " = " << visible_facet_ << std::endl;
                std::cout << " for each visible facet neighbours:" << std::endl;
                points_type const & vertices_ = visible_facet_.vertices_;
                for (size_type const n : visible_facet_.neighbours_) {
                    if (visible_facets_.find(n) == vfend) { // facets intersection with keeping of points order
                        facet & horizon_facet_ = facets_.at(n);
                        std::cout << "  beyond the horizon neighbouring facet #" << n << " = " << horizon_facet_ << std::endl;
                        point_set horizon_(horizon_facet_.vertices_.cbegin(),
                                           horizon_facet_.vertices_.cend()); // n * log(n) +
                        auto const hend = horizon_.end();
                        points_type ridge_;
                        for (size_type const p : vertices_) { // n *
                            auto const h = horizon_.find(p); // (log(n) +
                            if (h == hend) {
                                ridge_.push_back(apex);
                            } else {
                                ridge_.push_back(p);
                                horizon_.erase(h); // const)
                            }
                        }std::cout << std::flush;
                        assert(horizon_.size() == 1); // horizon_ contains the only invisible point beyond the horizon
                        assert(ridge_.size() == dimension_); // ridge_ contains newfacet vertices (ridge + current furthest point)
                        { // replace visible facet became internal with newly created facet
                            std::cout << "  for this neighbouring facet repalce old (visible) neighbouring facet #" << v << " by #" << facet_key << std::endl;
                            horizon_facet_.neighbours_.erase(v);
                            horizon_facet_.neighbours_.insert(horizon_facet_.neighbours_.cend(), facet_key);
                        }
                        newfacets_.push_back(facet_key);
                        auto const f = facets_.emplace_hint(fend, facet_key, facet(std::move(ridge_), visible_facet_.outward_, n));
                        ++facet_key;
                        std::cout << "  new facet #" << f->first << " = " << f->second << std::endl;
                    }
                }
            }
            std::cout << std::endl;
            std::cout << " full list of newly created facets: ";
            std::copy(newfacets_.cbegin(), newfacets_.cend(), std::ostream_iterator< size_type >(std::cout, ";"));
            std::cout << std::endl;
            std::cout << " for newly created facets add neighbouring facets:" << std::endl;
            {
                auto const nend = newfacets_.end();
                for (auto first = newfacets_.begin(); first != nend; ++first) {
                    size_type const f = *first;
                    facet & first_ = facets_.at(f);
                    points_type const & ofirst_ = first_.points_;
                    for (auto second = std::next(first); second != nend; ++second) {
                        size_type const s = *second;
                        facet & second_ = facets_.at(s);
                        points_type const & osecond_ = second_.points_;
                        size_type count_ = 0;
                        std::set_difference(ofirst_.cbegin(), ofirst_.cend(),
                                            osecond_.cbegin(), osecond_.cend(),
                                            counter{count_});
                        if (count_ == 1) {
                            first_.neighbours_.insert(s);
                            second_.neighbours_.insert(f);
                            std::cout << "new neighbours: #" << f << " + #" << s << std::endl;
                        }
                    }
                }
            }
            std::cout << " for each visible facet steal its outside set into single pool" << std::endl;
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
                std::cout << "  remove visible facet #" << v << std::endl;
            }
            std::cout << " full list of possible outside points for newly created facets: ";
            std::copy(outside_set_.cbegin(), outside_set_.cend(), std::ostream_iterator< size_type >(std::cout, ";"));
            std::cout << std::endl;
            for (size_type const n : newfacets_) {
                facet & newfacet_ = facets_.at(n);
                G const v_ = partition(newfacet_, outside_set_);
                rank(v_, n);
                //assert(outward_ == !(G(0.0L) < volume(facet_.vertices_, inner_point_)));
            }
            std::cout << furthest << " ? " << facet_key << std::endl;
            std::cout << "facets count: " << facets_.size() << std::endl << std::endl;
        }
        std::cout << "result:" << std::endl;
        for (auto const & facet_ : facets_) {
            std::cout << "result #" << facet_.first << " = " << facet_.second << std::endl;
        }
    }

};
