#include <boost/program_options.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <valarray>
#include <deque>
#include <map>
#include <random>
#include <limits>
#include <chrono>

#include <cmath>
#include <cstdlib>

template< typename G >
struct randombox
{

    using size_type = std::size_t;

    G const eps = std::numeric_limits< G >::epsilon();
    G const zero = G(0);
    G const one = G(1);

    using seed_type = typename std::random_device::result_type;
    seed_type seed_;
    std::mt19937_64 random_;

    void
    set_seed(seed_type const _seed)
    {
        seed_ = _seed;
        random_.seed(seed_);
    }

    void
    set_seed()
    {
#if 0
        std::random_device rd_;
        seed_ = rd_();
#else
        seed_ = std::chrono::high_resolution_clock::now().time_since_epoch().count();
#endif
        random_.seed(seed_);
    }

    using point_type = std::valarray< G >;
    using points_type = std::deque< point_type >;
    using mask_array_type = std::valarray< bool >;

    size_type dimension_ = 3;
    points_type points_;
    size_type count_;

    std::istream &
    operator () (std::istream & _in)
    {
        if (!!_in) {
            std::string line_;
            if (!std::getline(_in, line_)) {
                throw std::runtime_error("input: no 'dimensionality' value at first line");
            }
            std::istringstream iss_(line_);
            if (!(iss_ >> dimension_)) {
                throw std::runtime_error("input: bad 'dimensionality' value at first line");
            }
            if (!std::getline(_in, line_)) {
                throw std::runtime_error("input: no 'count' value at second line");
            }
            iss_.str(line_);
            size_type size_;
            if (!(iss_ >> size_)) {
                throw std::runtime_error("input: bad 'count' value at second line");
            }
            points_.resize(size_, point_type(dimension_));
            size_type i = 0;
            while (i < size_) {
                if (!std::getline(_in, line_) || line_.empty()) {
                    throw std::runtime_error("input: empty line or no 'count' lines with points coordinates");
                }
                if (line_.front() != '#') {
                    iss_.str(line_);
                    for (G & component_ : points_[i]) {
                        if (!(iss_ >> component_)) {
                            throw std::runtime_error("input: bad coordinate value");
                        }
                    }
                    ++i;
                }
            }
        }
        return _in;
    }

    std::ostream &
    operator () (std::ostream & _out) const
    {
        assert(0 < dimension_);
        std::ios state_(nullptr);
        state_.copyfmt(_out);
        {
            _out << dimension_ << '\n';
            size_type const size_ = points_.size();
            assert(0 < size_);
            _out << size_ << '\n';
            _out.precision(std::numeric_limits< G >::digits10);
            for (point_type const & point_ : points_) {
                assert(point_.size() == dimension_);
                auto const last = std::prev(std::end(point_));
                for (auto it = std::begin(point_); it != last; ++it) {
                    _out << *it << ' ';
                }
                _out << *last << '\n';
            }
        }
        _out.copyfmt(state_);
        return _out;
    }

    void
    set_dimension(size_type const _dimension)
    {
        if (0 < _dimension) {
            if (dimension_ != _dimension) {
                for (point_type & point_ : points_) {
                    point_type storage_ = std::move(point_);
                    point_.resize(_dimension, G(0));
                    std::copy_n(std::begin(storage_), std::min(dimension_, _dimension), std::begin(point_));
                }
                dimension_ = _dimension;
            }
        }
    }

    void
    set_count(size_type const _count)
    {
        if ((points_.size() == 0) && (_count == 0)) {
            count_ = dimension_ + 1;
        } else {
            count_ = _count;
        }
    }

    void
    add_sphere()
    {
        std::normal_distribution< G > N_; // N(0, 1) distribution
        point_type source_(dimension_);
        //points_.reserve(count_);
        while (points_.size() < count_) {
            for (size_type j = 0; j < dimension_; ++j) {
                source_[j] = N_(random_);
            }
            points_.push_back(source_);
            source_ *= source_;
            using std::sqrt;
            G norm_ = sqrt(source_.sum());
            if (norm_ < eps) {
                points_.pop_back();
            } else {
                points_.back() *= (one / std::move(norm_));
            }
        }
    }

    void
    add_ball()
    {
        add_sphere();
        std::uniform_real_distribution< G > UR_(zero,
                                                std::nextafter(one, std::numeric_limits< G >::infinity())); // uniform [0;1] ditribution
        for (size_type i = 0; i < count_; ++i) {
            point_type & destination_ = points_[i];
            using std::pow;
            destination_ *= pow(UR_(random_), one / G(dimension_));
        }
    }

    void
    add_unit_simplex_face()
    {
        std::uniform_real_distribution< G > UR_(std::nextafter(zero, one),
                                                std::nextafter(one, std::numeric_limits< G >::infinity())); // uniform (0;1] ditribution
        for (size_type i = 0; i < count_; ++i) {
            points_.emplace_back(dimension_);
            point_type & destination_ = points_.back();
            for (size_type j = 0; j < dimension_; ++j) {
                destination_[j] = UR_(random_);
            }
            destination_ = std::log(destination_);
            G norm_ = destination_.sum();
            if (norm_ == -std::numeric_limits< G >::infinity()) { // if some of logarithms of generated values is -HUGE_VAL, then the correspoinding nonnormalized value is one
                mask_array_type const ones_ = (destination_ == -std::numeric_limits< G >::infinity()); // store into std::valarray< bool > to prevent evaluations to being lazy
                destination_[ones_] = one;
                destination_[!ones_] = zero;
                norm_ = destination_.sum(); // number of close-to-zero generated values, can be zero (if there just an overflow)
            }
            if (-eps < norm_) { // if generated random point is too close to the origin, then assume, that origin is good choise
                destination_ = zero;
            } else {
                destination_ *= (one / std::move(norm_));
            }
        }
    }

};

template< typename G >
std::istream &
operator >> (std::istream & _in, randombox< G > & _randombox)
{
    return _randombox(_in);
}

template< typename G >
std::ostream &
operator << (std::ostream & _out, randombox< G > const & _randombox)
{
    return _randombox(_out);
}

int
main(int ac, char * av[])
{
    using size_type = std::size_t;
    using G = long double;

    //std::istream & in_ = std::cin;
    std::ostream & out_ = std::cout;

    using randombox_type = randombox< G >;
    using seed_type = typename randombox_type::seed_type;
    randombox_type randombox_;
    //in_ >> randombox_;

    enum class geometrical_object
    {
        sphere,
        ball,
        cube,
        diamond,
        unit_simplex_face,
    };

    std::map< std::string, geometrical_object > const gobject_map_{
        {"sphere",       geometrical_object::sphere},
        {"ball",         geometrical_object::ball},
        {"cube",         geometrical_object::cube},
        {"diamond",      geometrical_object::diamond},
        {"unit-simplex-face", geometrical_object::unit_simplex_face},
    };

    namespace po = boost::program_options;

    po::options_description options_("Information options");
    options_.add_options()
            ("help", "produce this help message")
            ("dimension,D", po::value< size_type >()->implicit_value(0), "dimensionality value")
            ("count,N", po::value< size_type >()->implicit_value(0), "count of points generated (can be specified without the key)")
            //("bounding-box,B", po::value< G >(), "bounding box coordinates")
            //("mesh,M", po::value< std::string >(), "lattice (mesh) rotated by {{n, -m, 0}, {m, n, 0}, {0, 0, r}}. Skipping the r, makes r = sqrt(n^2 + m^2). m = 1, n = 0 is orthogonal lattice")
            //("cospherical,S", "cospherical points randomly generated in a cube and projected to the unit sphere")
            //("in-simplex,X", "simplicial distribution, may be used with --width")
            //("simplex,Y", "simplicial distribution with simplex inclusive")
            //("width,W", "restrict points to distance of the surface")
            //("cube,C", "add a unit cube to the output")
            //("diamond,D", "add a unit diamond to the output")
            //("point,P", "add point with specified coordinates to the output")
            //("offset,O", "offset the data by adding specified component to each of the coordinates")
            //("rotate,R", po::value< std::string >(), "rotate the data around the specified axis by the specified angle")
            ("seed", po::value< seed_type >(), "use specified value as random number seed")
            //("integer", "generate integer coordinates in specified integer bounding box")
            ("add", po::value< std::string >(), "add specific object to the output. Possible object types: sphere, ball, cube, diamond, unit-simplex-face")
            ;

    po::positional_options_description positional_;
    positional_.add("count", 1);

    po::variables_map vm_;
    po::store(po::command_line_parser(ac, av).options(options_).positional(positional_).run(), vm_);
    po::notify(vm_);

    auto const vmend = vm_.cend();
    auto it = vm_.find("help");
    if (it != vmend) {
        std::clog << options_ << std::endl;
        return EXIT_SUCCESS;
    }
    it = vm_.find("seed");
    if (it != vmend) {
        randombox_.set_seed(it->second.template as< seed_type >());
    } else {
        randombox_.set_seed();
    }
    it = vm_.find("dimension");
    if (it != vmend) {
        randombox_.set_dimension(it->second.template as< size_type >());
    }
    it = vm_.find("count");
    if (it != vmend) {
        randombox_.set_count(it->second.template as< size_type >());
    }
    it = vm_.find("add");
    if (it != vmend) {
        auto const go = gobject_map_.find(it->second.template as< std::string >());
        if (go != gobject_map_.cend()) {
            switch (go->second) {
            case geometrical_object::sphere : {
                randombox_.add_sphere();
                break;
            }
            case geometrical_object::ball : {
                randombox_.add_ball();
                break;
            }
            case geometrical_object::cube : {
                break;
            }
            case geometrical_object::diamond : {
                break;
            }
            case geometrical_object::unit_simplex_face : {
                randombox_.add_unit_simplex_face();
                break;
            }
            default : {
                throw std::runtime_error("unsupported geometrical object");
            }
            }
        } else {
            throw std::runtime_error("bad geometrical object name");
        }
    }

    out_ << randombox_ << std::endl;
    return EXIT_SUCCESS;
}
