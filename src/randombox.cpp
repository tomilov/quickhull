#include <boost/program_options.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <valarray>
#include <deque>
#include <vector>
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

    using seed_type = typename std::mt19937_64::result_type;
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

    size_type const default_dimension = 3;
    size_type dimension_ = default_dimension;
    points_type source_points_;
    points_type separate_points_;
    points_type resulting_points_;
    size_type count_ = 0;
    G bounding_box_ = one;

    std::normal_distribution< G > N_; // standard normal distribution
    std::uniform_real_distribution< G > zero_to_one_; // uniform [0;1] ditribution

    randombox()
        : zero_to_one_(zero, std::nextafter(one, one + one)) // ? std::nextafter(zero, one)
    { ; }

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
            iss_.clear();
            iss_.str(line_);
            size_type size_;
            if (!(iss_ >> size_)) {
                throw std::runtime_error("input: bad 'count' value at second line");
            }
            iss_.clear();
            for (size_type i = 0; i < size_; ++i) {
                if (!std::getline(_in, line_) || line_.empty()) {
                    throw std::runtime_error("input: empty line or no 'count' lines with points coordinates");
                }
                if (line_.front() != '#') {
                    iss_.str(line_);
                    source_points_.emplace_back(dimension_);
                    for (G & component_ : source_points_.back()) {
                        if (!(iss_ >> component_)) {
                            throw std::runtime_error("input: bad point format");
                        }
                    }
                    iss_.clear();
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
            size_type const size_ = resulting_points_.size();
            //assert(0 < size_);
            _out << size_ << '\n';
            _out.precision(std::numeric_limits< G >::digits10);
            for (point_type const & point_ : resulting_points_) {
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
                size_type const subdimension_ = std::min(dimension_, _dimension);
                for (point_type & point_ : source_points_) {
                    point_type storage_ = std::move(point_);
                    point_.resize(_dimension, zero);
                    std::copy_n(std::begin(storage_), subdimension_, std::begin(point_));
                }
                dimension_ = _dimension;
            }
        }
    }

    point_type
    get_point(std::string const & _components)
    {
        point_type point_(zero, dimension_);
        if (!_components.empty()) { // implicit point is origin
            std::istringstream iss_(_components);
            for (G & component_ : point_) {
                if (!(iss_ >> component_)) {
                    throw std::runtime_error("input: bad coordinate value");
                }
            }
        }
        return point_;
    }

    void
    add_point(point_type && _point)
    {
        separate_points_.push_back(std::move(_point));
    }

    void
    add_point(point_type const & _point)
    {
        separate_points_.push_back(_point);
    }

    void
    add_point(std::string const & _components)
    {
        return add_point(get_point(_components));
    }

    void
    set_count(size_type const _count)
    {
        if (_count == 0) {
            count_ = dimension_ + 1;
        } else {
            count_ = _count;
        }
    }

    void
    set_bounding_box(G const & _bounding_box)
    {
        assert(eps < _bounding_box);
        bounding_box_ = _bounding_box;
    }

    void
    pick_unit_cube_point(point_type & _point)
    {
        for (G & component_ : _point) {
            component_ = zero_to_one_(random_);
        }
    }

    point_type
    pick_unit_cube_point(size_type const _dimension)
    {
        point_type point_(_dimension);
        pick_unit_cube_point(point_);
        return point_;
    }

    void
    add_unit_cube()
    {
        for (size_type i = 0; i < count_; ++i) {
            resulting_points_.push_back(pick_unit_cube_point(dimension_));
        }
    }

    void
    generate_parallelotope()
    {
        assert(!(dimension_ + 1 < source_points_.size()));
        assert(1 < source_points_.size());
        point_type const & vertex_ = separate_points_.front();
        auto const vbeg = std::next(separate_points_.cbegin());
        auto const vend = separate_points_.cend();
        point_type point_(separate_points_.size() - 1);
        for (size_type i = 0; i < count_; ++i) {
            pick_unit_cube_point(point_);
            resulting_points_.push_back(std::inner_product(vbeg, vend, std::begin(point_), vertex_));
        }
    }

    void
    add_diamond_surface()
    {
        add_unit_simplex();
        std::uniform_int_distribution< size_type > flip_sign_(0, 1);
        for (point_type & point_ : resulting_points_) {
            for (G & component_ : point_) {
                if (flip_sign_(random_) == 0) {
                    component_ = -component_;
                }
            }
        }
    }

    void
    add_diamond_solid()
    {
        std::uniform_int_distribution< size_type > flip_sign_(0, 1);
        point_type point_(dimension_ + 1);
        for (size_type i = 0; i < count_; ++i) {
            pick_uint_simplex_point(point_);
            resulting_points_.emplace_back(dimension_);
            point_type & destination_ = resulting_points_.back();
            for (size_type j = 0; j < dimension_; ++j) {
                if (flip_sign_(random_) == 0) {
                    destination_[j] =  point_[j];
                } else {
                    destination_[j] = -point_[j];
                }
            }
        }
    }

    void
    pick_uint_simplex_point(point_type & _point)
    {
        pick_unit_cube_point(_point);
        _point = -std::log(_point);
        G norm_ = _point.sum();
        if (norm_ == std::numeric_limits< G >::infinity()) { // if some of logarithms of generated values is -HUGE_VAL, then the correspoinding non-normalized value is one
            mask_array_type const ones_ = (_point == std::numeric_limits< G >::infinity()); // store into std::valarray< bool > to prevent evaluations to being lazy
            _point[ones_] = one;
            _point[!ones_] = zero;
            norm_ = _point.sum(); // number of close-to-zero generated values, can be zero (if there just an overflow)
        }
        if (eps < norm_) {
            _point *= (one / std::move(norm_));
        } else {
            _point = zero; // if generated random point is too close to the origin, then assume, that origin is good choise
        }
    }

    point_type
    pick_uint_simplex_point(size_type const _dimension)
    {
        point_type point_(_dimension);
        pick_uint_simplex_point(point_);
        return point_;
    }

    void
    add_unit_simplex()
    {
        for (size_type i = 0; i < count_; ++i) {
            resulting_points_.push_back(pick_uint_simplex_point(dimension_));
        }
    }

    void
    generate_simplex()
    {
        assert(!(dimension_ + 1 < source_points_.size()));
        assert(1 < source_points_.size());
        auto const pbeg = separate_points_.cbegin();
        auto const pend = separate_points_.cend();
        point_type point_(separate_points_.size());
        for (size_type i = 0; i < count_; ++i) {
            pick_uint_simplex_point(point_);
            resulting_points_.push_back(std::inner_product(pbeg, pend, std::begin(point_), point_type(zero, dimension_)));
        }
    }

    void
    add_sphere()
    {
        point_type source_(dimension_);
        while (resulting_points_.size() < count_) {
            for (size_type j = 0; j < dimension_; ++j) {
                source_[j] = N_(random_);
            }
            resulting_points_.push_back(source_);
            source_ *= source_;
            using std::sqrt;
            G norm_ = sqrt(source_.sum());
            if (norm_ < eps) {
                resulting_points_.pop_back();
            } else {
                resulting_points_.back() *= (one / std::move(norm_));
            }
        }
    }

    void
    add_ball()
    {
        add_sphere();
        G const power_ = (one / G(dimension_));
        for (point_type & destination_ : resulting_points_) {
            using std::pow;
            destination_ *= pow(zero_to_one_(random_), power_);
        }
    }

    void
    project_to_cylinder()
    {
        assert(separate_points_.size() == 1);
        point_type const & element_ = separate_points_.back();
        for (size_type i = 0; i < count_; ++i) {
            resulting_points_.push_back(source_points_[i] + element_ * zero_to_one_(random_));
        }
    }

    void
    project_to_cone()
    {
        assert(separate_points_.size() == 1);
        point_type const & peak_ = separate_points_.back();
        G const power_ = (one / G(dimension_));
        for (size_type i = 0; i < count_; ++i) {
            using std::pow;
            G const p_ = pow(zero_to_one_(random_), power_);
            resulting_points_.push_back(source_points_[i] * p_ + peak_ * (one - p_));
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

    using randombox_type = randombox< G >;
    using seed_type = typename randombox_type::seed_type;
    randombox_type randombox_;

    enum class geometrical_object
    {
        sphere,
        ball,
        cube,
        diamond_surface,
        diamond_solid,
        unit_simplex,
        simplex,
        parallelotope,
        cylinder,
        cone,
    };

    std::map< std::string,  geometrical_object > const gomap_{
        {"sphere",          geometrical_object::sphere},
        {"ball",            geometrical_object::ball},
        {"cube",            geometrical_object::cube},
        {"diamond-surface", geometrical_object::diamond_surface},
        {"diamond-solid",   geometrical_object::diamond_solid},
        {"unit-simplex",    geometrical_object::unit_simplex},
        {"simplex",         geometrical_object::simplex},
        {"parallelotope",   geometrical_object::parallelotope},
        {"cylinder",        geometrical_object::cylinder},
        {"cone",            geometrical_object::cone},
    };

    namespace po = boost::program_options;

    po::options_description options_("Information options");
    options_.add_options()
            ("help", "produce this help message")
            ("input,I", po::value< std::string >()->implicit_value(""), "input file name (nothing for stdin)")
            ("output,O", po::value< std::string >()->implicit_value(""), "output file name (stdout by default)")
            ("seed", po::value< seed_type >(), "use specified value as random number seed")
            ("dimension,D", po::value< size_type >()->default_value(0), "dimensionality value")
            ("count,N", po::value< size_type >(), "count of points generated (can be specified without the key)")
            //("bounding-box,B", po::value< G >(), "bounding box coordinates")
            //("mesh,M", po::value< std::string >(), "lattice (mesh) rotated by {{n, -m, 0}, {m, n, 0}, {0, 0, r}}. Skipping the r, makes r = sqrt(n^2 + m^2). m = 1, n = 0 is orthogonal lattice")
            //("cospherical,S", "cospherical points randomly generated in a cube and projected to the unit sphere")
            //("in-simplex,X", "simplicial distribution, may be used with --width")
            //("simplex,Y", "simplicial distribution with simplex inclusive")
            //("width,W", "restrict points to distance of the surface")
            //("cube,C", "add a unit cube to the output")
            //("diamond,D", "add a unit diamond to the output")
            ("point,P", po::value< std::vector< std::string > >(), "add point with specified coordinates to the input")
            //("offset,O", "offset the data by adding specified component to each of the coordinates")
            //("rotate,R", po::value< std::string >(), "rotate the data around the specified axis by the specified angle")
            ("add", po::value< std::string >(), "Minkowski sum of input set and generated body or surface. Possible object types: sphere, ball, cube, diamond-surface, diamond-solid, unit-simplex")
            ("generate", po::value< std::string >(), "generate points inside embeddable simplex or parallelotope")
            ("project-to", po::value< std::string >(), "project from affine subspace to cone or cylinder. Possible projections types: cone, cylinder")
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
    it = vm_.find("input");
    if (it != vmend) {
        std::string const input_file_name_ = it->second.template as< std::string >();
        if (input_file_name_.empty()) {
            std::cin >> randombox_;
        } else {
            std::ifstream ifs_(input_file_name_);
            if (!ifs_) {
                std::cerr << "can't open input file" << std::endl;
                return EXIT_FAILURE;
            }
            ifs_ >> randombox_;
        }
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
    it = vm_.find("point");
    if (it != vmend) {
        for (std::string const & components_ : it->second.template as< std::vector< std::string > >()) {
            randombox_.add_point(components_);
        }
    }
    it = vm_.find("count");
    if (it != vmend) {
        randombox_.set_count(it->second.template as< size_type >());
    }
    it = vm_.find("bounding-box");
    if (it != vmend) {
        randombox_.set_bounding_box(it->second.template as< G >());
    }
    it = vm_.find("add");
    if (it != vmend) {
        auto const gom = gomap_.find(it->second.template as< std::string >());
        if (gom != gomap_.cend()) {
            switch (gom->second) {
            case geometrical_object::sphere : {
                randombox_.add_sphere();
                break;
            }
            case geometrical_object::ball : {
                randombox_.add_ball();
                break;
            }
            case geometrical_object::cube : {
                randombox_.add_unit_cube();
                break;
            }
            case geometrical_object::diamond_surface : {
                randombox_.add_diamond_surface();
                break;
            }
            case geometrical_object::diamond_solid : {
                randombox_.add_diamond_solid();
                break;
            }
            case geometrical_object::unit_simplex : {
                randombox_.add_unit_simplex();
                break;
            }
            case geometrical_object::simplex : {
                randombox_.generate_simplex();
                break;
            }
            case geometrical_object::parallelotope : {
                randombox_.generate_parallelotope();
                break;
            }
            case geometrical_object::cylinder : {
                randombox_.project_to_cylinder();
                break;
            }
            case geometrical_object::cone : {
                randombox_.project_to_cone();
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
    it = vm_.find("output");
    if (it != vmend) {
        std::string const output_file_name_ = it->second.template as< std::string >();
        if (output_file_name_.empty()) {
            std::cout << randombox_ << std::flush;
        } else {
            std::ofstream ofs_(output_file_name_);
            if (!ofs_) {
                std::cerr << "can't open output file" << std::endl;
                return EXIT_FAILURE;
            }
            ofs_ << randombox_ << std::flush;
        }
    } else {
        std::cout << randombox_ << std::flush;
    }
    return EXIT_SUCCESS;
}
