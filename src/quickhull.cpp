#include "quickhull.hpp"

#include <cstdlib>

int main()
{
    using H = convex_hull< double >;
    using point = typename H::point_type;
    using points = typename H::points_type;
    H convex_hull_;
    points points_{{0.0, 2.0},
                   {0.0, -1.0}};
    point point_{-10.0, 2.0};
    std::cout << convex_hull_.signed_distance_to_hyperplane(points_, point_) << std::endl;
    return EXIT_SUCCESS;
}
