*NOTE: This library is header-only.*

Implementation of the **Quickhull algorithm** (Barber et al) for the convex hulls finding in arbitrary dimension (>1) space. Also implemented the **Mehlhorn algorithm** (Mehlhorn et al) for checking convexity of resulting geometric structure.

- [Barber, C. B., D.P. Dobkin, and H.T. Huhdanpaa, 1995. "The Quickhull Algorithm for Convex Hulls", ACM Transactions on Mathematical Software.](https://www.cs.princeton.edu/~dpd/Papers/BarberDobkinHuhdanpaa.pdf)
- [Kurt Mehlhorn, Stefan Näher, Thomas Schilz, Stefan Schirra, Michael Seel, Raimund Seidel, and Christian Uhrig. "Checking geometric programs or verification of geometric structures", In Proc. 12th Annu. ACM Sympos. Comput. Geom., pages 159–165, 1996](https://people.mpi-inf.mpg.de/~mehlhorn/ftp/programc.ps)

Example:
```cpp
#include <quickhull.hpp>

#include <array>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

#include <cstdlib>

int main()
{
    using F = float;
    constexpr std::size_t dim = 3;
    using Points = std::vector<std::array<F, dim>>;

    Points points(10); // input

    { // fill it somehow (use real data)
        std::mt19937 gen;
        for (auto & [x, y, z] : points) {
            x = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
            y = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
            z = std::generate_canonical<F, std::numeric_limits<F>::digits>(gen);
        }
    }

    const auto eps = std::numeric_limits<F>::epsilon();
    quick_hull<typename Points::const_iterator> qh{dim, eps};
    qh.add_points(std::cbegin(points), std::cend(points));
    auto initial_simplex = qh.get_affine_basis();
    if (initial_simplex.size() < dim + 1) {
        return EXIT_FAILURE; // degenerated input set
    }
    qh.create_initial_simplex(std::cbegin(initial_simplex), std::prev(std::cend(initial_simplex)));
    qh.create_convex_hull();
    if (!qh.check()) {
        return EXIT_FAILURE; // resulted structure is not convex (generally due to precision errors)
    }

    qh.facets_; // use as result
}
```
