#include "quickhull.hpp"

#include <cstdlib>

int main()
{
    boost::numeric::ublas::matrix< double > m_(3, 3);
    m_(0, 0) = 1;
    m_(0, 1) = 2;
    m_(0, 2) = 3;
    m_(1, 0) = 4;
    m_(1, 1) = 5;
    m_(1, 2) = 6;
    m_(2, 0) = 7;
    m_(2, 1) = 8;
    m_(2, 2) = 0;
    std::cout << determinant(m_) << std::endl;
    return EXIT_SUCCESS;
}
