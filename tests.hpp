//v10 13.06.2015

#ifndef TESTS_HEADER
#define TESTS_HEADER

#include "constants.hpp"

void test_maketimes(void);
void test_read_single(void);
void test_read_consecutive(void);
void test_pi_size(void);
void test_mixing_parameter(void);
void test_pi(void);
void test_n_topologcial(void);
void test_n_metric(void);
void test_n_voronoi(void);
void test_cint(void);
void test_average(void);
void test_distance_periodic(void);
void test_scalar_product(void);
void test_matrix_add(void);
void test_get_average_number_of_neighbours(void);

bool is_number_of_neighbours_equal(const MATRIX &n);
bool is_symmetric(const MATRIX &n);

template <typename T> bool is_weakly_diagonally_dominant_and_non_negativ_diagonal_entries(const T &m)
{
    for (unsigned int r = 0; r < m.rows(); ++r)
    {
        real absrowsumwithoutdiag = 0;
        for (unsigned int c = 0; c < m.cols(); ++c)
        {
            if (r != c)
            {
                absrowsumwithoutdiag += std::abs(m(r, c));
            }
        }

        if (absrowsumwithoutdiag > std::abs(m(r,r)) || m(r,r) < 0)
        {
            return false;
        }
    }
    return true;
}

#endif //TESTS_HEADER
