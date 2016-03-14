//v10 13.06.2015

#include "tests.hpp"

void test_pi_size(void)
{
    std::vector<real> angles = {0., 0., 0., 0., 0.};
    Angdict angdict;
    angdict[1337] = angles;
    Pi pi(angdict);
    assert(pi.size() == 5);
}

void test_mixing_parameter(void)
{
    MATRIX n(2, 2);
    n.Set(0, 0, 1); n.Set(0, 1, 0);
    n.Set(1, 0, 0); n.Set(1, 1, 1);

    MATRIX m(2, 2);
    m.Set(0, 0, 1); m.Set(0, 1, 1);
    m.Set(1, 0, 0); m.Set(1, 1, 1);

    real nc = get_average_number_of_neighbours(n);

    assert_almost_equal(mixing_parameter(0.01, n, m, nc), real(50.));
}

void test_pi(void)
{
    std::vector<real> angles;
    std::vector<Vector> pos;
    Angdict angdict;
    Posdict posdict;

    Vector a, b, c, d;
    a.x = 1.; a.y = 2.;
    b.x = 2.; b.y = 2.;
    c.x = 2.; c.y = 0.;
    d.x = 4.; d.y = 4.;

    angles = {real(0), real(PI/4), real(-PI/4), real(0)};
    pos = {a, b, c, d};

    angdict[0] = angles;
    posdict[0] = pos;

    Pi pi(angdict);
    assert(pi.Get(0, 0).Compare(Vector(0., 0.)));
    assert(pi.Get(0, 1).Compare(Vector(0., std::sqrt(real(2.))/2.)));
    assert(pi.Get(0, 2).Compare(Vector(0., -std::sqrt(real(2.))/2.)));
    assert(pi.Get(0, 3).Compare(Vector(0., 0.)));

    angles = {PI/2., 3. * PI/4., PI/4., PI/2.};

    angdict[0] = angles;
    posdict[0] = pos;

    pi = Pi(angdict);
    assert(pi.Get(0, 0).Compare(Vector(0., 0.)));
    assert(pi.Get(0, 1).Compare(Vector(-std::sqrt(real(2.))/2., 0.)));
    assert(pi.Get(0, 2).Compare(Vector(std::sqrt(real(2.))/2., 0.)));
    assert(pi.Get(0, 3).Compare(Vector(0., 0.)));


    angles = {PI/4, PI/2, 0, PI/4};

    angdict[0] = angles;
    posdict[0] = pos;

    pi = Pi(angdict);
    assert(pi.Get(0, 0).Compare(Vector(0., 0.)));
    assert(pi.Get(0, 1).Compare(Vector(-0.5, 0.5)));
    assert(pi.Get(0, 2).Compare(Vector(0.5, -0.5)));
    assert(pi.Get(0, 3).Compare(Vector(0., 0.)));
}

void test_n_topologcial(void)
{
    real L = 32.;
    std::vector<Vector> pos;

    Vector a, b, c, d;
    a.x = 1.; a.y = 2.;
    b.x = 2.; b.y = 2.;
    c.x = 2.; c.y = 0.;
    d.x = 4.; d.y = 4.;

    pos = {a, b, c, d};

    MATRIX n(4, 4);
    n.Zero();
    calculate_n_topological_unsymmetrical(n, pos, L, 2);

    MATRIX m(4, 4);
    /*
    0 1 1 0
    1 0 1 0
    1 1 0 0
    1 1 0 0
    */
    m.Set(0, 0, 0); m.Set(0, 1, 1); m.Set(0, 2, 1); m.Set(0, 3, 0);
    m.Set(1, 0, 1); m.Set(1, 1, 0); m.Set(1, 2, 1); m.Set(1, 3, 0);
    m.Set(2, 0, 1); m.Set(2, 1, 1); m.Set(2, 2, 0); m.Set(2, 3, 0);
    m.Set(3, 0, 1); m.Set(3, 1, 1); m.Set(3, 2, 0); m.Set(3, 3, 0);

    assert(n == m);
}

void test_n_metric(void)
{
    real L = 32;
    std::vector<Vector> pos;

    Vector a, b, c, d;
    a.x = 1.; a.y = 2.;
    b.x = 2.; b.y = 2.;
    c.x = 2.; c.y = 0.;
    d.x = 4.; d.y = 4.;

    pos = {a, b, c, d};

    MATRIX n(4, 4);
    n.Zero();
    calculate_n_metric(n, pos, L, 2);
    MATRIX m(4, 4);
    /*
    0 1 0 0
    1 0 1 0
    0 1 0 0
    0 0 0 0
    */
    m.Set(0, 0, 0); m.Set(0, 1, 1); m.Set(0, 2, 0); m.Set(0, 3, 0);
    m.Set(1, 0, 1); m.Set(1, 1, 0); m.Set(1, 2, 1); m.Set(1, 3, 0);
    m.Set(2, 0, 0); m.Set(2, 1, 1); m.Set(2, 2, 0); m.Set(2, 3, 0);
    m.Set(3, 0, 0); m.Set(3, 1, 0); m.Set(3, 2, 0); m.Set(3, 3, 0);

    assert(n == m);
}

void test_n_voronoi(void)
{
    std::vector<Vector> pos;
    real L = 32;

    Vector a, b, c, d;
    a.x = 1.; a.y = 2.;
    b.x = 2.; b.y = 2.;
    c.x = 2.; c.y = 0.;
    d.x = 4.; d.y = 4.;

    pos = {a, b, c, d};

    MATRIX n(4, 4);
    n.Zero();
    calculate_n_voronoi(n, pos, L);

    MATRIX m(4, 4);
    /*
    0 1 1 1
    1 0 1 1
    1 1 0 1
    1 1 1 0
    */
    m.Set(0, 0, 0); m.Set(0, 1, 1); m.Set(0, 2, 1); m.Set(0, 3, 1);
    m.Set(1, 0, 1); m.Set(1, 1, 0); m.Set(1, 2, 1); m.Set(1, 3, 1);
    m.Set(2, 0, 1); m.Set(2, 1, 1); m.Set(2, 2, 0); m.Set(2, 3, 1);
    m.Set(3, 0, 1); m.Set(3, 1, 1); m.Set(3, 2, 1); m.Set(3, 3, 0);

    assert(n == m);

    /*9 particles Voronoi*/

    MATRIX n2(9, 9);
    n2.Zero();

    Vector v1(10,30), v2(20,30), v3(30,30), v4(10,20), v5(20,20), v6(30,20), v7(10,10), v8(20,10), v9(30,10);
    pos = {v1, v2, v3, v4, v5, v6, v7, v8, v9};
    calculate_n_voronoi(n2, pos, L);

    MATRIX m2(9, 9);
    /*
    0 1 1 1 0 0 1 0 0
    1 0 1 0 1 0 0 1 0
    1 1 0 0 0 1 0 0 1
    1 0 0 0 1 1 1 0 0
    0 1 0 1 0 1 0 1 0
    0 0 1 1 1 0 0 0 1
    1 0 0 1 0 0 0 1 1
    0 1 0 0 1 0 1 0 1
    0 0 1 0 0 1 1 1 0
    */
    m2.Set(0, 0, 0);m2.Set(0, 1, 1);m2.Set(0, 2, 1);m2.Set(0, 3, 1);m2.Set(0, 4, 0);m2.Set(0, 5, 0);m2.Set(0, 6, 1);m2.Set(0, 7, 0);m2.Set(0, 8, 0);
    m2.Set(1, 0, 1);m2.Set(1, 1, 0);m2.Set(1, 2, 1);m2.Set(1, 3, 0);m2.Set(1, 4, 1);m2.Set(1, 5, 0);m2.Set(1, 6, 0);m2.Set(1, 7, 1);m2.Set(1, 8, 0);
    m2.Set(2, 0, 1);m2.Set(2, 1, 1);m2.Set(2, 2, 0);m2.Set(2, 3, 0);m2.Set(2, 4, 0);m2.Set(2, 5, 1);m2.Set(2, 6, 0);m2.Set(2, 7, 0);m2.Set(2, 8, 1);
    m2.Set(3, 0, 1);m2.Set(3, 1, 0);m2.Set(3, 2, 0);m2.Set(3, 3, 0);m2.Set(3, 4, 1);m2.Set(3, 5, 1);m2.Set(3, 6, 1);m2.Set(3, 7, 0);m2.Set(3, 8, 0);
    m2.Set(4, 0, 0);m2.Set(4, 1, 1);m2.Set(4, 2, 0);m2.Set(4, 3, 1);m2.Set(4, 4, 0);m2.Set(4, 5, 1);m2.Set(4, 6, 0);m2.Set(4, 7, 1);m2.Set(4, 8, 0);
    m2.Set(5, 0, 0);m2.Set(5, 1, 0);m2.Set(5, 2, 1);m2.Set(5, 3, 1);m2.Set(5, 4, 1);m2.Set(5, 5, 0);m2.Set(5, 6, 0);m2.Set(5, 7, 0);m2.Set(5, 8, 1);
    m2.Set(6, 0, 1);m2.Set(6, 1, 0);m2.Set(6, 2, 0);m2.Set(6, 3, 1);m2.Set(6, 4, 0);m2.Set(6, 5, 0);m2.Set(6, 6, 0);m2.Set(6, 7, 1);m2.Set(6, 8, 1);
    m2.Set(7, 0, 0);m2.Set(7, 1, 1);m2.Set(7, 2, 0);m2.Set(7, 3, 0);m2.Set(7, 4, 1);m2.Set(7, 5, 0);m2.Set(7, 6, 1);m2.Set(7, 7, 0);m2.Set(7, 8, 1);
    m2.Set(8, 0, 0);m2.Set(8, 1, 0);m2.Set(8, 2, 1);m2.Set(8, 3, 0);m2.Set(8, 4, 0);m2.Set(8, 5, 1);m2.Set(8, 6, 1);m2.Set(8, 7, 1);m2.Set(8, 8, 0);
    assert(n2 == m2);
}

void test_cint(void)
{
    real L = 32;
    std::vector<real> angles;
    std::vector<Vector> pos;
    Angdict angdict;
    //Posdict posdict;

    Vector a, b, c, d;
    a.x = 1.; a.y = 2.;
    b.x = 2.; b.y = 2.;
    c.x = 2.; c.y = 0.;
    d.x = 4.; d.y = 4.;

    angles = {0., PI/4., -PI/4., 0.};
    pos = {a, b, c, d};

    angdict[0] = angles;
    //posdict[0] = pos;

    Pi pi(angdict);

    /*pi:
    Vector(0., 0.)
    Vector(-0.70710678118, 0.)
    Vector(0.70710678118, 0.)
    Vector(0., 0.)
    */

    std::vector<real> ni(4);
    real nc;
    MATRIX n(4, 4);
    n.Zero();

    /* topological */

    calculate_n_topological(n, pos, L, 2);
    create_ni(n, ni);
    nc = average(ni);
    assert_almost_equal(nc, real(2));
    assert_almost_equal(Cint(0, pi, nc, n), real(-0.125));

    /* Metric */

    calculate_n_metric(n, pos, L, 2);
    create_ni(n, ni);
    nc = average(ni);
    assert_almost_equal(nc, real(1));
    assert_almost_equal(Cint(0, pi, nc, n), real(-0.25));

    /* Voronoi */

    calculate_n_voronoi(n, pos, L);
    create_ni(n, ni);
    nc = average(ni);
    assert_almost_equal(nc, real(3));
    assert_almost_equal(Cint(0, pi, nc, n), real(-1./12.));
}

void test_cts(void)
{
    real L = 32;
    std::vector<real> angles;
    std::vector<Vector> pos;
    Angdict angdict;
    //Posdict posdict;

    Vector a, b, c, d;
    a.x = 1.; a.y = 2.;
    b.x = 2.; b.y = 2.;
    c.x = 2.; c.y = 0.;
    d.x = 4.; d.y = 4.;

    angles = {0., PI/4., -PI/4., 0.};
    pos = {a, b, c, d};

    angdict[0] = angles;
    //posdict[0] = pos;

    Pi pi(angdict);

    /*pi:
    Vector(0., 0.)
    Vector(-0.70710678118, 0.)
    Vector(0.70710678118, 0.)
    Vector(0., 0.)
    */

    std::vector<real> ni(4);
    real nc;
    MATRIX n(4, 4);
    n.Zero();

    /* topological */

    calculate_n_topological(n, pos, L, 2);
    create_ni(n, ni);
    nc = average(ni);
    assert_almost_equal(nc, real(2));
    assert_almost_equal(CTs(0, pi, nc, ni), real(0.25));

    /* Metric */

    calculate_n_metric(n, pos, L, 2);
    create_ni(n, ni);
    nc = average(ni);
    assert_almost_equal(nc, real(1));
    assert_almost_equal(CTs(0, pi, nc, ni), real(0.375));

    /* Voronoi */

    calculate_n_voronoi(n, pos, L);
    create_ni(n, ni);
    nc = average(ni);
    assert_almost_equal(nc, real(3));
    assert_almost_equal(CTs(0, pi, nc, ni), real(0.25));
}

void test_average(void)
{
    //real average(const std::vector<real> &ni)
    std::vector<real> nums = {1., 2., 3., 4., 5., 6.};
    assert_almost_equal(average(nums), real(3.5));
}

void test_get_average_number_of_neighbours(void)
{
    MATRIX n(4, 4);
    n.Set(0, 0, 0); n.Set(0, 1, 1); n.Set(0, 2, 1); n.Set(0, 3, 1);
    n.Set(1, 0, 1); n.Set(1, 1, 0); n.Set(1, 2, 1); n.Set(1, 3, 1);
    n.Set(2, 0, 0); n.Set(2, 1, 1); n.Set(2, 2, 0); n.Set(2, 3, 1);
    n.Set(3, 0, 0); n.Set(3, 1, 1); n.Set(3, 2, 1); n.Set(3, 3, 0);

    std::vector<real> ni(4);
    create_ni(n, ni);

    assert_almost_equal(get_average_number_of_neighbours(n), average(ni));
}

void test_distance_periodic(void)
{
    //real distance_periodic(Vector p1, Vector p2)
    real L = 32;
    Vector a, b;

    a.x = 1;  a.y = 1;
    b.x = 31; b.y = 1;
    assert_almost_equal(distance_periodic(a, b, L), real(2));

    a.x = 1; a.y = 1;
    b.x = 1; b.y = 31;
    assert_almost_equal(distance_periodic(a, b, L), real(2));

    a.x = 1; a.y = 1;
    b.x = 31; b.y = 31;
    assert_almost_equal(distance_periodic(a, b, L), real(2)*std::sqrt(real(2)));

    a.x = 1; a.y = 1;
    b.x = 3; b.y = 1;
    assert_almost_equal(distance_periodic(a, b, L), real(2));

    a.x = 1; a.y = 1;
    b.x = 1; b.y = 3;
    assert_almost_equal(distance_periodic(a, b, L), real(2));

    a.x = 1; a.y = 1;
    b.x = 3; b.y = 3;
    assert_almost_equal(distance_periodic(a, b, L), real(2)*std::sqrt(real(2)));
}

void test_scalar_product(void)
{
    Vector a, b;
    a.x = 1.; a.y = 2.;
    b.x = 3.; b.y = 4.;
    assert_almost_equal(a*b, real(11));
}

void test_matrix_add(void)
{
    Vector a, b;
    a.x = 1; a.y = 2;
    b.x = 3; b.y = 4;
    Vector c;
    c.x = 4;
    c.y = 6;
    assert(a+b == c);
}

bool is_number_of_neighbours_equal(const MATRIX &n)
{
    real ns_old, ns_new;

    ns_old = get_number_of_neighbours(n, 0);

    for (unsigned int i = 1; i < n.GetWidth(); ++i)
    {
        ns_new = get_number_of_neighbours(n, i);
        if (ns_old != ns_new)
        {
            return false;
        }
    }
    return true;
}

bool is_symmetric(const MATRIX &n)
{
    unsigned int x, y;

    for (y = 0; y < n.GetHeight(); ++y)
    {
        for (x = 0; x < n.GetWidth(); ++x)
        {
            if (*n.Get(y, x) != *n.Get(x, y))
            {
                return false;
            }
        }
    }
    return true;
}

void test_maketimes(void)
{
    std::vector<unsigned int> maketimes, realtimes;

    maketimes = MakeTimes(1, 3, 2);
    realtimes = {1};
    assert(maketimes == realtimes);

    maketimes = MakeTimes(1, 4, 2);
    realtimes = {1, 3};
    assert(maketimes == realtimes);

    maketimes = MakeTimes(1, 5, 2);
    realtimes = {1, 3};
    assert(maketimes == realtimes);

}

void test_read_single(void)
{
    Posdict positions;
    Angdict angles;
    std::vector<unsigned int> times = {1, 2, 3};
    read_data_single("test.txt", positions, angles, times);
    assert(positions.size() == 3);
}


void test_read_consecutive(void)
{
    Posdict positions;
    Angdict angles;
    std::vector<unsigned int> times;

    times = {1, 2};
    read_data_consecutive("test.txt", positions, angles, times);
    assert(positions.size() == 3);
    assert(angles[1][0] == 1.);
    assert(angles[2][0] == 2.);
    assert(angles[3][0] == 3.);

    positions.clear();
    angles.clear();

    times = {1, 3};
    read_data_consecutive("test.txt", positions, angles, times);
    assert(positions.size() == 4);
    assert(angles[1][0] == 1.);
    assert(angles[2][0] == 2.);
    assert(angles[3][0] == 3.);
    assert(angles[4][0] == 4.);

    positions.clear();
    angles.clear();

    times = {1, 4};
    read_data_consecutive("test.txt", positions, angles, times);
    assert(positions.size() == 4);
    assert(angles[1][0] == 1.);
    assert(angles[2][0] == 2.);
    assert(angles[4][0] == 4.);
    assert(angles[5][0] == 5.);
}
