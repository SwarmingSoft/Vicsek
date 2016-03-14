//v14 19.08.2015

//Vicsek.exe datatest 0 0 1 0 linearvicsek 110000 256 16 0.01 1. 1. 0.12 1401 metric 0.45

//CARE: linearvicseknoisetransformation
//CARE!!!: eigen's .inverse() does not work with "-o3"
//confer: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=424

#include <random>
#include <cmath>
#include <sstream>

//only for linearvicseknoisetransform
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
//#include <unsupported/Eigen/MatrixFunctions>


#include "constants.hpp"
#include "tests.hpp"

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixReal;
typedef Eigen::Matrix<std::complex<real>, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixComplex;

template <typename T> findmin(T list)
{
     unsigned int pos = 0;
     typename T::Scalar elm = list(0);
     for (int i = 1; i < list.size(); ++i)
     {
          if (std::abs(elm) >= std::abs(list(i)))
          {
               elm = list(i);
               pos = i;
          }
     }
     //std::cout << elm << std::endl;
     return pos;
}

enum class Model {Vicsek, LinearVicsek, VicsekNoiseTransform, LinearVicsekNoiseTransform};

class CMDParameterGenerate
{

public:

    std::string out_file, out_extrasfile, out_neighboursfile;
    int out_opt;
    unsigned int out_skip; // 0
    unsigned int out_every; // 1
    bool out_consecutive; // false

    Model sim_model; //vicsek, linearvicsek, linearvicseknoisetransform
    unsigned int sim_timesteps; // 110000
    unsigned int sim_N; // 256
    real sim_L; // 16
    real sim_deltat; // 0.01
    real sim_v0; // 1.

    real sim_JV; // 1.
    real sim_eta; // 0.12 //care: linearvicsek's eta not same variance as vicsek's

    unsigned int sim_seed;

    Geometry sim_geometry; //metric, topological, voronoi
    real sim_param;        //r0    , neighbours , [gets discarded]

    /*param: r0 = 0.5 for metric; metric interaction radius //0.45 interpolation nc=6 for eta=0.3 paramter set
      param: NEIGHBOURS = 6 only for topological
    */

    CMDParameterGenerate(int argc, char** argv)
    {
        if (argc == 1)
        {
            std::cerr << "filename(string)\toutoption(0:data+extras,1:data,2:extras)\tSkip(int)\tEvery(int)\tConsecutive(int)\tModel(vicsek,linearvicsek,linearvicseknoisetransform)\tTimesteps(int)\tN(int)\tL(real)\tdeltat(real)\tv0(real)\tJV(real)\teta(real)\tseed(int)\tGeometry(metric,topological,voronoi)\tParam(real)" << std::endl;
            throw std::runtime_error("Parameters missing");
        }
        assert(argc == 1 + 16);

        out_file = std::string(argv[1]) + ".txt";
        out_extrasfile = std::string(argv[1]) + ".extras.txt";
        out_neighboursfile = std::string(argv[1]) + ".neighbours.txt";

        switch (std::atoi(argv[2]))
        {
        case 0: //7
            out_opt = OutOpt::Main|OutOpt::Extra|OutOpt::Neighbours;
            break;
        case 1:
            out_opt = OutOpt::Main;
            break;
        case 2:
            out_opt = OutOpt::Extra;
            break;
        case 3:
            out_opt = OutOpt::Main|OutOpt::Extra;
            break;
        case 4:
            out_opt = OutOpt::Neighbours;
            break;
        case 5:
            out_opt = OutOpt::Main|OutOpt::Neighbours;
            break;
        case 6:
            out_opt = OutOpt::Extra|OutOpt::Neighbours;
            break;
        case 7:
            out_opt = OutOpt::Main|OutOpt::Extra|OutOpt::Neighbours;
            break;
        default:
            throw std::runtime_error("Invalid output option parameter");
        }

        out_skip = std::atoi(argv[3]);
        out_every = std::atoi(argv[4]);
        out_consecutive = std::atoi(argv[5]);

        std::string _model = argv[6];

        if (_model == "vicsek")
        {
            sim_model = Model::Vicsek;
        }
        else if (_model == "linearvicsek")
        {
            sim_model = Model::LinearVicsek;
        }
        else if (_model == "linearvicseknoisetransform")
        {
            sim_model = Model::LinearVicsekNoiseTransform;
        }
        else
        {
            throw std::runtime_error("Invalid model parameter");
        }

        sim_timesteps = std::atoi(argv[7]);
        sim_N = std::atoi(argv[8]);
        sim_L = std::atof(argv[9]);
        sim_deltat = std::atof(argv[10]);
        sim_v0 = std::atof(argv[11]);

        sim_JV = std::atof(argv[12]);
        sim_eta = std::atof(argv[13]);

        sim_seed = std::atoi(argv[14]);

        std::string _geometry = argv[15];

        if (_geometry == "topological")
        {
            sim_geometry = Geometry::topological;
        }
        else if (_geometry == "metric")
        {
            sim_geometry = Geometry::metric;
        }
        else if (_geometry == "voronoi")
        {
            sim_geometry = Geometry::voronoi;
        }
        else
        {
            throw std::runtime_error("Invalid geometry parameter");
        }
        sim_param = std::atof(argv[16]);
    }

    std::string header(void)
    {
        std::stringstream ss;

        ss << "Skip: " << out_skip << " Every: " << out_every << " Consecutive: " << out_consecutive << " Model: " << int(sim_model) << " Timesteps: " << sim_timesteps << " N: " << sim_N << " L: " << sim_L << " deltat: " << sim_deltat << " v0: " << sim_v0 << " JV: " << sim_JV << " eta: " << sim_eta << " seed: " << sim_seed << " Geometry: " << int(sim_geometry) << " Param: " << sim_param;

        return ss.str();
    }

};

bool test_output(const CMDParameterGenerate &args, unsigned int step)
{
    return step >= args.out_skip && (step % args.out_every == 0 || (step % args.out_every == 1 && args.out_consecutive));
}

void OutputNeighbours(std::ostream &outneighbours, const CMDParameterGenerate &args, const std::vector<real> &ni)
{
    if (args.out_opt & OutOpt::Neighbours)
    {
        for (unsigned int i = 0; i < ni.size()-1; ++i)
        {
            outneighbours << ni[i] << " ";
        }
        outneighbours << ni.back() << std::endl;
    }
}

void OutputPositionsAndAngles(std::ostream &out, std::ostream &outextras, const CMDParameterGenerate &args, unsigned int step, const std::vector<Vector> &positions, const std::vector<real> &angles, const MATRIX &n, MATRIX &n_old, real nc)
{
    unsigned int i;
    if (test_output(args, step))
    {
        if (args.out_opt & OutOpt::Main)
        {
            for (i = 0; i < positions.size()-1; ++i)
            {
                out << positions[i].x << " " << positions[i].y << " ";
            }
            out << positions.back().x << " " << positions.back().y << "\n";
            for (i = 0; i < angles.size()-1; ++i)
            {
                out << angles[i] << " ";
            }
            out << angles.back() << "\n";
        }
        if (args.out_opt & OutOpt::Extra)
        {
            outextras  << nc << " " << mixing_parameter(args.sim_deltat, n, n_old, nc) << " " << order_parameter_angle(angles) << std::endl;
        }
    }
}

void OutputPositionsAndAngleV(std::ostream &out, std::ostream &outextras, const CMDParameterGenerate &args, unsigned int step, const std::vector<Vector> &positions, const std::vector<Vector> &st1, const MATRIX &n, const MATRIX &n_old, real nc)
{
    unsigned int i;
    if (test_output(args, step))
    {
        if (args.out_opt & OutOpt::Main)
        {
            for (i = 0; i < positions.size()-1; ++i)
            {
                out << positions[i].x << " " << positions[i].y << " ";
            }
            out << positions.back().x << " " << positions.back().y << "\n";
            for (i = 0; i < st1.size()-1; ++i)
            {
                out << std::atan2(st1[i].y, st1[i].x) << " ";
            }
            out << std::atan2(st1.back().y, st1.back().x) << "\n";
        }
        if (args.out_opt & OutOpt::Extra)
        {
            outextras  << nc << " " << mixing_parameter(args.sim_deltat, n_old, n, nc) << " " << order_parameter(st1) << std::endl;
        }
    }
}

void FixPeriodicBoundaries(std::vector<Vector> &positions, unsigned int i, real L)
{
    //for periodic boundary conditions
    if (positions[i].x < 0.) { positions[i].x += L; }
    else if (positions[i].x >= L) { positions[i].x -= L; }
    if (positions[i].y < 0.) { positions[i].y += L; }
    else if (positions[i].y >= L) { positions[i].y -= L; }

    assert(positions[i].x >= 0 && positions[i].x <= L);
    assert(positions[i].y >= 0 && positions[i].y <= L);
}

int main(int argc, char** argv)
{
    CMDParameterGenerate args(argc, argv);

    std::ofstream out, outextras, outneighbours;
    unsigned int i, j, step;

    //write parameters in first line of file
    if (args.out_opt & OutOpt::Main)
    {
        out.open(args.out_file);
        //out.precision(4);
        out << args.header() << std::endl;
    }
    if (args.out_opt & OutOpt::Extra)
    {
        outextras.open(args.out_extrasfile);
        outextras << args.header() << std::endl;
        outextras << "nc" << " " << "mixing_parameter" << " " << "order_parameter" << std::endl;
    }
    if (args.out_opt & OutOpt::Neighbours)
    {
        outneighbours.open(args.out_neighboursfile);
        outneighbours << args.header() << std::endl;
    }

    MATRIX n(args.sim_N, args.sim_N); //neighborhood matrix
    n.Zero();

    MATRIX n_old(args.sim_N, args.sim_N); //neighborhood matrix
    n_old.Zero();

    real sqrtdeltat = std::sqrt(args.sim_deltat);
    std::default_random_engine gen(args.sim_seed);

    real nc = 0;

    std::vector<Vector> positions(args.sim_N);

    std::cout << "Initialization complete" << std::endl;

    if (args.sim_model == Model::Vicsek || args.sim_model == Model::VicsekNoiseTransform)
    {
        std::vector<real> angles(args.sim_N);

        std::uniform_real_distribution<real> distn(-args.sim_eta*PI, args.sim_eta*PI);
        std::uniform_real_distribution<real> dista(-PI, PI);
        std::uniform_real_distribution<real> distp(0, args.sim_L);

        for (i = 0; i < positions.size(); ++i)
        {
            angles[i] = dista(gen);
            positions[i] = Vector(distp(gen), distp(gen));
        }

        Vector w;
        std::vector<Vector> xy(args.sim_N);

        std::vector<real> ni(args.sim_N);
        for (step = 0; step < args.sim_timesteps; ++step)
        {
            switch (args.sim_geometry)
            {
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions, args.sim_L);
                    assert(is_symmetric(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions, args.sim_L, args.sim_param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::topological:
                    calculate_n_topological(n, positions, args.sim_L, int(args.sim_param));
                    assert(is_number_of_neighbours_equal(n));
                    break;
            }

            create_ni(n, ni);
            nc = average(ni); //nc = get_average_number_of_neighbours(n);

            for (i = 0; i < n.GetWidth(); ++i)
            {
                w.x = 0.;
                w.y = 0.;
                for (j = 0; j < n.GetHeight(); ++j)
                {
                    if (*n.Get(i, j) != 0)
                    {
                        w += *n.Get(i, j)*Vector(std::cos(angles[j]), std::sin(angles[j]));
                    }
                }
                xy[i] = Vector(std::cos(angles[i]), std::sin(angles[i])) + w*args.sim_JV*args.sim_deltat;
            }

            for (i = 0; i < n.GetWidth(); ++i)
            {
                angles[i] = std::atan2(xy[i].y, xy[i].x) + sqrtdeltat*distn(gen);
                positions[i] += Vector(std::cos(angles[i]), std::sin(angles[i])) * args.sim_v0*args.sim_deltat;

                FixPeriodicBoundaries(positions, i, args.sim_L);
            }

            OutputPositionsAndAngles(out, outextras, args, step, positions, angles, n, n_old, nc);
            OutputNeighbours(outneighbours, args, ni);
            n_old = n;
        }
    }
    else if (args.sim_model == Model::LinearVicsek || args.sim_model == Model::LinearVicsekNoiseTransform)
    {
        std::vector<Vector> st0(args.sim_N); //s_i^t im paper
        std::vector<Vector> st1(args.sim_N); //s_i^{t+1} im paper

        std::vector<Vector> pit0(args.sim_N); //pi_i^t im paper
        std::vector<Vector> pit1(args.sim_N); //pi_i^{t+1} im paper

        std::vector<Vector> slt0(args.sim_N); //sl_i^t im paper
        std::vector<Vector> slt1(args.sim_N); //sl_i^{t+1} im paper

        std::uniform_real_distribution<real> distn(-args.sim_eta, args.sim_eta);
        std::uniform_real_distribution<real> dista(-PI, PI);
        std::uniform_real_distribution<real> distp(0, args.sim_L);

        real angle;
        for (i = 0; i < positions.size(); ++i)
        {
            angle = dista(gen);
            st0[i] = Vector(std::cos(angle), std::sin(angle));
            positions[i] = Vector(distp(gen), distp(gen));
        }

        std::vector<real> ni(args.sim_N);
        for (step = 0; step < args.sim_timesteps; ++step)
        {
            switch (args.sim_geometry)
            {
                case Geometry::voronoi:
                    calculate_n_voronoi(n, positions, args.sim_L);
                    assert(is_symmetric(n));
                    break;
                case Geometry::metric:
                    calculate_n_metric(n, positions, args.sim_L, args.sim_param);
                    assert(is_symmetric(n));
                    break;
                case Geometry::topological:
                    calculate_n_topological(n, positions, args.sim_L, int(args.sim_param));
                    assert(is_number_of_neighbours_equal(n));
                    break;
            }

            create_ni(n, ni);
            nc = average(ni); //nc = get_average_number_of_neighbours(n);

            // average flight direction n
            Vector normn(0, 0);
            for (i = 0; i < n.GetWidth(); ++i)
            {
                normn += st0[i];
            }
            normn.Normalize();

            //pit0 as orthogonal projection of st0 on n
            for (i = 0; i < n.GetWidth(); ++i)
            {
                slt0[i] = normn*(st0[i]*normn);
                pit0[i] = st0[i] - slt0[i];
            }

            //real newJ = 0.5;
            real newJ = args.sim_deltat*args.sim_JV/(1.+args.sim_deltat*args.sim_JV*nc);

            //sum over neighbourhoods
            Vector sum_pit0, sum_slt0;
            for (i = 0; i < n.GetWidth(); ++i)
            {
                sum_pit0.x = 0; sum_pit0.y = 0;
                sum_slt0.x = 0; sum_slt0.y = 0;

                for (j = 0; j < n.GetWidth(); ++j)
                {
                    if (*n.Get(i, j) != 0) // optional, but maybe faster
                    {
                        sum_pit0 += *n.Get(i, j)*pit0[j];
                        sum_slt0 += *n.Get(i, j)*slt0[j];
                    }
                }

                // D6 without noise
                pit1[i] = (1. - newJ*ni[i]) * pit0[i] + newJ*sum_pit0;
                slt1[i] = (1. - newJ*ni[i]) * slt0[i] + newJ*sum_slt0;
                assert_almost_equal(pit1[i]*slt1[i], real(0));
                //std::cout << i << ": " << pit1[i].x << ", " << pit1[i].y << "\n";
            }

            //std::cout << pit0[0].x << "," << pit0[0].y << " | " << sum_pit0.x << "," << sum_pit0.y << std::endl;

            //add noise only to pit1
            Vector pitabs = pit1[0] / pit1[0].Norm(); // test if all pi do indeed point in the same direction
            //Vector pitabs1 = pit1[1] / pit1[1].Norm(); // test if all pi do indeed point in the same direction
            //std::cout << pitabs.x << "," << pitabs.y << " | " << pitabs1.x << "," << pitabs1.y << std::endl;


            if (args.sim_model == Model::LinearVicsekNoiseTransform)
            {
                //noise; still scalar
                EigenMatrixComplex noise(args.sim_N, 1);
                for (i = 0; i < n.GetWidth(); ++i)
                {
                    noise(i) = distn(gen);
                }

                // Eigenvectors
                EigenMatrixReal lambda(args.sim_N, args.sim_N);
                for (i = 0; i < n.GetWidth(); ++i)
                {
                    for (j = 0; j < n.GetHeight(); ++j)
                    {
                        lambda(i, j) = ni[i]*KroneckerDelta(i, j) - real(*n.Get(i, j));
                    }
                }

                assert(is_weakly_diagonally_dominant_and_non_negativ_diagonal_entries(lambda));

                Eigen::EigenSolver<EigenMatrixReal> es;
                es.compute(lambda, /* computeEigenvectors = */ Eigen::ComputeEigenvectors);
                assert(es.info() == Eigen::ComputationInfo::Success);

                auto ev_vector = es.eigenvalues();
                EigenMatrixComplex trafo = es.eigenvectors(); // consists of eigen vectors

                /*std::cout << noise(0) << std::endl;
                real noisesum = 0.;
                for (i = 0; i < n.GetWidth(); ++i)
                {
                    noisesum += noise(i).real();
                }
                std::cout << noisesum << std::endl;*/

                // Transform
                noise = trafo.inverse() * noise; //trafo or trafo.inverse() ???

                // Cancel noise
                //std::cout << ev_vector << std::endl;
                int ewPos = findmin(ev_vector); //find position of eigenvalue == 0
                noise(ewPos) = 0.;

                /*std::cout << "ewPos: " << ewPos << std::endl;
                EigenMatrixComplex trafoinv = trafo.inverse();
                for (i = 0; i < n.GetWidth(); ++i)
                {
                    std::cout << trafoinv(i, ewPos) << std::endl;
                    std::cout << trafoinv(ewPos, i) << std::endl;
                }*/

                // Transform back
                noise = trafo * noise;

                /*std::cout << noise(0) << std::endl;
                noisesum = 0;
                for (i = 0; i < n.GetWidth(); ++i)
                {
                    noisesum += noise(i).real();
                }
                std::cout << noisesum << std::endl;*/

                for (i = 0; i < n.GetWidth(); ++i)
                {
                    //std::cout << noise(i).imag() << std::endl;
                    assert(std::abs(noise(i).imag()) < 0.01);
                    pit1[i] = pit1[i] + sqrtdeltat * noise(i).real() * pitabs;
                }
            }
            else // normal linear model
            {
                for (i = 0; i < n.GetWidth(); ++i)
                {
                    pit1[i] = pit1[i] + sqrtdeltat * distn(gen) * pitabs;
                }
            }


            //real sl, pit1norm;
            //real st1avg = 0;
            for (i = 0; i < n.GetWidth(); ++i)
            {
                /*pit1norm = pit1[i].Norm();
                if (pit1norm > 1.)
                {
                    sl = 0.;
                    pit1[i] /= pit1[i].Norm();
                }
                else
                {
                    sl = std::sqrt(1 - pit1norm*pit1norm);
                }
                slt1[i] = sl*normn;*/ //old n??
                st1[i] = pit1[i] + slt1[i];

                /*if (st1[i].Norm() > 1.01)
                {
                    std::cout << "st1: " << st1[i].Norm() << "\n";
                }*/

                //DESTROYS TRUE MODEL: normalize st1
                st1[i] = st1[i]/st1[i].Norm();
                //st1avg += st1[i].Norm();

                //Vector slt1alt = std::sqrt(1-pit1[i]*pit1[i])*normn;
                //std::cout << i << ":  "<< pit1[i]*pit1[i] << "  " << slt1[i].x << "," << slt1[i].y << "  " << slt1alt.x << "," << slt1alt.y << "->"<< (slt1alt/slt1alt.Norm())*(slt1[i]/slt1[i].Norm()) << "\n";
                //std::cout << st1[i].Norm() << std::endl;
            }
            //std::cout << st1avg/st1.size() << std::endl;

            for (i = 0; i < n.GetWidth(); ++i)
            {
                // D4
                positions[i] += args.sim_v0*args.sim_deltat*st1[i];

                FixPeriodicBoundaries(positions, i, args.sim_L);
            }

            OutputPositionsAndAngleV(out, outextras, args, step, positions, st1, n, n_old, nc);
            OutputNeighbours(outneighbours, args, ni);
            n_old = n;
            st0 = st1;
        }
    }

    std::cout << "Evaluation complete" << std::endl;

    if (args.out_opt & OutOpt::Main)
    {
        out.close();
    }
    if (args.out_opt & OutOpt::Extra)
    {
        outextras.close();
    }
    if (args.out_opt & OutOpt::Neighbours)
    {
        outneighbours.close();
    }

    return EXIT_SUCCESS;
}
