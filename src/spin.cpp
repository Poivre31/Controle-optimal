#include <solver/coordinates.h>
#include <solver/solver.h>
#include <math/gradient_descent.h>
#include <matplot/matplot.h>
#include <functional>
#include <chrono>
#include <iostream>
#include <Eigen/Core>
#include <complex>
#include <Eigen/Eigen>

using C = Eigen::Vector<std::complex<double>, 2>;
using namespace std::complex_literals;
namespace plt = matplot;

enum
{
    c1,
    c2,
};

double dt = 0.1;
double tf = 1;

C propagate(double, C Cj, double Bj)
{
    Eigen::Matrix2d sig_x{{0., 1.}, {1., 0.}};
    C Cj_1 = Cj;
    Cj_1[0] += dt * 1i * Bj * Cj[1] * .5;
    std::cout << "1 " << Cj_1[0] << Cj_1[1] << " " << dt * 1i * Bj * Cj[1] * .5 << std::endl;
    Cj_1[1] += dt * 1i * Bj * Cj[0] * .5;
    std::cout << "2 " << Cj_1[0] << Cj_1[1] << std::endl;

    return Cj_1;
}

int main()
{
    solver_degree_I<C> solver;
    solver.set_initial_state(0, C(1 / sqrt(2), 1 / sqrt(2)));
    auto a = C(1 / sqrt(2), 1 / sqrt(2));
    solver.set_timestep(dt);
    size_t N = solver.get_timeline(tf).size();
    std::vector<double> B_t(N);
    for (size_t i = 0; i < N; i++)
    {
        B_t[i] = 0.1;
    }

    auto dcdt = [B_t](double t, C Cj)
    {
        return propagate(0, Cj, B_t[round(t * dt)]);
    };

    solver.solve_euler(tf, dcdt);
    auto result = solver.get_positions();

    std::vector<double> fidelity(N);
    for (size_t i = 0; i < N; i++)
    {
        fidelity[i] = pow(abs(result[i][0]), 2);
    }
    plt::plot(solver.get_timeline(tf), fidelity);
    plt::show();
}