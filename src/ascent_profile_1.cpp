#include <solver/coordinates.h>
#include <solver/solver.h>
#include <math/gradient_descent.h>
#include <matplot/matplot.h>
#include <functional>
#include <chrono>
#include <iostream>

using namespace std::chrono;

namespace plt = matplot;

using state = coordinates<4>;
using parameters = coordinates<2>;

enum
{
    X,
    Z,
    Vx,
    Vz,
};

enum
{
    C1,
    C2,
};

double dt = 0.01;

double M = 10e3;
double k = 20 * M;
double A = k / 3000;
double tf = 0.9 * M / A;

double g = 9.81;

state dsdt(double t, state state_i, parameters c)
{

    state state_change;
    double m;

    m = M - t * A;

    double vx = state_i[Vx];
    double vz = state_i[Vz];

    double α = M_PI_2;
    if (c[C2] != t)
        α = atan(c[C1] / (c[C2] - t));
    if (α < 0)
        α += M_PI;

    state_change[X] = vx;
    state_change[Z] = vz;
    state_change[Vx] = k / m * sin(α);
    state_change[Vz] = k / m * cos(α) - g;

    return state_change;
}

std::vector<state> propagate(solver_degree_I<state> &solver, parameters c)
{
    state X0 = {0., 0., 0., 0.};

    solver.set_initial_state(0, X0);
    solver.set_timestep(dt);
    solver.solve_RK4(tf, [c](double t, state state_i)
                     { return dsdt(t, state_i, c); });
    return solver.get_positions();
}

int main()
{
    auto t1 = high_resolution_clock::now();
    solver_degree_I<state> solver;

    parameters c = {1, tf};

    size_t outside = 0;
    auto loss = [&solver, &outside](parameters c)
    {
        auto result = propagate(solver, c).back();
        outside++;
        if (outside == 2 * 2)
        {
            fmt::println("Z: {:.5g}, Vx: {:.5g}, Vz: {:.5g}, c1: {:.3g}, c2: {:.3g}", result[Z], result[Vx], result[Vz], c[C1], c[C2]);
            outside = 0;
        }
        double target = 50000;
        double loss = (result[Z] / target - 1) * (result[Z] / target - 1) + result[Vz] * result[Vz] / (8000 * 8000);
        return loss;
    };

    c = gradient_descent<2>(c, loss, 1000, 0, 1e-10);
    std::cout << c << std::endl;

    auto result = stride(solver.get_positions(), 100); ///> Stride to visualize only one tenth of datapoints to reduce load on matplot
    auto time = stride(solver.get_timeline(tf), 100);

    auto t2 = high_resolution_clock::now();
    double runtime = duration_cast<milliseconds>(t2 - t1).count();
    fmt::println("Runtime: {}ms", runtime);
    fmt::println("Final height: {:.1f}m, velocity: ({:.1f}, {:.1f})m/s, time: {:.1f}s", parse(result, Z).back(), parse(result, Vx).back(), parse(result, Vz).back(), time.back());

    if (true)
    {
        auto fig = plt::figure();
        fig->size(1920, 1080);
        plt::subplot(2, 2, 0);
        plt::plot(parse(result, X), parse(result, Z));
        plt::title("Trajectory z(x)");

        plt::subplot(2, 2, 1);
        plt::plot(time, apply_element_wise<double, double>(time, [c](double t)
                                                           {
            double α = M_PI_2;
            if (c[C2] != t)
                α = atan(c[C1] / (c[C2] - t));
            if (α < 0)
                α += M_PI;
            return α; }));

        plt::title("Angle α(t)");

        plt::subplot(2, 2, 2);
        plt::plot(time, parse(result, Z));
        plt::title("Alitutde z(t)");

        plt::subplot(2, 2, 3);
        plt::plot(time, parse(result, Vx));
        plt::title("Ortho velocity vx(t)");

        plt::show();
    }
}
