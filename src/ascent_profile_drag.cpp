#include <solver/coordinates.h>
#include <solver/solver.h>
#include <math/gradient_descent.h>
#include <matplot/matplot.h>
#include <functional>
#include <chrono>
#include <iostream>

using namespace std::chrono;

namespace plt = matplot;

using state = coordinates<7>;
using parameters = coordinates<3>;

enum
{
    X,
    Z,
    Vx,
    Vz,
    λ2,
    λ3,
    λ4
};

enum
{
    λ2_0,
    λ3_0,
    λ4_0
};

double dt = 0.01;

double M = 10e3;
double k = 20 * M;
double A = k / 3000;
double tf = 0.95 * M / A;

double g = 9.81;
double Cd = 0.2;
double H = 3000;

double η(double z, double m)
{
    return Cd * exp(-z / H) / m;
}

double dηdz(double z, double m)
{
    return -η(z, m) / H;
}

state dsdt(double t, state state_i, parameters c)
{

    state state_change;
    double m;

    m = M - t * A;

    double vx = state_i[Vx];
    double vz = state_i[Vz];
    double v = sqrt(vx * vx + vz * vz);
    double z = state_i[Z];

    double α = M_PI_2;
    if (state_i[λ4] != 0)
        α = atan(state_i[λ3] / state_i[λ4]);
    if (α < 0)
        α += M_PI;
    if (std::isnan(α))
    {
        std::cout << state_i << std::endl;
        exit(1);
    }

    state_change[X] = vx;
    state_change[Z] = vz;
    state_change[Vx] = k / m * sin(α) - η(z, m) * vx * v;
    state_change[Vz] = k / m * cos(α) - η(z, m) * vz * v - g;

    if (v != 0)
    {
        state_change[λ2] = dηdz(z, m) * v * (state_i[λ3] * vx + state_i[λ4] * vz);
        state_change[λ3] = η(z, m) * v * (state_i[λ3] * (1 + vx * vx / (v * v)) + state_i[λ4] * vx * vz / (v * v));
        state_change[λ4] = -state_i[λ2] + η(z, m) * v * (state_i[λ4] * (1 + vz * vz / (v * v)) + state_i[λ3] * vx * vz / (v * v));
    }
    else
    {
        state_change[λ2] = 0;
        state_change[λ3] = 0;
        state_change[λ4] = -state_i[λ2];
    }

    return state_change;
}

std::vector<state> propagate(solver_degree_I<state> &solver, parameters c)
{
    state X0 = {0., 0., 0., 0., c[λ2_0], c[λ3_0], c[λ4_0]};

    solver.set_initial_state(0, X0);
    solver.set_timestep(dt);
    solver.solve_RK4(tf, [c](double t, state state_i)
                     { return dsdt(t, state_i, c); });
    return solver.get_positions();
}

matrix<3, 4> AA{
    {0., 2., 1., -2.},
    {4., -7., 0., 0.},
    {3., 9., 0., 1.}};

coordinates<4> f(coordinates<3> x)
{
    return {x[0] * x[1], x[0] * x[1] * x[2], cos(x[2]), sin(x[1])};
}

int main()
{
    auto t1 = high_resolution_clock::now();
    solver_degree_I<state> solver;

    parameters c = {1., 12., 60.};

    size_t outside = 0;
    auto loss = [&solver, &outside](parameters c)
    {
        auto result = propagate(solver, c).back();
        outside++;
        if (outside == 2 * 2 * 10)
        {
            fmt::println("Z: {:.5g}, Vx: {:.5g}, Vz: {:.5g}, c2: {:.5g}, c3: {:.5g}, c4: {:.5g}", result[Z], result[Vx], result[Vz], c[λ2_0], c[λ3_0], c[λ4_0]);
            outside = 0;
        }
        double target = 100000;
        double loss = (result[Z] / target - 1) * (result[Z] / target - 1) + result[Vz] * result[Vz] / (2000 * 2000);
        return loss;
    };

    c = gradient_descent<3>(c, loss, 10000, 0, 1e-12);
    std::cout << c << std::endl;

    auto result = stride(solver.get_positions(), 100); ///> Stride to visualize only one tenth of datapoints to reduce load on matplot
    auto time = stride(solver.get_timeline(tf), 100);

    auto t2 = high_resolution_clock::now();
    double runtime = duration_cast<milliseconds>(t2 - t1).count();
    fmt::println("Runtime: {}ms", runtime);
    fmt::println("Final height: {:.1f}m, velocity: ({:.1f}, {:.1f})m/s, time: {:.1f}s", parse(result, Z).back(), parse(result, Vx).back(), parse(result, Vz).back(), time.back());
    fmt::println("Final covector: {:.3g}, {:.3g}, {:.3g}", parse(result, λ2).back(), parse(result, λ3).back(), parse(result, λ4).back());

    if (true)
    {
        auto fig = plt::figure();
        fig->size(1920, 1080);
        plt::subplot(2, 2, 0);
        plt::plot(parse(result, X), parse(result, Z));
        plt::title("Trajectory z(x)");

        plt::subplot(2, 2, 1);
        plt::plot(time, apply_element_wise<state, double>(result, [c](state state_i)
                                                          {
        double α = M_PI_2;
        if (state_i[λ4] != 0)
            α = atan(state_i[λ3] / state_i[λ4]);
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
