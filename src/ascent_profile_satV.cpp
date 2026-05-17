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
using parameters = coordinates<8>;

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
    C3,
    C4,
    C1b,
    C2b,
    C3b,
    C4b
};

double dt = 1;

double M_1 = 3e6;
double A_1 = 5 * 2661;
double k_1 = 2596 * A_1;
double t_1 = 159;

double M_2 = 8e5;
double A_2 = 5 * 250;
double k_2 = 4130 * A_2;
double t_2 = 395;

double tf = t_1 + t_2;

double g = 9.81;
double γ = 0.0e-5;
double H = 8000.;

state dsdt(double t, state state_i, parameters c)
{

    state state_change;
    double m;
    double k;
    if (t < t_1)
    {
        m = M_1 - t * A_1;
        k = k_1;
    }
    else
    {
        m = M_2 - t * A_2;
        k = k_2;
    }
    double drag = γ * exp(-state_i[Z] / H);
    double vx = state_i[Vx];
    double vz = state_i[Vz];
    double v = sqrt(vx * vx + vz * vz);

    double α = M_PI_2;
    if (t < t_1)
    {
        if (c[C4] != c[C2] * t)
            α = atan((0 - c[C1] * t) / (c[C4] - c[C2] * t));
        if (α < 0)
            α += M_PI;
    }
    else
    {
        if (c[C4b] != c[C2b] * t)
            α = atan((c[C3b] - c[C1b] * t) / (c[C4b] - c[C2b] * t));
        if (α < 0)
            α += M_PI;
    }

    state_change[X] = vx;
    state_change[Z] = vz;
    state_change[Vx] = k / m * sin(α) - drag * vx * v;
    state_change[Vz] = k / m * cos(α) - drag * vz * v - g;

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

    parameters c = {-1, 1, 0., t_1, -10, 1, 0, tf};

    size_t outside = 0;
    auto loss = [&solver, &outside](parameters c)
    {
        auto result = propagate(solver, c).back();
        outside++;
        if (outside == 100 * 2)
        {
            fmt::println("Z: {:.5g}, Vx: {:.5g}, Vz: {:.5g}", result[Z], result[Vx], result[Vz]);
            outside = 0;
        }
        double target = 185000;
        return (result[Z] / target - 1) * (result[Z] / target - 1) + result[Vz] * result[Vz] / 4000000 + (result[Vx] / 7815 - 1) * (result[Vx] / 7815 - 1);
    };

    c = gradient_descent<8>(c, loss, 100000, 0, 1e-10);
    std::cout << c << std::endl;

    auto result = solver.get_positions();
    auto time = solver.get_timeline(tf);

    auto t2 = high_resolution_clock::now();
    double runtime = duration_cast<milliseconds>(t2 - t1).count();
    fmt::println("Runtime: {}ms", runtime);
    fmt::println("Final height: {:.1f}m, velocity: ({:.1f}, {:.1f})m/s, time: {:.1f}s", parse(result, Z).back(), parse(result, Vx).back(), parse(result, Vz).back(), time.back());
    fmt::println("Final horizontal distance: {:.1f}km", parse(result, X).back() / 1000);

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
            double α;
            if (t < t_1)
            {
                if (c[C4] != c[C2] * t)
                    α = atan((0 - c[C1] * t) / (c[C4] - c[C2] * t));
                if (α < 0)
                    α += M_PI;
            }
            else
            {
                if (c[C4b] != c[C2b] * t)
                    α = atan((c[C3b] - c[C1b] * t) / (c[C4b] - c[C2b] * t));
                if (α < 0)
                    α += M_PI;
            }
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
