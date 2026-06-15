#include <solver/coordinates.h>
#include <solver/solver.h>
#include <math/gradient_descent.h>
#include <matplot/matplot.h>
#include <functional>
#include <chrono>
#include <iostream>

/// PROBLEM CONDITIONS

double dt = 0.01;
double g = 9.81;

double target_altitude = 185000;
size_t max_steps = 1000;
double target_precision = 1e-15;

double isp = 3000;
double Mf = 140000; ///> Initial mass of the rocket
// double Mf = Mi / 20;        ///> Final mass, at end of burn
double Mi = Mf * exp(8000 / isp);
double thrust = 1.54 * g * Mi; ///> Thrust (here initial TWR of 2)
double A = thrust / isp;       ///> Fuel consumption rate (here with an ISP of 3000m/s)
double tf = (Mi - Mf) / A;     ///> Burn time

/// PRO

using namespace std::chrono;
namespace plt = matplot;

using state = coordinates<4>; ///> We define our state vector...
enum
{
    X,
    Z,
    Vx,
    Vz,
};

using parameters = coordinates<2>; ///> And our parameters
enum
{
    C1,
    C2,
};

/// @brief Here we define the evolution function of our state variables X,Y,Vx,Vz according to Newton's 2nd law, using the
/// current estimate of the parameters c.
state dsdt(double t, state state_i, parameters c)
{

    state state_change;

    double M = Mi - t * A;

    double α = M_PI_2;
    if (c[C2] != t)
        α = atan(c[C1] / (c[C2] - t));
    if (α < 0)
        α += M_PI;

    state_change[X] = state_i[Vx];
    state_change[Z] = state_i[Vz];
    state_change[Vx] = thrust / M * sin(α);
    state_change[Vz] = thrust / M * cos(α) - g;

    return state_change;
}

/// @brief Here we compute the final state for a given set of parameters using successive Runge-Kutta iterations.
/// @return A vector containing all the intermediate states.
std::vector<state> propagate(solver_degree_I<state> &solver, parameters c)
{
    solver.set_initial_state(0, {0., 0., 0., 0.});
    solver.set_timestep(dt);
    solver.solve_RK4(tf, [c](double t, state state_i)
                     { return dsdt(t, state_i, c); });
    return solver.get_positions();
}

int main()
{
    ///> Highest reachable altitude, with the rocket pointing straight up for the duration of the burn. Target altitude can't be higher than this.
    double max_reachable_altitude = thrust * Mi / (A * A) * (1 - A * tf / Mi) * log(1 - A * tf / Mi) + thrust * tf / A - g * tf * tf / 2;
    fmt::println("Highest reachable altitude {}", max_reachable_altitude);
    if (target_altitude > max_reachable_altitude)
        fmt::println("WARNING: target altitude {} is higher than the highest reachable altitude {}", target_altitude, max_reachable_altitude);
    auto t1 = high_resolution_clock::now();

    solver_degree_I<state> solver;

    parameters c = {1, tf};

    auto loss = [&solver](parameters c)
    {
        state result = propagate(solver, c).back();
        double loss = (result[Z] / target_altitude - 1) * (result[Z] / target_altitude - 1) + result[Vz] * result[Vz] / (3000 * 3000);
        return loss;
    };

    c = gradient_descent_BB<2>(c, loss, max_steps, target_precision, 100);

    auto result = stride(solver.get_positions(), 100); ///> To only visualise a portion of the points, reduces load on matplot
    auto time = stride(solver.get_timeline(tf), 100);

    auto t2 = high_resolution_clock::now();
    double runtime = duration_cast<milliseconds>(t2 - t1).count();
    fmt::println("Runtime: {}ms", runtime);
    fmt::println("Final height: {:.1f}m, velocity: ({:.1f}, {:.1f})m/s, burn time: {:.1f}s", parse(result, Z).back(), parse(result, Vx).back(), parse(result, Vz).back(), time.back());

    if (true)
    {
        auto fig = plt::figure();
        fig->size(1920, 1080);
        plt::subplot(2, 2, 0);
        plt::plot(parse(result, X), parse(result, Z));
        plt::title("Trajectoire Z(X)");

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
        plt::title("Altitude Z(t)");

        plt::subplot(2, 2, 3);
        plt::plot(time, parse(result, Vx));
        plt::title("Vitesse horizontale Vx(t)");

        plt::show();
    }
}
