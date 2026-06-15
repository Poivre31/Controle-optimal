#pragma once
#include <Eigen/Core>
#include <fmt/format.h>
#include "euler.h"
#include "midpoint.h"
#include "RK4.h"
#include "solver/coordinates.h"
#include <iostream>
#include "solver/verlet.h"
#include "solver/yoshida.h"

size_t get_number_of_steps(double t0, double tf, double dt)
{
  if (tf <= t0)
    fmt::println("ERROR: end time must be greater than initial time");
  return (size_t)round((tf - t0) / dt);
}

template <typename T>
class solver_degree_I
{
public:
  solver_degree_I() {};

  void set_initial_state(double t0, T x0)
  {
    _t0 = t0;
    _x0 = x0;
  }

  void set_timestep(double dt) { _dt = dt; }

  T solve_euler(double tf, std::function<T(double t, T x)> dxdt)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    _positions = euler_explicit(_t0, tf, N, _x0, dxdt);
    return _positions[N];
  }

  T solve_midpoint(double tf, std::function<T(double t, T x)> dxdt)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    _positions = midpoint(_t0, tf, N, _x0, dxdt);
    return _positions[N];
  }

  T solve_RK4(double tf, std::function<T(double t, T x)> dxdt)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    _positions = RK4_explicit(_t0, tf, N, _x0, dxdt);
    return _positions[N];
  }

  std::vector<T> get_positions() { return _positions; }

  std::vector<double> get_timeline(double tf)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::vector<double> t(N + 1);
    for (size_t i = 0; i < N + 1; i++)
    {
      t[i] = _dt * i;
    }
    return t;
  }

private:
  double _t0;
  double _dt;
  T _x0;
  std::vector<T> _positions;
};

template <typename T>
class solver_degree_II
{
public:
  solver_degree_II() {};

  void set_initial_state(double t0, T x0, T v0)
  {
    _t0 = t0;
    _x0 = x0;
    _v0 = v0;
  }

  void set_timestep(double dt) { _dt = dt; }

  T solve_euler(double tf, std::function<T(double t, T x, T v)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::tie(_positions, _velocities) = euler_explicit(_t0, tf, N, _x0, _v0, a);
    return _positions[N];
  }

  T solve_euler(double tf, std::function<T(double t, T x)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::function<T(double, T, T)> b = [a](double t, T x, T v)
    { return a(t, x); };
    std::tie(_positions, _velocities) = euler_explicit(_t0, tf, N, _x0, _v0, b);
    return _positions[N];
  }

  T solve_midpoint(double tf, std::function<T(double t, T x, T v)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::tie(_positions, _velocities) = midpoint(_t0, tf, N, _x0, _v0, a);
    return _positions[N];
  }

  T solve_midpoint(double tf, std::function<T(double t, T x)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::function<T(double, T, T)> b = [a](double t, T x, T v)
    { return a(t, x); };
    std::tie(_positions, _velocities) = midpoint(_t0, tf, N, _x0, _v0, b);
    return _positions[N];
  }

  T solve_verlet(double tf, std::function<T(double t, T x)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::tie(_positions, _velocities) = verlet_velocity(_t0, tf, N, _x0, _v0, a);
    return _positions[N];
  }

  T solve_yoshida_4th(double tf, std::function<T(double t, T x)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::tie(_positions, _velocities) = yoshida_4th(_t0, tf, N, _x0, _v0, a);
    return _positions[N];
  }

  T solve_RK4(double tf, std::function<T(double t, T x)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::function<T(double, T, T)> b = [a](double t, T x, T v)
    { return a(t, x); };
    std::tie(_positions, _velocities) = RK4_explicit(_t0, tf, N, _x0, _v0, b);
    return _positions[N];
  }

  T solve_RK4(double tf, std::function<T(double t, T x, T v)> a)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::tie(_positions, _velocities) = RK4_explicit(_t0, tf, N, _x0, _v0, a);
    return _positions[N];
  }

  std::vector<T> get_positions() { return _positions; }
  std::vector<T> get_velocities() { return _velocities; }
  std::vector<double> get_timeline(double tf)
  {
    size_t N = get_number_of_steps(_t0, tf, _dt);
    std::vector<double> t(N + 1);
    for (size_t i = 0; i < N + 1; i++)
    {
      t[i] = _dt * i;
    }
    return t;
  }

  std::vector<double> get_energy(double m, std::function<double(double t, T x)> potential)
  {
    std::vector<double> Em(_positions.size());
    double t = _t0;
    for (size_t i = 0; i < _positions.size(); i++)
    {
      Em[i] = m / 2 * _velocities[i].dot(_velocities[i]) + potential(t, _positions[i]);
      t += _dt;
    }
    return Em;
  }
  std::vector<double> get_energy(double m, std::function<double(T x, T v)> kinetic, std::function<double(double t, T x)> potential)
  {
    std::vector<double> Em(_positions.size());
    double t = _t0;
    for (size_t i = 0; i < _positions.size(); i++)
    {
      Em[i] = m / 2 * kinetic(_positions[i], _velocities[i]) + potential(t, _positions[i]);
      t += _dt;
    }
    return Em;
  }

private:
  double _t0;
  double _dt;
  T _x0;
  T _v0;
  std::vector<T> _positions;
  std::vector<T> _velocities;
};