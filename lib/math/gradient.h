#pragma once

#include "../solver/coordinates.h"
#include <functional>

template <size_t N, size_t M>
using matrix = Eigen::Matrix<double, N, M>;

template <size_t N>
coordinates<N> gradient(coordinates<N> X0, std::function<double(coordinates<N>)> f)
{
    double dX = 1e-8;
    coordinates<N> grad;

    for (size_t i = 0; i < N; i++)
    {
        auto X = X0;
        X[i] -= dX;
        auto Ym = f(X);
        X[i] = X0[i] + dX;
        auto Yp = f(X);
        grad[i] = (Yp - Ym) / (2 * dX);
    }

    return grad;
}

template <size_t N, size_t M>
matrix<M, N> jacobian(coordinates<N> X0, std::function<coordinates<M>(coordinates<N>)> f)
{
    double dX = 1e-8;
    matrix<M, N> J;

    for (size_t i = 0; i < N; i++)
    {
        auto X = X0;
        X[i] -= dX;
        auto Ym = f(X);
        X[i] = X0[i] + dX;
        auto Yp = f(X);
        J.col(i) = (Yp - Ym) / (2 * dX);
    }

    return J;
}