#pragma once

#include "gradient.h"
#include <Eigen/Core>
#include <fmt/format.h>
#include "momentum_descent.h"
#include "lma.h"

template <size_t N, coordinates<N> target>
double L2_loss(coordinates<N> X)
{
    return (X - target).norm();
}

template <size_t N>
coordinates<N> gradient_descent(
    coordinates<N> X0, std::function<double(coordinates<N>)> cost, size_t i_max, double parameter = 0, double threshold = 1e-5)
{
    coordinates<N> X = X0;
    size_t i = 0;
    if (parameter != 0)
    {
        while (i < i_max)
        {
            coordinates<N> grad = gradient<N>(X, cost);
            X -= parameter * grad;
            i++;
            if (cost(X) < threshold)
                break;
        }
    }
    else
    {
        coordinates<N> grad = gradient<N>(X, cost);
        coordinates<N> X_n = X - grad * 1e-6;
        while (i < i_max)
        {
            coordinates<N> grad_n = gradient<N>(X_n, cost);
            coordinates<N> ΔX = X_n - X;
            coordinates<N> Δgrad = grad_n - grad;
            double η = 1e-6;
            if (Δgrad.dot(Δgrad) != 0)
                η = abs(ΔX.dot(Δgrad) / Δgrad.dot(Δgrad));
            grad = grad_n;
            X = X_n;
            X_n -= grad * η;
            i++;
            if (cost(X_n) < threshold)
                break;
        }
        X = X_n;
    }
    auto grad = gradient<N>(X, cost);
    fmt::println("Completed gradient descent in {} iterations with final cost {}", i, cost(X));
    return X;
}