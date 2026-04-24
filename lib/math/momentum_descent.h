#pragma once

#include "gradient.h"
#include <Eigen/Core>
#include <fmt/format.h>

template <size_t N>
coordinates<N> gradient_descent_heavy_ball(
    coordinates<N> X0, std::function<double(coordinates<N>)> cost, size_t i_max, double learning_rate = 1, double momentum_parameter = 0.1, double threshold = 1e-5)
{
    coordinates<N> X = X0;
    coordinates<N> M = gradient<N>(X0, cost);
    double θ = 1;
    size_t i = 0;

    while (i < i_max)
    {
        M = momentum_parameter * M + (1 - momentum_parameter) * gradient<N>(X, cost);
        X = X - learning_rate * M;
        i++;
        if (cost(X) < threshold)
            break;
    }

    auto grad = gradient<N>(X, cost);
    fmt::println("Completed gradient descent in {} iterations with final cost {}", i, cost(X));
    return X;
}