#pragma once

#include "gradient.h"
#include <Eigen/Core>
#include <fmt/format.h>

template <size_t N>
coordinates<N> gradient_descent_ogm1(
    coordinates<N> X0, std::function<double(coordinates<N>)> cost, size_t i_max, double L = 1, double threshold = 1e-5)
{
    coordinates<N> X = X0;
    coordinates<N> Y = X0;
    double θ = 1;
    size_t i = 0;

    while (i < i_max)
    {
        coordinates<N> grad = gradient<N>(X, cost);
        coordinates<N> Y_n = X - 1 / L * grad;
        double θ_n = .5 * (1 + sqrt(1 + 4 * θ * θ));
        X = Y_n + (θ - 1) / θ_n * (Y_n - Y) + θ / θ_n * (Y_n - X);

        Y = Y_n;
        θ = θ_n;
        i++;
        if (cost(X) < threshold)
            break;
    }

    auto grad = gradient<N>(X, cost);
    fmt::println("Completed gradient descent in {} iterations with final cost {}", i, cost(X));
    return X;
}