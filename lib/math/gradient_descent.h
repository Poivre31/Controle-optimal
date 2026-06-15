#pragma once

#include "gradient.h"
#include <Eigen/Core>
#include <fmt/format.h>
#include "momentum_descent.h"
#include "lma.h"

/// @brief Find a parameter vector `p` that minimises the function `loss(p)` using constant learning rate gradient descent
/// @param p0 Initial guess
/// @param loss The function from R^N to R to minimise
/// @param i_max Maximum number of iterations
/// @param learning_rate Learning rate, use small values for stability
/// @param threshold Iterations stop if loss is smaller than this. If set to 0, always performs `i_max` iterations
/// @return The optimised parameter vector
template <size_t N>
coordinates<N> gradient_descent(
    coordinates<N> p0, std::function<double(coordinates<N>)> loss, size_t i_max, double learning_rate, double threshold = 0, size_t update_rate = 0)
{
    coordinates<N> p = p0;
    size_t i = 0;

    while (i < i_max)
    {
        coordinates<N> grad = gradient<N>(p, loss);
        p -= learning_rate * grad;
        i++;
        double current_loss = loss(p);
        if (current_loss < threshold)
            break;
        if (current_loss < threshold)
            break;
        else if (update_rate != 0 && i % update_rate == 0)
        {
            fmt::print("Iteration {}/{}, loss: {:.3g}, parameters: ", i, i_max, current_loss);
            for (size_t k = 0; k < N; k++)
            {
                std::cout << " " << p[k];
            }
            std::cout << "\n";
        }
    }

    fmt::print("Completed gradient descent in {} iterations with final loss {:.3g} and optimised parameters ", i, loss(p));
    for (size_t k = 0; k < N; k++)
    {
        std::cout << " " << p[k];
    }
    std::cout << "\n";
    return p;
}

/// @brief Find a parameter vector `p` that minimises the function `loss(p)` using Barzilai-Borwein gradient descent
/// @param p0 Initial guess
/// @param loss The function from R^N to R to minimise
/// @param i_max Maximum number of iterations
/// @param learning_rate Learning rate, use small values for stability
/// @param threshold Iterations stop if loss is smaller than this. If set to 0, always performs `i_max` iterations
/// @return The optimised parameter vector
template <size_t N>
coordinates<N> gradient_descent_BB(
    coordinates<N> p0, std::function<double(coordinates<N>)> loss, size_t i_max, double threshold = 0, size_t update_rate = 0)
{
    coordinates<N> p = p0;
    size_t i = 0;

    coordinates<N> grad = gradient<N>(p, loss);
    coordinates<N> X_n = p - grad * 1e-6;
    while (i < i_max)
    {
        coordinates<N> grad_n = gradient<N>(X_n, loss);
        coordinates<N> ΔX = X_n - p;
        coordinates<N> Δgrad = grad_n - grad;
        double η = 1e-6;
        if (Δgrad.dot(Δgrad) != 0)
            η = abs(ΔX.dot(Δgrad) / Δgrad.dot(Δgrad));
        grad = grad_n;
        p = X_n;
        X_n -= grad * η;
        i++;
        double current_loss = loss(X_n);
        if (current_loss < threshold)
            break;
        else if (update_rate != 0 && i % update_rate == 0)
        {
            fmt::print("Iteration {}/{}, loss: {:.3g}, parameters: ", i, i_max, current_loss);
            for (size_t k = 0; k < N; k++)
            {
                std::cout << " " << X_n[k];
            }
            std::cout << "\n";
        }
    }

    fmt::print("Completed gradient descent in {} iterations with final loss {:.3g} and optimised parameters ", i, loss(X_n));
    for (size_t k = 0; k < N; k++)
    {
        std::cout << " " << X_n[k];
    }
    std::cout << "\n";
    return X_n;
}