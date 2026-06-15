// #pragma once

// #include "gradient.h"
// #include <Eigen/Core>
// #include <Eigen/SVD>
// #include <fmt/format.h>

// template <size_t N, size_t M>
// coordinates<N> gradient_descent_lma(
//     coordinates<N> X0, std::function<coordinates<M>(coordinates<N>)> f, coordinates<M> target, size_t i_max, double L = 1)
// {
//     coordinates<N> X = X0;
//     matrix<M, N> J;
//     size_t i = 0;
//     double λ = 1e5;
//     double ε4 = .1;

//     while (i < i_max)
//     {
//         double χ2_1 = pow((target - f(X)).norm(), 2);
//         J = jacobian<N, M>(X, f);
//         auto b = J.transpose() * (target - f(X));
//         matrix<N, N> U = J.transpose() * J + λ * matrix<N, N>::Identity();
//         Eigen::JacobiSVD<matrix<N, N>> svd(U);
//         double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
//         while (cond > 100)
//         {
//             U += 1e-4 * U.trace() * matrix<N, N>::Identity();
//             Eigen::JacobiSVD<matrix<N, N>> svd2(U);
//             cond = svd2.singularValues()(0) / svd2.singularValues()(svd2.singularValues().size() - 1);
//         }

//         std::cout << "condtion: " << cond << std::endl;
//         auto δ = U.colPivHouseholderQr().solve(b);
//         double χ2_2 = pow((target - f(X + δ)).norm(), 2);
//         double ρ = (χ2_1 - χ2_2) / abs(δ.transpose() * (λ * δ + b));
//         if (ρ > ε4 || λ >= 1e5)
//         {
//             X += δ;
//             λ = std::max(λ / 9, 1e-5);
//         }
//         else
//         {
//             λ = std::min(11 * λ, 1e5);
//         }
//         i++;
//     }

//     fmt::println("Completed lma in {} iterations with final cost {}", i, 0);
//     return X;
// }