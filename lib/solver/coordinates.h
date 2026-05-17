#pragma once

#include <Eigen/Dense>
#include <fmt/format.h>

template <size_t N>
using coordinates = Eigen::Vector<double, N>;

using xyz = coordinates<3>;
using vec3 = coordinates<3>;
using xy = coordinates<2>;
using vec2 = coordinates<2>;
using polar = coordinates<2>;
using cyl = coordinates<3>;
using sph = coordinates<3>;

double sign(double x)
{
    if (x > 0)
        return 1;
    else if (x == 0)
        return 0;
    else
        return -1;
}

vec3 dim_2_to_dim_3(vec2 v)
{
    return {v[0], v[1], 0.};
}
vec2 dim_3_to_dim_2(vec3 v)
{
    return {v[0], v[1]};
}
vec2 cross_2D(vec2 v, vec3 outside)
{
    return dim_3_to_dim_2(dim_2_to_dim_3(v).cross(outside));
}

template <typename T>
std::vector<double> parse(std::vector<T> data, size_t index)
{
    std::vector<double> result(data.size());
    for (size_t i = 0; i < data.size(); i++)
    {
        result[i] = data[i][index];
    }
    return result;
}

template <typename T>
std::vector<T> stride(std::vector<T> data, size_t stride)
{
    if (stride < 1)
        fmt::println("ERROR: stride must be positive");
    std::vector<T> result;
    for (size_t i = 0; i < data.size(); i += stride)
    {
        result.push_back(data[i]);
    }
    return result;
}

template <typename T, typename U>
std::vector<U> apply_element_wise(std::vector<T> X, std::function<U(T)> f)
{
    std::vector<U> Y(X.size());
    for (size_t i = 0; i < X.size(); i++)
    {
        Y[i] = f(X[i]);
    }
    return Y;
}

template <typename T, typename U>
std::vector<U> apply_element_wise(std::vector<T> X, std::vector<T> Y, std::function<U(T, T)> f)
{
    std::vector<U> Z(X.size());
    for (size_t i = 0; i < X.size(); i++)
    {
        Z[i] = f(X[i], Y[i]);
    }
    return Z;
}

vec3 get_radial_vector(xyz r)
{
    return r.normalized();
}

vec2 get_radial_vector(xy r)
{
    return r.normalized();
}

std::vector<polar> cart_to_polar(std::vector<xy> XY)
{
    std::vector<polar> pol(XY.size());
    for (size_t i = 0; i < XY.size(); i++)
    {
        pol[i] = {XY[i].norm(), atan2(XY[i][1], XY[i][0])};
    }

    return pol;
}

std::vector<xy> polar_to_cart(std::vector<polar> pol)
{
    std::vector<xy> XY(pol.size());
    for (size_t i = 0; i < pol.size(); i++)
    {
        XY[i] = {cos(pol[i][1]) * pol[i][0], sin(pol[i][1]) * pol[i][0]};
    }

    return XY;
}

size_t get_maximum_index(std::vector<double> data)
{
    size_t index = data.size() - 1;
    double value = 0;
    for (size_t i = 0; i < data.size() - 1; i++)
    {
        if (data[i + 1] < data[i] && data[i] > value)
        {
            value = data[i];
            index = i;
        }
    }
    return index;
}

template <typename T>
std::vector<T> truncate_vector(std::vector<T> data, size_t max_index)
{
    std::vector<T> result(max_index + 1);
    for (size_t i = 0; i < result.size(); i++)
    {
        result[i] = data[i];
    }
    return result;
}