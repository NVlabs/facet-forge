/*
 * Copyright (c) <2023> NVIDIA CORPORATION & AFFILIATES. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <bsdfs/NDFs/ShapeInvariantNDF.h>

//////////////////////////////////////////////////////////////////////////////////
// GGXNDF
//////////////////////////////////////////////////////////////////////////////////

class GGXNDF : public ShapeInvariantNDF
{
public:
    GGXNDF(const BSDF *bsdf, const double roughness_x, const double roughness_y)
        : ShapeInvariantNDF(bsdf, roughness_x, roughness_y)
    {
    }

    // distribution of slopes
    virtual double P22(const double slope_x, const double slope_y) const;
    // cross section
    virtual double sigma(const Vector3 &wi) const;
    // sample the distribution of visible slopes with roughness=1.0
    virtual Vector2 sampleP22_11(const double theta_i) const;
};

double GGXNDF::P22(const double slope_x, const double slope_y) const
{
    const double tmp = 1.0 + slope_x * slope_x / (m_roughness_x * m_roughness_x) + slope_y * slope_y / (m_roughness_y * m_roughness_y);
    return 1.0 / (M_PI * m_roughness_x * m_roughness_y) / (tmp * tmp);
}

Vector2 GGXNDF::sampleP22_11(const double theta_i) const
{
    Vector2 slope;

    const double U = RandomReal();
    const double U_2 = RandomReal();

    if (theta_i < 0.0001)
    {
        const double r = sqrt(U / (1 - U));
        const double phi = 2 * Pi * U_2;
        slope.x = r * cos(phi);
        slope.y = r * sin(phi);
        return slope;
    }

    // constant
    const double sin_theta_i = sin(theta_i);
    const double cos_theta_i = cos(theta_i);
    const double tan_theta_i = sin_theta_i / cos_theta_i;

    // slope associated to theta_i
    const double slope_i = cos_theta_i / sin_theta_i;

    // projected area
    const double sigma = 0.5 * (cos_theta_i + 1);
    if (sigma < 0.0001f || sigma != sigma)
        return Vector2(0, 0);

    // normalization coefficient
    const double c = 1.0 / sigma;

    const double A = 2 * U / cos_theta_i / c - 1;
    const double B = tan_theta_i;
    const double tmp = 1 / (A * A - 1);

    const double D = sqrt(std::max(0.0, B * B * tmp * tmp - (A * A - B * B) * tmp));
    const double slope_x_1 = B * tmp - D;
    const double slope_x_2 = B * tmp + D;
    slope.x = (A < 0.0 || slope_x_2 > 1.0 / tan_theta_i) ? slope_x_1 : slope_x_2;
    slope.y = sqrt(-1 - slope.x * slope.x + (1 + slope.x * slope.x) / pow(1 - U_2, 2.0 / 3.0)) * sin(2 * Pi * RandomReal());

    return slope;
}

double GGXNDF::sigma(const Vector3 &wi) const
{
    if (wi.z > 0.9999)
        return 1.0;
    if (wi.z < -0.9999)
        return 0.0;

    const double theta_i = acos(wi.z);
    const double sin_theta_i = sin(theta_i);
    const double roughnessi = roughness_i(wi);

    return 0.5 * (wi.z + sqrt(wi.z * wi.z + sin_theta_i * sin_theta_i * roughnessi * roughnessi));
}