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
// BeckmannNDF
//////////////////////////////////////////////////////////////////////////////////

class BeckmannNDF : public ShapeInvariantNDF
{
public:
    BeckmannNDF(const BSDF *bsdf, const double roughness_x, const double roughness_y)
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

double BeckmannNDF::P22(const double slope_x, const double slope_y) const
{
    return 1.0 / (M_PI * m_roughness_x * m_roughness_y) * exp(-slope_x * slope_x / (m_roughness_x * m_roughness_x) - slope_y * slope_y / (m_roughness_y * m_roughness_y));
}

double BeckmannNDF::sigma(const Vector3 &wi) const
{
    if (wi.z > 0.9999)
        return 1.0;
    if (wi.z < -0.9999)
        return 0.0;

    const double roughnessi = roughness_i(wi);
    const double theta_i = acos(wi.z);
    const double sin_theta_i = sin(theta_i);
    const double a = 1.0 / tan(theta_i) / roughnessi;

    return 0.5 * (erf(a) + 1.0) * wi.z + INV_2_SQRT_M_PI * roughnessi * sin(theta_i) * exp(-a * a);
}

Vector2 BeckmannNDF::sampleP22_11(const double theta_i) const
{
    Vector2 slope;

    const double U = RandomReal();
    const double U_2 = RandomReal();

    if (theta_i < 0.00001)
    {
        const double r = sqrt(-log(U));
        const double phi = 2 * Pi * U_2;
        slope.x = r * cos(phi);
        slope.y = r * sin(phi);
        return slope;
    }

    // constant
    const double sin_theta_i = sin(theta_i);
    const double cos_theta_i = cos(theta_i);

    // slope associated to theta_i
    const double slope_i = cos_theta_i / sin_theta_i;

    // projected area
    const double a = cos_theta_i / sin_theta_i;
    const double sigma = 0.5 * (erf(a) + 1.0) * cos_theta_i + INV_2_SQRT_M_PI * sin_theta_i * exp(-a * a);

    // VNDF normalization factor
    const double c = 1.0 / sigma;

    // search
    double erf_min = -0.9999;
    double erf_max = std::max(erf_min, erf(slope_i));
    double erf_current = 0.5 * (erf_min + erf_max);

    while (erf_max - erf_min > 0.000001)
    {
        if (!(erf_current >= erf_min && erf_current <= erf_max))
            erf_current = 0.5 * (erf_min + erf_max);

        // evaluate slope
        const double slope = erfinv(erf_current);

        // CDF
        const double CDF = (slope >= slope_i) ? 1.0 : c * (INV_2_SQRT_M_PI * sin_theta_i * exp(-slope * slope) + cos_theta_i * (0.5 + 0.5 * erf(slope)));
        const double diff = CDF - U;

        // test estimate
        if (std::abs(diff) < 0.000001)
            break;

        // update bounds
        if (diff > 0.0)
        {
            if (erf_max == erf_current)
                break;
            erf_max = erf_current;
        }
        else
        {
            if (erf_min == erf_current)
                break;
            erf_min = erf_current;
        }

        // update estimate
        const double derivative = 0.5 * c * cos_theta_i - 0.5 * c * sin_theta_i * slope;
        erf_current -= diff / derivative;
    }

    slope.x = erfinv(std::min(erf_max, std::max(erf_min, erf_current)));
    slope.y = erfinv(2.0 * U_2 - 1.0);

    const double u = cos_theta_i;

    if (u < 0.0)
    {
        const double m = u / sqrt(1.0 - u * u);
        double xx;
        if (RandomReal() < powf(-u, 1.3))
        {
            xx = RandomReal() * RandomReal();
        }
        else
        {
            xx = 1.0 - erf(sqrt(-log(RandomReal())));
        }
        slope.x = erfinv(-1.0 + xx * (1.0 + erf(m)));
    }

    if (u < -0.9)
    {
        const double m = u / sqrt(1.0 - u * u);
        const double x = RandomReal() * pow(RandomReal(), -u);
        slope.x = -(sqrt(2 * Power(m, 2) + log(8) - 2 * log((x - 2 * pow(m, 2) * x) / pow(m, 3)) -
                         log(2 * pow(m, 2) + log(8) - 2 * log((x - 2 * pow(m, 2) * x) / pow(m, 3)))) /
                    sqrt(2));
    }

    return slope;
}