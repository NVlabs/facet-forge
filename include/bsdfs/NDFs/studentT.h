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

#include <bsdfs/NDFs/beckmann.h>

//////////////////////////////////////////////////////////////////////////////////
// StudentTNDF
//////////////////////////////////////////////////////////////////////////////////

class StudentTNDF : public ShapeInvariantNDF
{
public:
    double m_gamma; // shape parameter
    StudentTNDF(const BSDF *bsdf, const double roughness_x, const double roughness_y, const double gamma)
        : ShapeInvariantNDF(bsdf, roughness_x, roughness_y), m_gamma(gamma){};

    // distribution of slopes
    virtual double P22(const double slope_x, const double slope_y) const;
    // cross section
    virtual double sigma(const Vector3 &wi) const;
    // sample the distribution of visible slopes with roughness=1.0
    virtual Vector2 sampleP22_11(const double theta_i) const;

    virtual Vector3 sampleD_wi(const Vector3 &wi) const;
};

//////////////////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////////////////

// probability mass of the first term in the 3-term m' Beckmann superposition
double p_term1(const double u, const double gamma)
{
    return (3 * Sqrt(3 - 3 * Power(u, 2))) /
           (Power(1 + Power(u, 2) / (3 - 3 * Power(u, 2)), 2.5) *
            (8 * u - (81 * Power(-1 + Power(u, 2), 3)) / Power(3 - 2 * Power(u, 2), 2.5) -
             (Power(u, 2) * Sqrt(3 - 2 * Power(u, 2)) * Sqrt((-1 + Power(u, 2)) / (-3 + 2 * Power(u, 2))) *
              (135 - 210 * Power(u, 2) + 83 * Power(u, 4))) /
                 (Power(-1 + Power(u, 2), 3) * Power((-3 + 2 * Power(u, 2)) / (-1 + Power(u, 2)), 2.5))));
}

// probability mass of the second term in the 3-term m' Beckmann superposition
double p_term2(const double u, const double x)
{
    return Power(u, 0.809494 + 0.170783 / (-1.1224 + x)) /
           (1 + Power(u, 0.145598 + 0.000805627 * x + atan(1.05504 * (-2.91109 + Power(x, 2)))));
}

// Beckmann superposition sampling of Student-T: sample angle-dependent Beckmann roughness
double sample_m_prime(const double u, const double gamma)
{
    if (u < 0.0)
    {
        // approximate method using fitted gen gamma distribution

        const double y = gamma;
        const double a = -1.49293 + y - (0.0655156 * (0.0442664 + u)) / (1.20697 + cos(0.633638 + y));
        const double c = 1.00448 + (0.0138041 + u) / (-0.850096 - u + 0.516877 * pow(y, 2) - asinh(u));

        // {2,2}-Pade approximant of full m1 at u = 0
        const double m1 = (Sqrt(-1 + y) * (-3 + 2 * y) * (-6 * Power(Pi, 1.5) * u * Power(-1 + y, 3.5) * (-4 + 3 * y) * Power(mygamma(-1 + y), 3) + 4 * (4 * Power(-1 + y, 2) * (-3 + Power(u, 2) + (3 + Power(u, 2)) * y) * Power(mygamma(-0.5 + y), 3) + Power(2, 3 - 2 * y) * Pi * u * Power(-1 + y, 2.5) * (-14 + 13 * y) * mygamma(-0.5 + y) * mygamma(-2 + 2 * y) + (Power(Pi, 1.5) * (-6 * (-1 + y) * (-4 + 3 * y) + Power(u, 2) * (8 - 5 * Power(y, 2))) * mygamma(y) * mygamma(-1 + 2 * y)) / Power(4, y)))) /
                          (2. * (16 * (Power(u, 2) * (-2 + y) + 3 * (-1 + y)) * Power(-1 + y, 2.5) *
                                     Power(mygamma(-0.5 + y), 3) +
                                 Power(2, 5 - 2 * y) * Pi * u * Power(-1 + y, 3) * (-20 + 13 * y) *
                                     mygamma(-0.5 + y) * mygamma(-2 + 2 * y) +
                                 Power(Pi, 1.5) * mygamma(y) * (-3 * u * (-3 + 2 * y) * (-4 + 3 * y) * Power(mygamma(y), 2) - Power(2, 3 - 2 * y) * Power(-1 + y, 1.5) * (6 * (-1 + y) * (-4 + 3 * y) + Power(u, 2) * (-2 + y) * (-6 + 5 * y)) * mygamma(-2 + 2 * y))));

        const double b = m1 * mygamma(a) / mygamma(a + 1 / c);

        return pow(RandomGamma(a), 1.0 / c) * b;
    }

    assert(gamma > 2.0);

    double xi1 = RandomReal();
    const double p1 = p_term1(u, gamma);
    const double p2 = p_term2(u, gamma);

    const double gamma_width = 1 / (1 - Power(u, 2) / ((-1 + gamma) * (-1 + Power(u, 2))));

    if (xi1 < p1)
    {
        // term 1
        return RandomGamma(-1.5 + gamma) * gamma_width;
    }
    else
    {
        if (xi1 < p1 + p2)
        {
            // term 2
            return RandomGamma(gamma - 1.0);
        }
        else
        {
            // term 3
            double m = RandomGamma(gamma - 1.0);
            while (RandomReal() > erf(u * Sqrt(-(m / ((-1 + gamma) * (-1 + Power(u, 2)))))))
            {
                m = RandomGamma(gamma - 1.0);
            }
            return m;
        }
    }
}

// auxiliary functions for approximating sigma/Lambda:
double auxF(const double u, const double g)
{
    return atan(2.00141 - 1.6253863790572571 * g) * sin(0.993127 *
                                                        (-1.00658 + u - (0.0209307 * (-2.63062 + g) * u) / (2.19417 + g)) * tan(u));
}

double auxF2(const double x, const double g)
{
    return 1 + 1 / erf(auxF(x / Sqrt(1 + Power(x, 2)), g) / (1 - x / Sqrt(1 + Power(x, 2))));
}

//////////////////////////////////////////////////////////////////////////////////
// implementation
//////////////////////////////////////////////////////////////////////////////////

double StudentTNDF::P22(const double p, const double q) const
{
    return pow((-1 + m_gamma) /
                   (-1 + pow(p, 2) / pow(m_roughness_x, 2) +
                    pow(q, 2) / pow(m_roughness_y, 2) + m_gamma),
               m_gamma) /
           (m_roughness_x * m_roughness_y * Pi);
}

// use an approximation of the cross section to avoid having to use 2F1
double StudentTNDF::sigma(const Vector3 &wi) const
{
    if (wi.z > 0.9999)
        return 1.0;
    if (wi.z < -0.9999)
        return 0.0;

    const double roughnessi = roughness_i(wi);
    const double theta_i = acos(wi.z);
    const double x = 1.0 / tan(theta_i) / roughnessi;

    const double u = wi.z;

    if (u > 0.0)
    {
        return 0.5 * u * auxF2(x, m_gamma);
    }
    else
    {
        return -0.5 * u * auxF2(-x, m_gamma) + u;
    }
}

Vector2 StudentTNDF::sampleP22_11(const double theta_i) const
{
    assert(false); // handled by sampleD_wi()
    return Vector2(0, 0);
}

// vNDF sampling using Beckmann superpositions
Vector3 StudentTNDF::sampleD_wi(const Vector3 &wi) const
{
    const double m_prime = sample_m_prime(wi.z, m_gamma);
    const double beck_rough = 1.0 / sqrt(m_prime / (m_gamma - 1.0));
    BeckmannNDF beck(0, beck_rough * m_roughness_x, beck_rough * m_roughness_y);
    return beck.sampleD_wi(wi);
}