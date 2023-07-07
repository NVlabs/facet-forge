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

#include <vector.h>
#include <math_functions.h>
#include <random.h>
#include <bsdf.h>
#include <bsdfs/NDF.h>

//////////////////////////////////////////////////////////////////////////////////
// Shape Invariant NDF
//////////////////////////////////////////////////////////////////////////////////

class ShapeInvariantNDF : public NDF
{
public:
    ShapeInvariantNDF(const BSDF *bsdf, const double roughness_x, const double roughness_y)
        : NDF(bsdf), m_roughness_x(roughness_x), m_roughness_y(roughness_y)
    {
    }

public:
    // roughness
    double m_roughness_x, m_roughness_y;
    // projected roughness in wi
    double roughness_i(const Vector3 &wi) const;

public:
    // distribution of normals (NDF)
    virtual double D(const Vector3 &wm) const;
    // sample the VNDF
    virtual Vector3 sampleD_wi(const Vector3 &wi) const;

public:
    // distribution of slopes
    virtual double P22(const double slope_x, const double slope_y) const = 0;
    // cross section
    virtual double sigma(const Vector3 &wi) const = 0;
    // sample the distribution of visible slopes with roughness=1.0
    virtual Vector2 sampleP22_11(const double theta_i) const = 0;

    // sample a free-path length along direction wr from starting height hr
    // if a collision occurs before escape, return the normal (out_wm) and BSDF (out_bsdf) of the sampled facet
    virtual double sampleHeight(const Vector3 &wr, const double hr, const bool outside,
                                Vector3 &out_wm, const BSDF *&out_bsdf) const;
};

double ShapeInvariantNDF::sampleHeight(const Vector3 &wr, const double hr, const bool outside,
                                       Vector3 &out_wm, const BSDF *&out_bsdf) const
{
    const double sigma_t = sigma(-wr);

    if (sigma_t < 0.00001)
        return (wr.z < 0.0) ? hr : 0.0;

    const double dh = -log(RandomReal()) * wr.z / sigma_t;

    const double h = std::min(0.0, hr) + dh;

    if (h < 0.0)
    {
        out_wm = sampleD_wi(-wr);
        out_bsdf = m_bsdf; // uniform microsurface
    }

    return h;
}

double ShapeInvariantNDF::D(const Vector3 &wm) const
{
    if (wm.z <= 0.0)
        return 0.0;

    // slope of wm
    const double slope_x = -wm.x / wm.z;
    const double slope_y = -wm.y / wm.z;

    return P22(slope_x, slope_y) / (wm.z * wm.z * wm.z * wm.z);
}

Vector3 ShapeInvariantNDF::sampleD_wi(const Vector3 &wi) const
{

    // stretch to match configuration with roughness=1.0
    const Vector3 wi_11 = normalize(Vector3(m_roughness_x * wi.x, m_roughness_y * wi.y, wi.z));

    // sample visible slope with roughness=1.0
    Vector2 slope_11 = sampleP22_11(acos(wi_11.z));

    // align with view direction
    const double phi = atan2(wi_11.y, wi_11.x);
    Vector2 slope(cos(phi) * slope_11.x - sin(phi) * slope_11.y, sin(phi) * slope_11.x + cos(phi) * slope_11.y);

    // stretch back
    slope.x *= m_roughness_x;
    slope.y *= m_roughness_y;

    // if numerical instability
    if ((slope.x != slope.x) || !IsFiniteNumber(slope.x))
    {
        if (wi.z > 0)
            return Vector3(0, 0, 1);
        else
            return normalize(Vector3(wi.x, wi.y, 0));
    }

    // compute normal
    return normalize(Vector3(-slope.x, -slope.y, 1.0));
}

double ShapeInvariantNDF::roughness_i(const Vector3 &wi) const
{
    const double invSinTheta2 = 1.0 / (1.0 - wi.z * wi.z);
    const double cosPhi2 = wi.x * wi.x * invSinTheta2;
    const double sinPhi2 = wi.y * wi.y * invSinTheta2;
    return sqrt(cosPhi2 * m_roughness_x * m_roughness_x + sinPhi2 * m_roughness_y * m_roughness_y);
}