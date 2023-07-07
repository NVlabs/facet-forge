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
#include <bsdf.h>

//////////////////////////////////////////////////////////////////////////////////
// NDF
//////////////////////////////////////////////////////////////////////////////////

class NDF
{
public:
    NDF(const BSDF *bsdf)
        : m_bsdf(bsdf)
    {
    }

public:
    // BSDF on the microfacets:
    const BSDF *m_bsdf;

public:
    // distribution of normals (NDF)
    virtual double D(const Vector3 &wm) const = 0;
    // distribution of visible normals (vNDF)
    virtual double D_wi(const Vector3 &wi, const Vector3 &wm) const;

    // only used for ShapeInvariant NDF - and included in NullNDF for debugging purposes
    virtual Vector3 sampleD_wi(const Vector3 &wi) const = 0;

public:
    // cross section (projected area) sigma_t when moving in direction wi
    virtual double sigma(const Vector3 &wi) const = 0;

    // sample a free-path length along direction wr from starting height hr
    // if a collision occurs before escape, return the normal (out_wm) and BSDF (out_bsdf) of the sampled facet
    virtual double sampleHeight(const Vector3 &wr, const double hr, const bool outside,
                                Vector3 &out_wm, const BSDF *&out_bsdf) const = 0;

    virtual double evalPhaseFunctionSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo, const bool wi_outside, const bool wo_outside) const;
    // transmittance from a depth h0 in the half space along a ray with direction wi
    virtual double G_1(const Vector3 &wi, const double h0) const;
};

double NDF::G_1(const Vector3 &wi, const double h0) const
{
    if (wi.z <= 0.0)
        return 0.0;

    if (h0 >= 0.0)
        return 1.0;

    const double a = sigma(-wi);
    return exp(h0 / wi.z * a); // exponential transmittance
}

double NDF::D_wi(const Vector3 &wi, const Vector3 &wm) const
{
    // normalization coefficient
    const double l_sigma = sigma(wi);
    if (l_sigma == 0)
        return 0;
    const double c = 1.0 / l_sigma;

    return c * std::max(0.0, dot(wi, wm)) * D(wm);
}

double NDF::evalPhaseFunctionSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo, const bool wi_outside, const bool wo_outside) const
{
    const double etaRatio = ior_t / ior_i;
    double eta = wi_outside ? etaRatio : 1.0 / etaRatio;

    if (wi_outside == wo_outside) // reflection
    {
        // half vector
        const Vector3 wh = normalize(wi + wo);
        // value
        const double value = (wi_outside) ? (0.25 * D_wi(wi, wh) / dot(wi, wh) * m_bsdf->evalSingular(1.0, eta, wi, wo, wh)) : (0.25 * D_wi(-wi, -wh) / dot(-wi, -wh) * m_bsdf->evalSingular(1.0, eta, -wi, -wo, -wh));
        return value;
    }
    else // transmission
    {
        Vector3 wh = -normalize(wi + wo * eta);
        wh *= (wi_outside) ? (sign(wh.z)) : (-sign(wh.z));

        if (dot(wh, wi) < 0)
            return 0;

        float value;
        if (wi_outside)
        {
            value = m_bsdf->evalSingular(1.0, eta, wi, wo, wh) *
                    D_wi(wi, wh) * std::max(0.0, -dot(wo, wh)) *
                    1.0 / pow(dot(wi, wh) + eta * dot(wo, wh), 2.0);
        }
        else
        {
            value = m_bsdf->evalSingular(wi_outside ? eta : 1.0, wi_outside ? 1.0 : eta, -wi, -wo, -wh) *
                    D_wi(-wi, -wh) * std::max(0.0, -dot(-wo, -wh)) *
                    1.0 / pow(dot(-wi, -wh) + eta * dot(-wo, -wh), 2.0);
        }

        return value;
    }
}