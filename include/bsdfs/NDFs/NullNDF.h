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
// General NDF implementation using null scattering
//////////////////////////////////////////////////////////////////////////////////

class NullNDF : public NDF
{
public:
    NullNDF(const BSDF *bsdf, const double majorant)
        : NDF(bsdf), m_majorant(majorant){};
    double m_majorant;

public:
    // distribution of normals (NDF)
    virtual double D(const Vector3 &wm) const = 0;
    virtual double D_wi(const Vector3 &wi, const Vector3 &wm) const;
    // sample the VNDF - for debugging purposes
    virtual Vector3 sampleD_wi(const Vector3 &wi) const;

public:
    // cross section
    virtual double sigma(const Vector3 &wi) const;

    // sample a free-path length along direction wr from starting height hr
    // if a collision occurs before escape, return the normal (out_wm) and BSDF (out_bsdf) of the sampled facet
    virtual double sampleHeight(const Vector3 &wr, const double hr, const bool outside,
                                Vector3 &out_wm, const BSDF *&out_bsdf) const;
};

// projected Area (or sigma_t) for singular NDF with all microfacets having cosine un
double diracSigma(const double u, const double un)
{
    double result = 0.0;
    if (u * un > 0.0)
    {
        result += 2 * Pi * u * un;
    }
    if (std::abs(u) < sqrt(1.0 - un * un))
    {
        result += 2.0 * sqrt(-un * un - u * u + 1);
        const double acosarg = -u * un / (sqrt((1 - u * u) * (1 - un * un)));
        result += 2 * u * un * acos(acosarg);
        if (u * un > 0.0)
        {
            result -= 2 * Pi * u * un;
        }
    }
    return result;
}

double NullNDF::sigma(const Vector3 &wi) const
{
    // quadrature integration over D(wm) using the Dirac NDF sigma() Green's function:
    double result = 0.0;
    for (int i = 0; i < 100; ++i)
    {
        const double u = Gauss100xs[i] * 2 - 1; // rescale abscissas into [-1,1] from [0,1]
        Vector3 wm(sqrt(1.0 - u * u), 0, u);
        result += Gauss100ws[i] * D(wm) * diracSigma(wi.z, u);
    }

    return result * 2 / Pi / m_majorant; // adjust for 2x change of interval length
}

double NullNDF::D_wi(const Vector3 &wi, const Vector3 &wm) const
{

    // normalization coefficient
    const double l_sigma = sigma(wi);
    if (l_sigma == 0)
        return 0;
    const double c = 1.0 / l_sigma;
    return c * std::max(0.0, dot(wi, wm)) * D(wm) / Pi / m_majorant;
}

double NullNDF::sampleHeight(const Vector3 &wr, const double hr, const bool outside,
                             Vector3 &out_wm, const BSDF *&out_bsdf) const
{
    double dh = -log(RandomReal()) * wr.z;
    double h = std::min(0.0, hr) + dh;
    out_bsdf = 0;

    while (h <= 0.0)
    {
        const double u = -wr.z;
        Vector2 diskOffset = diskSample2D(0.999999);

        // microfacet normal / sphere position - pre rotation
        Vector3 mPR(diskOffset.x, diskOffset.y, sqrt(1.0 - diskOffset.x * diskOffset.x - diskOffset.y * diskOffset.y));
        // rotate to the same cos(theta)
        Vector3 mPR2(mPR.x * u + sqrt(1.0 - u * u) * mPR.z, mPR.y, u * mPR.z - mPR.x * sqrt(1.0 - u * u));
        // rotate to match azimuth
        const double phi = atan2(-wr.x, -wr.y);
        const double cosphi = cos(phi);
        const double sinphi = sin(phi);
        // this is where we strike the unit sphere of the microsurface NDFs - this is then the microfacet normal
        Vector3 microspherePos(-cosphi * mPR2.y + sinphi * mPR2.x, sinphi * mPR2.y + cosphi * mPR2.x, mPR2.z);

        // null scattering unless there is enough density for this microfacet normal - cosine factor is already accounted
        // for in the disk sampling above
        if (RandomReal() < D(microspherePos) / m_majorant)
        {
            out_wm = microspherePos;
            assert(dot(wr, out_wm) <= 0.0);
            out_bsdf = m_bsdf;
            return h;
        }
        h += -log(RandomReal()) * wr.z;
    }

    return h;
}

Vector3 NullNDF::sampleD_wi(const Vector3 &wi) const
{
    while (true)
    {
        const double u = wi.z;
        Vector2 diskOffset = diskSample2D(0.999999);

        // microfacet normal / sphere position - pre rotation
        Vector3 mPR(diskOffset.x, diskOffset.y, sqrt(1.0 - diskOffset.x * diskOffset.x - diskOffset.y * diskOffset.y));
        // rotate to the same cos(theta)
        Vector3 mPR2(mPR.x * u + sqrt(1.0 - u * u) * mPR.z, mPR.y, u * mPR.z - mPR.x * sqrt(1.0 - u * u));
        // rotate to match azimuth
        const double phi = atan2(wi.x, wi.y);
        const double cosphi = cos(phi);
        const double sinphi = sin(phi);
        // this is where we strike the unit sphere of the microsurface NDFs - this is then the microfacet normal
        Vector3 microspherePos(-cosphi * mPR2.y + sinphi * mPR2.x, sinphi * mPR2.y + cosphi * mPR2.x, mPR2.z);

        // null scattering unless there is enough density for this microfacet normal - cosine factor is already accounted
        // for in the disk sampling above
        if (RandomReal() < D(microspherePos) / m_majorant)
        {
            return microspherePos;
        }
    }
}