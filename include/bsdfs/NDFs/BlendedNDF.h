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
// BlendedNDF
//////////////////////////////////////////////////////////////////////////////////

class BlendedNDF : public NDF
{
public:
    BlendedNDF(const NDF *ndf1, const NDF *ndf2, const double w1) // mix between NDFs ndf1 and ndf2 where ndf1 has a weight of w1 in the mixture
        : NDF(ndf1->m_bsdf), m_w1(w1), m_ndf1(ndf1), m_ndf2(ndf2){};
    double m_w1;
    const NDF *m_ndf1;
    const NDF *m_ndf2;

public:
    // distribution of normals (NDF)
    virtual double D(const Vector3 &wm) const;
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

double BlendedNDF::D(const Vector3 &wm) const
{
    return m_w1 * m_ndf1->D(wm) + (1.0 - m_w1) * m_ndf2->D(wm);
}

double BlendedNDF::sigma(const Vector3 &wi) const
{
    return m_w1 * m_ndf1->sigma(wi) + (1.0 - m_w1) * m_ndf2->sigma(wi);
}

double BlendedNDF::sampleHeight(const Vector3 &wr, const double hr, const bool outside,
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

Vector3 BlendedNDF::sampleD_wi(const Vector3 &wi) const
{
    double sigma1 = m_ndf1->sigma(wi);
    double sigma2 = m_ndf2->sigma(wi);
    double p1 = m_w1 * sigma1 / (m_w1 * sigma1 + (1.0 - m_w1) * sigma2);
    if (RandomReal() < p1)
    {
        return m_ndf1->sampleD_wi(wi);
    }
    else
    {
        return m_ndf2->sampleD_wi(wi);
    }
}