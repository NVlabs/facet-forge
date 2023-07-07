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

#include <bsdfs/NDF.h>
#include <bsdf.h>
#include <random.h>

#define MAX_WALK_LENGTH 10

class Microsurface : public BSDF
{
public:
    const NDF *m_ndf;

    Microsurface(NDF *ndf) : m_ndf(ndf)
    {
    }

    virtual Vector3 sample(const double ior_i, const double ior_t, const Vector3 &wi, double &io_weight) const;
    virtual double eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const;
    virtual double evalSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        return 0.0; // TODO - special case this when both roughnesses are 0
    }
};

Vector3 Microsurface::sample(const double ior_i, const double ior_t, const Vector3 &wi, double &io_weight) const
{
    if (wi.z < 0)
    {
        io_weight = 0;
        return Vector3(0, 0, 1);
    }

    Vector3 wr = -wi;
    double hr = 0.0;
    bool outside = true;

    // random walk
    size_t collision_count = 0;
    while (collision_count < MAX_WALK_LENGTH + 1)
    {
        // next height
        Vector3 wm; // microfacet normal
        const BSDF *microfacet_bsdf = 0;
        hr = m_ndf->sampleHeight(wr, hr, outside, wm, microfacet_bsdf);

        // leave the microsurface?
        if (hr >= 0.0)
            break;

        assert(0 != microfacet_bsdf);

        // next direction
        collision_count++;
        wr = microfacet_bsdf->sample(outside ? ior_i : ior_t, outside ? ior_t : ior_i, -wr, io_weight, wm);
        if (dot(wr, wm) < 0.0)
        {
            outside = !outside;
            wr = -wr;
            hr = log(1.0 - exp(hr));
        }

        // if NaN (should not happen, just in case)
        if ((hr != hr) || (wr.z != wr.z))
        {
            io_weight = 0.0;
            return Vector3(0, 0, 1);
        }
    }

    if (hr < 0.0)
    {
        io_weight = 0.0;
    }

    return outside ? wr : -wr;
}

double Microsurface::eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
{
    if (wi.z < 0)
        return 0;

    Vector3 wr = -wi; // direction of the ray
    double hr = 0.0;  // height of the ray
    bool outside = true;
    double sum = 0;
    double weight = 1.0;

    // random walk
    size_t collision_count = 0;
    while (collision_count < MAX_WALK_LENGTH)
    {
        // next height
        Vector3 wm; // microfacet normal
        const BSDF *microfacet_bsdf = 0;
        hr = m_ndf->sampleHeight(wr, hr, outside, wm, microfacet_bsdf);

        // leave the microsurface?
        if (hr >= 0.0)
            break;

        assert(0 != microfacet_bsdf);

        // next event estimation
        const double phaseFunctionSingular = m_ndf->evalPhaseFunctionSingular(ior_i, ior_t, outside ? -wr : wr, wo, outside, (wo.z > 0));
        const double phaseFunction = microfacet_bsdf->eval(outside ? ior_i : ior_t, outside ? ior_t : ior_i, -wr, wo, wm);
        const double hr_inside = outside ? log(1.0 - exp(hr)) : hr;
        const double hr_outside = outside ? hr : log(1.0 - exp(hr));
        const double shadowingSingular = (wo.z > 0) ? m_ndf->G_1(wo, hr_outside) : m_ndf->G_1(-wo, hr_inside);

        // apply the shadowing that aligns with what side we started on and whether the facet reflected or not
        bool facet_reflected = dot(wo, wm) >= 0.0;
        assert(dot(-wr, wm) >= 0.0); // assuming the facet always faces the ray
        const double shadowing = (facet_reflected == outside) ? m_ndf->G_1(wo, hr_outside) : m_ndf->G_1(-wo, hr_inside);
        const double I = weight * (phaseFunctionSingular * shadowingSingular + phaseFunction * shadowing);

        if (IsFiniteNumber(I))
            sum += I;

        // next direction
        collision_count++;
        wr = microfacet_bsdf->sample(outside ? ior_i : ior_t, outside ? ior_t : ior_i, -wr, weight, wm);
        if (dot(wr, wm) < 0.0)
        {
            outside = !outside;
            wr = -wr;
            hr = log(1.0 - exp(hr));
        }

        // if NaN (should not happen, just in case)
        if ((hr != hr) || (wr.z != wr.z))
        {
            return 0.0;
        }
    }

    return sum;
}
