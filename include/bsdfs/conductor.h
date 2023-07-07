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

#include <bsdf.h>
#include <fresnel.h>

// smooth conductor
class ConductorBRDF : public BSDF
{
public:
    double m_eta; // real part of ior
    double m_k;   // imaginary part of ior

    ConductorBRDF(const double eta, const double k) : m_eta(eta), m_k(k){};

    virtual Vector3 sample(const double ior_i, const double ior_t, const Vector3 &wi, double &weight) const
    {
        // NB: ior_t is ignored, since this comes from the conductor properties directly on creation
        const double FR = ConductorR(wi.z, ior_i, m_eta, m_k);
        weight *= FR;
        return reflect(wi, Vector3(0, 0, 1));
    }

    virtual double eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        return 0.0;
    }

    virtual double evalSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        return ConductorR(wi.z, ior_i, m_eta, m_k);
    }
};