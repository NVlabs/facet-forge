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
#include <random.h>

// smooth mirror
class BlendBSDF : public BSDF
{
public:
    const BSDF *m_bsdfA;
    const BSDF *m_bsdfB;
    const double m_mix_A;

    BlendBSDF(const BSDF *bsdfA, const BSDF *bsdfB, double mix_A) : m_bsdfA(bsdfA), m_bsdfB(bsdfB), m_mix_A(mix_A){};

    virtual Vector3 sample(const double ior_i, const double ior_t, const Vector3 &wi, double &weight) const
    {
        if (RandomReal() < m_mix_A)
        {
            return m_bsdfA->sample(ior_i, ior_t, wi, weight);
        }
        else
        {
            return m_bsdfB->sample(ior_i, ior_t, wi, weight);
        }
    }

    virtual double eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        return m_mix_A * m_bsdfA->eval(ior_i, ior_t, wi, wo) + (1.0 - m_mix_A) * m_bsdfB->eval(ior_i, ior_t, wi, wo);
    }

    virtual double evalSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        return m_mix_A * m_bsdfA->evalSingular(ior_i, ior_t, wi, wo) + (1.0 - m_mix_A) * m_bsdfB->evalSingular(ior_i, ior_t, wi, wo);
    }
};