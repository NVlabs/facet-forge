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
#include <random.h>

// smooth dielectric
class DielectricBSDF : public BSDF
{
public:
    DielectricBSDF(){};

    virtual Vector3 sample(const double ior_i, const double ior_t, const Vector3 &wi, double &weight) const
    {
        const double FR = DielectricR(wi.z, ior_t / ior_i);
        if (RandomReal() <= FR)
        {
            return reflect(wi, Vector3(0, 0, 1));
        }
        else
        {
            return refract(-wi, Vector3(0, 0, 1), ior_i, ior_t);
        }
    }

    virtual double eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        return 0.0;
    }

    virtual double evalSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const
    {
        double eta = ior_t / ior_i;
        if (wo.z >= 0.0)
        {
            return DielectricR(wi.z, eta);
        }
        else
        {
            return (1.0 - DielectricR(wi.z, eta)) * eta * eta;
        }
    }
};