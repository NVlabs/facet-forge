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

class BSDF
{
public:
    // NB: this function takes an input weight and MODIFIES it
    virtual Vector3 sample(const double ior_i, const double ior_t, const Vector3 &wi, double &io_weight) const = 0;

    // sample in a space oriented to normal wm
    Vector3 sample(const double ior_i, const double ior_t, const Vector3 &wi, double &io_weight, const Vector3 &wm) const
    {
        Vector3 w1(0, 0, 0);
        Vector3 w2(0, 0, 0);
        buildOrthonormalBasis(w1, w2, wm);

        Vector3 wi_local(dot(wi, w1), dot(wi, w2), dot(wi, wm));
        Vector3 wo_local = sample(ior_i, ior_t, wi_local, io_weight);

        return wo_local.x * w1 + wo_local.y * w2 + wo_local.z * wm;
    }

    // BRDF * cos(theta_o) evaluation
    virtual double eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const = 0;
    double eval(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo, const Vector3 &wm) const
    {
        Vector3 w1(0, 0, 0);
        Vector3 w2(0, 0, 0);
        buildOrthonormalBasis(w1, w2, wm);

        Vector3 wi_local(dot(wi, w1), dot(wi, w2), dot(wi, wm));
        Vector3 wo_local(dot(wo, w1), dot(wo, w2), dot(wo, wm));

        return eval(ior_i, ior_t, wi_local, wo_local);
    }

    virtual double evalSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo) const = 0;
    // sample in a space oriented to normal wm
    double evalSingular(const double ior_i, const double ior_t, const Vector3 &wi, const Vector3 &wo, const Vector3 &wm) const
    {
        Vector3 w1(0, 0, 0);
        Vector3 w2(0, 0, 0);
        buildOrthonormalBasis(w1, w2, wm);

        Vector3 wi_local(dot(wi, w1), dot(wi, w2), dot(wi, wm));
        Vector3 wo_local(dot(wo, w1), dot(wo, w2), dot(wo, wm));

        return evalSingular(ior_i, ior_t, wi_local, wo_local);
    }
};