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

#include <bsdfs/NDFs/NullNDF.h>

// null-collision von-Mises Fischer (spherical Gaussian/vMF) NDF

class NullvMFNDF : public NullNDF
{
public:
    NullvMFNDF(BSDF *bsdf, double roughness)
        : NullNDF(bsdf, 1), m_roughness(roughness){};

    double m_roughness;

    virtual double D(const Vector3 &wm) const
    {
        // vMF matched to Beckmann roughness, normalized to 1.0 at normal incidence
        const double u = wm.z;
        return exp((2.0 * (-1.0 + u)) / pow(m_roughness, 2));
    };
};
