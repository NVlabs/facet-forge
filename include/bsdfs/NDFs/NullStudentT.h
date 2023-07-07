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

// null-collision student-T NDF

class NullStudentTNDF : public NullNDF
{
public:
    NullStudentTNDF(BSDF *bsdf, double roughness, double gamma, double majorant)
        : NullNDF(bsdf, majorant), m_roughness(roughness), m_gamma(gamma){};

    double m_gamma, m_roughness;

    virtual double D(const Vector3 &wm) const
    {
        const double u = wm.z;
        if (u > 0.0 && u <= 1.0)
        {
            return 1 / (Pi * Power(u, 4) * Power(m_roughness, 2) * Power(1 + (1 - Power(u, 2)) / (Power(u, 2) * Power(m_roughness, 2) * (-1 + m_gamma)), m_gamma));
        }
        elseß
        {
            return 0.0;
        }
    };
};