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

// Compare eval() and sample() in a lat-long histogram

#pragma once

#include <util.h>
#include <bsdfs/NDF.h>

const int numtheta = 100;
const int numphi = 101;
const int numOrdinates = numphi * numtheta;

int thetaIndex(const Vector3 &o)
{
    const double theta_o = asin(o.z);
    return int(floor((theta_o + 0.5 * M_PI) * numtheta / M_PI));
}

int phiIndex(const Vector3 &o, const double azimuth)
{
    if (o.x == double(0) && o.y == double(0))
        return 0;
    double phi = atan2(o.x, o.y) - azimuth;
    while (phi > M_PI)
    {
        phi -= 2.0 * M_PI;
    }
    while (phi < -M_PI)
    {
        phi += 2.0 * M_PI;
    }
    return int(floor((phi + M_PI) * numphi / (2.0 * M_PI)));
}

int oIndex(const Vector3 &o, const float azimuth)
{
    const int thetai = std::max(0, std::min(numtheta - 1, thetaIndex(o)));
    const int phii = std::max(0, std::min(numphi - 1, phiIndex(o, azimuth)));
    return numphi * thetai + phii;
}

void writeHistogram(double *histogram, const size_t numsamples)
{
    int index = 0;
    for (int i = 0; i < numtheta; ++i)
    {
        for (int j = 0; j < numphi; ++j)
        {
            std::cout << histogram[index] / double(numsamples) << " ";
            index++;
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
}

// BSDF buffers
double g_bsdfsampled[numOrdinates];

const double evalFactor = 2.0 * M_PI * M_PI / (numtheta * numphi);

////////////////////////////////////////////////////////////////////////////
// compareEvalSample:
//
// numsamplesSample: number of samples to test sample() with
// numsamplesEval: number of samples to test eval() with (run for every pixel in the output)
////////////////////////////////////////////////////////////////////////////

void compareEvalSample(BSDF &bsdf, const double theta_i, const size_t numsamplesSample, const size_t numsamplesEval, const double ior_i, const double ior_t)
{
    for (size_t i = 0; i < numOrdinates; ++i)
    {
        g_bsdfsampled[i] = 0.0;
    }

    const double phi = -M_PI * 0.5;
    const Vector3 wi = Vector3(sin(theta_i) * cos(phi), sin(theta_i) * sin(phi), cos(theta_i));

    for (size_t i = 0; i < numsamplesSample; ++i)
    {
        double w(1.0);
        Vector3 wo = bsdf.sample(ior_i, ior_t, wi, w);
        int oi = oIndex(wo, 0.0);
        g_bsdfsampled[oi] += w;
    }

    std::cout << "BSDF.sample():\n";
    writeHistogram(g_bsdfsampled, numsamplesSample);

    std::cout << "BSDF.eval():\n";
    for (int theta_index = 0; theta_index < numtheta; ++theta_index)
    {
        for (int phi_i = 0; phi_i < numphi; ++phi_i)
        {
            double meanEval(0.f);
            for (size_t i = 0; i < numsamplesEval; ++i)
            {
                const double theta = -0.5 * M_PI + double(theta_index + RandomReal()) * M_PI / double(numtheta);
                const double phi = -M_PI + double(phi_i + RandomReal()) * 2.0 * M_PI / double(numphi);
                Vector3 wo = Vector3(cosf(theta) * sinf(phi), cosf(theta) * cosf(phi), sinf(theta));
                meanEval += bsdf.eval(ior_i, ior_t, wi, wo) * cosf(theta);
            }
            std::cout << evalFactor * meanEval / double(numsamplesEval) << " ";
        }
        std::cout << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////
// testVNDF() - samples VNDF and compares to an explicit evaluation
//
// numsamplesSample: number of samples to test sample() with
// numsamplesEval: number of samples to test eval() with (run for every pixel in the output)
////////////////////////////////////////////////////////////////////////////

void testVNDF(NDF &ndf, const double theta_i, const double phi, const size_t numsamplesSample, const size_t numsamplesEval)
{
    for (size_t i = 0; i < numOrdinates; ++i)
    {
        g_bsdfsampled[i] = 0.0;
    }

    const Vector3 wi = Vector3(sin(theta_i) * cos(phi), sin(theta_i) * sin(phi), cos(theta_i));

    for (size_t i = 0; i < numsamplesSample; ++i)
    {
        double w(1.0);
        Vector3 wo = ndf.sampleD_wi(wi);
        int oi = oIndex(wo, 0.0);
        g_bsdfsampled[oi] += w;
    }

    std::cout << "vNDF.sample():\n";
    writeHistogram(g_bsdfsampled, numsamplesSample);

    std::cout << "vNDF.eval():\n";
    for (int theta_index = 0; theta_index < numtheta; ++theta_index)
    {
        for (int phi_i = 0; phi_i < numphi; ++phi_i)
        {
            double meanEval(0.f);
            for (size_t i = 0; i < numsamplesEval; ++i)
            {
                const double theta = -0.5 * M_PI + double(theta_index + RandomReal()) * M_PI / double(numtheta);
                const double phi = -M_PI + double(phi_i + RandomReal()) * 2.0 * M_PI / double(numphi);
                Vector3 wo = Vector3(cosf(theta) * sinf(phi), cosf(theta) * cosf(phi), sinf(theta));
                meanEval += ndf.D_wi(wi, wo) * cosf(theta);
                // meanEval += bsdf.eval(ior_i, ior_t, wi, wo) * cosf(theta);
            }
            std::cout << evalFactor * meanEval / double(numsamplesEval) << " ";
        }
        std::cout << std::endl;
    }
}