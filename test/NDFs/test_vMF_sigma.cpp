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

// test vMF NDF sigma (projected area or cross-section)

#include <bsdfs/NDFs/NullvMF.h>
#include <testing/compare_eval_sample.h>

int main(int argc, char **argv)
{
    srand48(time(NULL));

    if (argc != 3)
    {
        std::cout << "usage: test roughness du\n";
        exit(-1);
    }

    const double roughness = StringToNumber<double>(std::string(argv[1]));
    const double du = StringToNumber<double>(std::string(argv[2]));

    NullvMFNDF ndf(0, roughness);

    for (double u = -1.0; u < 1.0; u += du)
    {
        Vector3 wi(sqrt(1 - u * u), 0.0, u);
        std::cout << ndf.sigma(wi) << " ";
    }
    std::cout << std::endl;

    const size_t NUMSAMPLES = 100000;

    for (double u = -1.0; u < 1.0; u += du)
    {
        double sigma = 0.0;
        for (size_t j = 0; j < NUMSAMPLES; j++)
        {
            Vector2 diskOffset = diskSample2D(0.999999);

            // microfacet normal / sphere position - pre rotation
            Vector3 mPR(diskOffset.x, diskOffset.y, sqrt(1.0 - diskOffset.x * diskOffset.x - diskOffset.y * diskOffset.y));
            // rotate to the same cos(theta)
            Vector3 mPR2(mPR.x * u + sqrt(1.0 - u * u) * mPR.z, mPR.y, u * mPR.z - mPR.x * sqrt(1.0 - u * u));
            // rotate to match azimuth
            const double phi = atan2(sqrt(1 - u * u), 0.0);
            const double cosphi = cos(phi);
            const double sinphi = sin(phi);
            // this is where we strike the unit sphere of the microsurface NDFs - this is then the microfacet normal
            Vector3 microspherePos(-cosphi * mPR2.y + sinphi * mPR2.x, sinphi * mPR2.y + cosphi * mPR2.x, mPR2.z);

            // null scattering unless there is enough density for this microfacet normal - cosine factor is already accounted
            // for in the disk sampling above
            if (RandomReal() < ndf.D(microspherePos))
            {
                sigma += 1.0;
            }
        }
        std::cout << sigma / double(NUMSAMPLES) << " ";
    }
    std::cout << std::endl;

    return 0;
}