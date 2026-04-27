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

#include <bsdfs/dielectric.h>
#include <bsdfs/microsurface.h>
#include <bsdfs/NDFs/beckmann.h>
#include <testing/compare_eval_sample.h>

int main(int argc, char **argv)
{
    srand48(time(NULL));

    if (argc != 9)
    {
        std::cout << "usage: testdiel inner_roughx inner_roughy outer_roughx outer_roughy theta_i ior numsamplesSample numsamplesEval \n";
        exit(-1);
    }

    const float inner_rough_x = StringToNumber<float>(std::string(argv[1]));
    const float inner_rough_y = StringToNumber<float>(std::string(argv[2]));
    const float outer_rough_x = StringToNumber<float>(std::string(argv[3]));
    const float outer_rough_y = StringToNumber<float>(std::string(argv[4]));
    const float theta_i = StringToNumber<float>(std::string(argv[5]));
    const float ior = StringToNumber<float>(std::string(argv[6]));
    size_t numsamplesSample = StringToNumber<size_t>(std::string(argv[7]));
    size_t numsamplesEval = StringToNumber<size_t>(std::string(argv[8]));

    DielectricBSDF micro_brdf;
    BeckmannNDF ndf(&micro_brdf, inner_rough_x, inner_rough_y);
    Microsurface brdf(&ndf);
    BeckmannNDF ndf2(&brdf, outer_rough_x, outer_rough_y);
    Microsurface brdf2(&ndf2);

    compareEvalSample(brdf2, theta_i, numsamplesSample, numsamplesEval, 1.0, ior);

    return 1;
}
