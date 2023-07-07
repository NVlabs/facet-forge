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

// test vNDF sampling for student-T NDF

#include <bsdfs/NDFs/studentT.h>
#include <testing/compare_eval_sample.h>

int main(int argc, char **argv)
{
    srand48(time(NULL));

    if (argc != 8)
    {
        std::cout << "usage: test roughx roughy gamma theta_i phi numsamplesEval numsamplesSample \n";
        exit(-1);
    }

    const double rough_x = StringToNumber<double>(std::string(argv[1]));
    const double rough_y = StringToNumber<double>(std::string(argv[2]));
    const double gamma = StringToNumber<double>(std::string(argv[3]));
    const double theta_i = StringToNumber<double>(std::string(argv[4]));
    const double phi = StringToNumber<double>(std::string(argv[5]));
    const size_t numsamplesEval = StringToNumber<size_t>(std::string(argv[6]));
    const size_t numsamplesSample = StringToNumber<size_t>(std::string(argv[7]));

    StudentTNDF ndf(0, rough_x, rough_y, gamma);

    testVNDF(ndf, theta_i, phi, numsamplesEval, numsamplesSample);

    return 0;
}