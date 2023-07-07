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

#include <cfloat>

const double sign(const double a)
{
    if (a > 0.0)
    {
        return 1.0;
    }
    else
    {
        return -1.0;
    }
}

bool IsFiniteNumber(float x)
{
    return (x <= FLT_MAX && x >= -FLT_MAX);
}

bool IsFiniteNumber(double x)
{
    return (x <= DBL_MAX && x >= -DBL_MAX);
}

double erfinv(double x)
{
    double w, p;
    w = -log((1.0 - x) * (1.0 + x));
    if (w < 5.0)
    {
        w = w - 2.500000;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p * w;
        p = -3.5233877e-06 + p * w;
        p = -4.39150654e-06 + p * w;
        p = 0.00021858087 + p * w;
        p = -0.00125372503 + p * w;
        p = -0.00417768164 + p * w;
        p = 0.246640727 + p * w;
        p = 1.50140941 + p * w;
    }
    else
    {
        w = sqrt(w) - 3.000000;
        p = -0.000200214257;
        p = 0.000100950558 + p * w;
        p = 0.00134934322 + p * w;
        p = -0.00367342844 + p * w;
        p = 0.00573950773 + p * w;
        p = -0.0076224613 + p * w;
        p = 0.00943887047 + p * w;
        p = 1.00167406 + p * w;
        p = 2.83297682 + p * w;
    }
    return p * x;
}

static double abgam(double x)
{
    double gam[10],
        temp;

    gam[0] = 1. / 12.;
    gam[1] = 1. / 30.;
    gam[2] = 53. / 210.;
    gam[3] = 195. / 371.;
    gam[4] = 22999. / 22737.;
    gam[5] = 29944523. / 19733142.;
    gam[6] = 109535241009. / 48264275462.;
    temp = 0.5 * log(2 * M_PI) - x + (x - 0.5) * log(x) + gam[0] / (x + gam[1] / (x + gam[2] / (x + gam[3] / (x + gam[4] / (x + gam[5] / (x + gam[6] / x))))));

    return temp;
}

static double mygamma(double x)
{
    double result;
    result = exp(abgam(x + 5)) / (x * (x + 1) * (x + 2) * (x + 3) * (x + 4));
    return result;
}

static double beta(double m, double n)
{
    return (mygamma(m) * mygamma(n) / mygamma(m + n));
}