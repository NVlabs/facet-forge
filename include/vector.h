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

#include <util.h>

//////////////////////////////////////////////////////////////////////////////////
// Vector2
//////////////////////////////////////////////////////////////////////////////////

struct Vector2
{
    Vector2(const double in_x, const double in_y)
        : x(in_x), y(in_y){};

    Vector2()
        : x(0.0), y(0.0){};

    double x, y;

    Vector2 operator*(const double c) const
    {
        return Vector2(c * this->x, c * this->y);
    };

    Vector2 operator/(const double c) const
    {
        return Vector2(this->x / c, this->y / c);
    };

    Vector2 operator+(Vector2 b) const
    {
        return Vector2(b.x + this->x, b.y + this->y);
    };

    Vector2 operator-() const
    {
        return Vector2(-this->x, -this->y);
    };

    Vector2 operator+=(Vector2 b)
    {
        this->x += b.x;
        this->y += b.y;
        return *this;
    };

    Vector2 operator-=(Vector2 b)
    {
        this->x -= b.x;
        this->y -= b.y;
        return *this;
    };
};

Vector2 operator*(const double c, const Vector2 &v)
{
    return Vector2(c * v.x, c * v.y);
}

Vector2 operator-(const Vector2 &a, const Vector2 &b)
{
    return Vector2(a.x - b.x, a.y - b.y);
}

bool operator==(const Vector2 &a, const Vector2 &b)
{
    return a.x == b.x && a.y == b.y;
}

bool operator!=(const Vector2 &a, const Vector2 &b)
{
    return a.x != b.x || a.y != b.y;
}

double dot(const Vector2 &a, const Vector2 &b)
{
    return a.x * b.x + a.y * b.y;
}

double Norm(const Vector2 &v)
{
    return sqrt(v.x * v.x + v.y * v.y);
}

Vector2 normalize(const Vector2 &v)
{
    return v * (1.0 / Norm(v));
}

//////////////////////////////////////////////////////////////////////////////////
// Vector3
//////////////////////////////////////////////////////////////////////////////////

struct Vector3
{
    Vector3(const double in_x, const double in_y, const double in_z)
        : x(in_x), y(in_y), z(in_z){};

    Vector3()
        : x(0.0), y(0.0), z(0.0){};
    double x, y, z;

    Vector3(const double *ptr)
    {
        x = ptr[0];
        y = ptr[1];
        z = ptr[2];
    }

    Vector3 operator*(const Vector3 &c) const
    {
        return Vector3(c.x * this->x, c.y * this->y, c.z * this->z);
    };

    double min() const
    {
        return std::min(x, std::min(y, z));
    }

    double max() const
    {
        return std::max(x, std::max(y, z));
    }

    double sum() const
    {
        return x + y + z;
    }

    Vector3 operator*(const double c) const
    {
        return Vector3(c * this->x, c * this->y, c * this->z);
    };

    Vector3 operator/(const double c) const
    {
        return Vector3(this->x / c, this->y / c, this->z / c);
    };

    Vector3 operator+(Vector3 b) const
    {
        return Vector3(b.x + this->x, b.y + this->y, b.z + this->z);
    };

    Vector3 operator-() const
    {
        return Vector3(-this->x, -this->y, -this->z);
    };

    Vector3 operator+=(Vector3 b)
    {
        this->x += b.x;
        this->y += b.y;
        this->z += b.z;

        return *this;
    };

    Vector3 operator*=(Vector3 b)
    {
        this->x *= b.x;
        this->y *= b.y;
        this->z *= b.z;

        return *this;
    };

    Vector3 operator*=(double b)
    {
        this->x *= b;
        this->y *= b;
        this->z *= b;

        return *this;
    };

    Vector3 operator-=(Vector3 b)
    {
        this->x -= b.x;
        this->y -= b.y;
        this->z -= b.z;

        return *this;
    };
};

Vector3 operator*(const double c, const Vector3 &v)
{
    return Vector3(c * v.x, c * v.y, c * v.z);
}

Vector3 operator/(const double c, const Vector3 &v)
{
    return Vector3(c / v.x, c / v.y, c / v.z);
}

Vector3 operator-(const Vector3 &a, const Vector3 &b)
{
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector3 operator/(const Vector3 &a, const Vector3 &b)
{
    return Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
}

bool operator==(const Vector3 &a, const Vector3 &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool operator!=(const Vector3 &a, const Vector3 &b)
{
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

double dot(const Vector3 &a, const Vector3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double Norm(const Vector3 &v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector3 normalize(const Vector3 &v)
{
    return v * (1.0 / Norm(v));
}

Vector3 Clamp(const Vector3 &v, const double min, const double max)
{
    return Vector3(
        Clamp(v.x, min, max),
        Clamp(v.y, min, max),
        Clamp(v.z, min, max));
}

Vector3 Cross(const Vector3 &a, const Vector3 &b)
{
    return Vector3(-b.y * a.z + a.y * b.z, b.x * a.z - a.x * b.z, -b.x * a.y + a.x * b.y);
}

Vector3 reflect(const Vector3 &in, const Vector3 &n)
{
    return -in + 2.0 * dot(in, n) * n;
}

// build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
void buildOrthonormalBasis(Vector3 &omega_1, Vector3 &omega_2, const Vector3 &omega_3)
{
    if (omega_3.z < -0.9999999f)
    {
        omega_1 = Vector3(0.0f, -1.0f, 0.0f);
        omega_2 = Vector3(-1.0f, 0.0f, 0.0f);
    }
    else
    {
        const float a = 1.0f / (1.0f + omega_3.z);
        const float b = -omega_3.x * omega_3.y * a;
        omega_1 = Vector3(1.0f - omega_3.x * omega_3.x * a, b, -omega_3.x);
        omega_2 = Vector3(b, 1.0f - omega_3.y * omega_3.y * a, -omega_3.y);
    }
}
