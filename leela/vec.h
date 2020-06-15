//
//  vec.h
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef vec_h
#define vec_h

#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>

//
// ---------- 1D vector struct ----------------
//

template <typename T>
struct vec1
{
    // --------------- Members ----------------
    union
    {
        std::array<T, 1> data {};                           // zero initialise
        struct
        {
            T x;
        };
    };
    
    // ------------ Constructor ---------------
    vec1() = default;
    
    vec1(T x) : data{x} {}
    
    // ------------ Helper methods ------------
    void zero()
    {
        x = 0;
    }
    
    void one()
    {
        x = 1;
    }
    
    // ------------ Overloaded operators ------
    vec1 operator+(const vec1 &v) const
    {
        return vec1<T>(x + v.x);
    }
    
    void operator+=(const vec1 &v)
    {
        x += v.x;
    }
    
    vec1 operator-(const vec1 &v) const
    {
        return vec1<T>(x - v.x);
    }
    
    void operator-=(const vec1 &v)
    {
        x -= v.x;
    }
    
    vec1 operator*(const T &value) const
    {
        return vec1<T>(x * value);
    }
    
    vec1 operator*(const vec1 &v) const
    {
        return vec1<T>(x * v.x);
    }
    
    void operator*=(const T &value)
    {
        x *= value;
    }
    
    void operator*=(const vec1 &v)
    {
        x *= v.x;
    }
    
    vec1 operator/(const T &value) const
    {
        T inv = 1 / value;
        return vec1<T>(x * inv);
    }
    
    vec1 operator/(const vec1 &v) const
    {
        return vec1<T>(x / v.x);
    }
    
    void operator/=(const T &value)
    {
        T inv = 1 / value;
        x *= inv;
    }
    
    void operator/=(const vec1 &v)
    {
        x /= v.x;
    }
    
    bool operator==(const vec1 &v) const
    {
        return (x == v.x);
    }
    
    bool operator!=(const vec1 &v) const
    {
        return !(*this == v);
    }
    
    vec1 operator-() const
    {
        return vec1<T>(-x);
    }
    
    friend std::ostream& operator<<(std::ostream& os, vec1& v) {
        return os << v.x;
    }
};

//
// ---------- 2D vector struct ---------------
//

template <typename T>
struct vec2
{
    union
    {
        std::array<T, 2> data {};                           // zero initialise
        struct
        {
            T x, y;
        };
    };
    
    // ------------ Constructor ---------------
    vec2() = default;
    
    vec2(T x, T y) : data{x, y} {}
    
    // ------------ Helper methods ------------
    void zero()
    {
        x = y = 0;
    }
    
    void one()
    {
        x = y = 1;
    }
    
    // ------------ Overloaded operators ------
    vec2 operator+(const vec2 &v) const
    {
        return vec2<T>(x + v.x, y + v.y);
    }
    
    void operator+=(const vec2 &v)
    {
        x += v.x;
        y += v.y;
    }
    
    vec2 operator-(const vec2 &v) const
    {
        return vec2<T>(x - v.x, y - v.y);
    }
    
    void operator-=(const vec2 &v)
    {
        x -= v.x;
        y -= v.y;
    }
    
    vec2 operator*(const T &value) const
    {
        return vec2<T>(x * value, y * value);
    }
    
    vec2 operator*(const vec2 &v) const
    {
        return vec2<T>(x * v.x, y * v.y);
    }
    
    void operator*=(const T &value)
    {
        x *= value;
        y *= value;
    }
    
    void operator*=(const vec2 &v)
    {
        x *= v.x;
        y *= v.y;
    }
    
    vec2 operator/(const T &value) const
    {
        T inv = 1 / value;
        return vec2<T>(x * inv, y * inv);
    }
    
    vec2 operator/(const vec2 &v) const
    {
        return vec2<T>(x / v.x, y / v.y);
    }
    
    void operator/=(const T &value)
    {
        T inv = 1 / value;
        x *= inv;
        y *= inv;
    }
    
    void operator/=(const vec2 &v)
    {
        x /= v.x;
        y /= v.y;
    }
    
    bool operator==(const vec2 &v) const
    {
        return ((x == v.x) && (y == v.y));
    }
    
    bool operator!=(const vec2 &v) const
    {
        return !(*this == v);
    }
    
    vec2 operator-() const
    {
        return vec2<T>(-x, -y);
    }
    
    friend std::ostream& operator<<(std::ostream& os, vec2& v) {
        return os << v.x << " " << v.y;
    }
};

//
// ---------- 3D vector struct ---------------
//

template <typename T>
struct vec3
{
    union
    {
        std::array<T, 3> data {};                           // zero initialise
        struct
        {
            T x, y, z;
        };
    };
    
    // ------------ Constructor ---------------
    vec3() = default;
    
    vec3(T x, T y, T z) : data{x, y, z} {}
    
    // ------------ Helper methods ------------
    void zero()
    {
        x = y = z = 0;
    }
    
    void one()
    {
        x = y = z = 1;
    }
    
    // ------------ Overloaded operators ------
    vec3 operator+(const vec3 &v) const
    {
        return vec3<T>(x + v.x, y + v.y, z + v.z);
    }
    
    void operator+=(const vec3 &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
    }
    
    vec3 operator-(const vec3 &v) const
    {
        return vec3<T>(x - v.x, y - v.y, z - v.z);
    }
    
    void operator-=(const vec3 &v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }
    
    vec3 operator*(const T &value) const
    {
        return vec3<T>(x * value, y * value, z * value);
    }
    
    vec3 operator*(const vec3 &v) const
    {
        return vec3<T>(x * v.x, y * v.y, z * v.z);
    }
    
    void operator*=(const T &value)
    {
        x *= value;
        y *= value;
        z *= value;
    }
    
    void operator*=(const vec3 &v)
    {
        x *= v.x;
        y *= v.y;
        z *= v.z;
    }
    
    vec3 operator/(const T &value) const
    {
        T inv = 1 / value;
        return vec3<T>(x * inv, y * inv, z * inv);
    }
    
    vec3 operator/(const vec3 &v) const
    {
        return vec3<T>(x / v.x, y / v.y, z / v.z);
    }
    
    void operator/=(const T &value)
    {
        T inv = 1 / value;
        x *= inv;
        y *= inv;
        z *= inv;
    }
    
    void operator/=(const vec3 &v)
    {
        x /= v.x;
        y /= v.y;
        z /= v.z;
    }
    
    bool operator==(const vec3 &v) const
    {
        return ((x == v.x) && (y == v.y) && (z == v.z));
    }
    
    bool operator!=(const vec3 &v) const
    {
        return !(*this == v);
    }
    
    vec3 operator-() const
    {
        return vec3<T>(-x, -y, -z);
    }
    
    friend std::ostream& operator<<(std::ostream& os, vec3& v) {
        return os << v.x << " " << v.y << " " << v.z;
    }
};

//
// ---------- vector functions ------------
//

// -------- perform lhs * operation -------
// scalar multiplication is commutative: s M = M s
template <typename T>
inline vec1<T> operator*(T const &value, const vec1<T> &rhs)
{
    return rhs * value;
}
template <typename T>
inline vec2<T> operator*(T const &value, const vec2<T> &rhs)
{
    return rhs * value;
}
template <typename T>
inline vec3<T> operator*(T const &value, const vec3<T> &rhs)
{
    return rhs * value;
}

// ------------ dot product ---------------
template <typename T>
inline T dot(const vec1<T> &v1, const vec1<T> &v2)
{
    return (v1.x * v2.x);
}
template <typename T>
inline T dot(const vec2<T> &v1, const vec2<T> &v2)
{
    return (v1.x * v2.x + v1.y * v2.y);
}
template <typename T>
inline T dot(const vec3<T> &v1, const vec3<T> &v2)
{
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// ------------ cross product -------------
template <typename T>
inline vec3<T> cross(const vec3<T> &v1, const vec3<T> &v2)
{
    return vec3<T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// ------------ distance ------------------
template <typename T>
inline T dist(const vec1<T> &v1, const vec1<T> &v2)
{
    T dx = v2.x - v1.x;
    
    return std::sqrt(dx * dx);
}
template <typename T>
inline T dist(const vec2<T> &v1, const vec2<T> &v2)
{
    T dx = v2.x - v1.x;
    T dy = v2.y - v1.y;
    
    return std::sqrt((dx * dx) + (dy * dy));
}
template <typename T>
inline T dist(const vec3<T> &v1, const vec3<T> &v2)
{
    T dx = v2.x - v1.x;
    T dy = v2.y - v1.y;
    T dz = v2.z - v1.z;
    
    return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

// ------------ distance square ------------
template <typename T>
inline T dist2(const vec1<T> &v1, const vec1<T> &v2)
{
    T dx = v2.x - v1.x;
    
    return (dx * dx);
}
template <typename T>
inline T dist2(const vec2<T> &v1, const vec2<T> &v2)
{
    T dx = v2.x - v1.x;
    T dy = v2.y - v1.y;
    
    return ((dx * dx) + (dy * dy));
}
template <typename T>
inline T dist2(const vec3<T> &v1, const vec3<T> &v2)
{
    T dx = v2.x - v1.x;
    T dy = v2.y - v1.y;
    T dz = v2.z - v1.z;
    
    return ((dx * dx) + (dy * dy) + (dz * dz));
}

// ------------ magnitude -------------------
template <typename T>
inline T mag(const vec1<T> &v)
{
    return std::sqrt(v.x * v.x);
}
template <typename T>
inline T mag(const vec2<T> &v)
{
    return std::sqrt((v.x * v.x) + (v.y * v.y));
}
template <typename T>
inline T mag(const vec3<T> &v)
{
    return std::sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

// ------------ magnitude square -----------
template <typename T>
inline T mag2(const vec1<T> &v)
{
    return (v.x * v.x);
}
template <typename T>
inline T mag2(const vec2<T> &v)
{
    return ((v.x * v.x) + (v.y * v.y));
}
template <typename T>
inline T mag2(const vec3<T> &v)
{
    return ((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

// ------------ normalise -------------------
template <typename T>
inline vec1<T> normalise(const vec1<T> &v)
{
    T magnitude = mag(v);
    
    if (magnitude == 0)
    {
        std::cerr << "ERROR: Vector normalise magnitude is zero!" << std::endl;
        exit(1);
    }
    return v / magnitude;
}
template <typename T>
inline vec2<T> normalise(const vec2<T> &v)
{
    T magnitude = mag(v);
    
    if (magnitude == 0)
    {
        std::cerr << "ERROR: Vector normalise magnitude is zero!" << std::endl;
        exit(1);
    }
    return v / magnitude;
}
template <typename T>
inline vec3<T> normalise(const vec3<T> &v)
{
    T magnitude = mag(v);
    
    if (magnitude == 0)
    {
        std::cerr << "ERROR: Vector normalise magnitude is zero!" << std::endl;
        exit(1);
    }
    return v / magnitude;
}

// ------------ minimum -------------------
template <typename T>
inline vec1<T> min(const vec1<T> &v1, const vec1<T> &v2)
{
    return vec1<T>(std::min(v1.x, v2.x));
}
template <typename T>
inline vec2<T> min(const vec2<T> &v1, const vec2<T> &v2)
{
    return vec2<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y));
}
template <typename T>
inline vec3<T> min(const vec3<T> &v1, const vec3<T> &v2)
{
    return vec3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
}

// ------------ maximum -------------------
template <typename T>
inline vec1<T> max(const vec1<T> &v1, const vec1<T> &v2)
{
    return vec1<T>(std::max(v1.x, v2.x));
}
template <typename T>
inline vec2<T> max(const vec2<T> &v1, const vec2<T> &v2)
{
    return vec2<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y));
}
template <typename T>
inline vec3<T> max(const vec3<T> &v1, const vec3<T> &v2)
{
    return vec3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
}

// ------------ ceiling -------------------
template <typename T>
inline vec1<T> ceil(const vec1<T> &v)
{
    return vec1<T>(std::ceil(v.x));
}
template <typename T>
inline vec2<T> ceil(const vec2<T> &v)
{
    return vec2<T>(std::ceil(v.x), std::ceil(v.y));
}
template <typename T>
inline vec3<T> ceil(const vec3<T> &v)
{
    return vec3<T>(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z));
}

// ------------ floor ---------------------
template <typename T>
inline vec1<T> floor(const vec1<T> &v)
{
    return vec1<T>(std::floor(v.x));
}
template <typename T>
inline vec2<T> floor(const vec2<T> &v)
{
    return vec2<T>(std::floor(v.x), std::floor(v.y));
}
template <typename T>
inline vec3<T> floor(const vec3<T> &v)
{
    return vec3<T>(std::floor(v.x), std::floor(v.y), std::floor(v.z));
}

// ------------ round ---------------------
template <typename T>
inline vec1<T> round(const vec1<T> &v)
{
    return vec1<T>(std::round(v.x));
}
template <typename T>
inline vec2<T> round(const vec2<T> &v)
{
    return vec2<T>(std::round(v.x), std::round(v.y));
}
template <typename T>
inline vec3<T> round(const vec3<T> &v)
{
    return vec3<T>(std::round(v.x), std::round(v.y), std::round(v.z));
}

// ------------ type casting --------------
template <typename T, typename U>
inline vec1<T> cast(const vec1<U> &v)
{
    return vec1<T>(static_cast<T>(v.x));
}
template <typename T, typename U>
inline vec2<T> cast(const vec2<U> &v)
{
    return vec2<T>(static_cast<T>(v.x), static_cast<T>(v.y));
}
template <typename T, typename U>
inline vec3<T> cast(const vec3<U> &v)
{
    return vec3<T>(static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z));
}

// ------------ abs ---------------------
template <typename T>
inline vec1<T> abs(const vec1<T> &v)
{
    return vec1<T>(std::abs(v.x));
}
template <typename T>
inline vec2<T> abs(const vec2<T> &v)
{
    return vec2<T>(std::abs(v.x), std::abs(v.y));
}
template <typename T>
inline vec3<T> abs(const vec3<T> &v)
{
    return vec3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

#endif /* vec_h */
