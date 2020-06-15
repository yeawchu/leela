//
//  ten.h
//  leela
//
//  Created by Yeaw Chu Lee on 23/06/2016.
//  Copyright © 2016 Yeaw Chu Lee. All rights reserved.
//

#ifndef ten_h
#define ten_h

#include <cmath>
#include <utility>
#include "vec.h"

//
// ---------- 2x2 tensor struct ---------------
//

template <typename T>
struct ten22
{
    // --------------- Members ----------------
    union
    {
        std::array<T, 4> data {};                           // zero initialise
        struct
        {
            T xx, xy,
            yx, yy;
        };
    };
    
    // ------------ Constructor ---------------
    ten22() = default;
    
    ten22(T xx, T xy, T yx, T yy) : data{xx, xy, yx, yy} {}
    
    // ------------ Helper methods ------------
    void zero()
    {
        xx = 0; xy = 0;
        yx = 0; yy = 0;
    }
    
    void identity()
    {
        xx = 1; xy = 0;
        yx = 0; yy = 1;
    }
    
    void transpose()
    {
    std:swap(xy, yx);
    }
    
    // ------------ Overloaded operators ------------
    ten22 operator+(const ten22 &t) const
    {
        return ten22<T>(xx + t.xx, xy + t.xy,
                        yx + t.yx, yy + t.yy);
    }
    
    void operator+=(const ten22 &t)
    {
        xx += t.xx; xy += t.xy;
        yx += t.yx; yy += t.yy;
    }
    
    ten22 operator-(const ten22 &t) const
    {
        return ten22<T>(xx - t.xx, xy - t.xy,
                        yx - t.yx, yy - t.yy);
    }
    
    void operator-=(const ten22 &t)
    {
        xx -= t.xx; xy -= t.xy;
        yx -= t.yx; yy -= t.yy;
    }
    
    ten22 operator*(const T &value) const
    {
        return ten22<T>(xx * value, xy * value,
                        yx * value, yy * value);
    }
    
    void operator*=(const T &value)
    {
        xx *= value; xy *= value;
        yx *= value; yy *= value;
    }
    
    ten22 operator/(const T &value) const
    {
        T inv = 1 / value;
        return ten22<T>(xx * inv, xy * inv,
                        yx * inv, yy * inv);
    }
    
    void operator/=(const T &value)
    {
        T inv = 1 / value;
        xx *= inv; xy *= inv;
        yx *= inv; yy *= inv;
    }
    
    bool operator==(const ten22 &t) const
    {
        return ((xx == t.xx) && (xy == t.xy) &&
                (yx == t.yx) && (yy == t.yy));
    }
    
    bool operator!=(const ten22 &t) const
    {
        return !(*this == t);
    }
    
    ten22 operator-() const
    {
        return ten22<T>(-xx, -xy,
                        -yx, -yy);
    }
    
    friend std::ostream& operator<<(std::ostream& os, ten22& t) {
        return os << "{xx:" << t.xx << " xy:" << t.xy << "}" << std::endl
        << "{yx:" << t.yx << " yy:" << t.yy << "}";
    }
};

//
// ---------- 3x3 tensor struct ---------------
//

template <typename T>
struct ten33
{
    // --------------- Members ----------------
    union
    {
        std::array<T, 9> data {};                          // zero initialise
        struct
        {
            T xx, xy, xz,
            yx, yy, yz,
            zx, zy, zz;
        };
    };
    
    // ------------ Constructor ---------------
    ten33() = default;
    
    ten33(T xx, T xy, T xz, T yx, T yy, T yz, T zx, T zy, T zz) : data{xx, xy, xz, yx, yy, yz, zx, zy, zz} {}
    
    // ------------ Helper methods ------------
    void zero()
    {
        xx = 0; xy = 0; xz = 0;
        yx = 0; yy = 0; yz = 0;
        zx = 0; zy = 0; zz = 0;
    }
    
    void identity()
    {
        xx = 1; xy = 0; xz = 0;
        yx = 0; yy = 1; yz = 0;
        zx = 0; zy = 0; zz = 1;
    }
    
    void transpose()
    {
        std::swap(xy, yx);
        std::swap(xz, zx);
        std::swap(zy, yz);
    }
    
    // ------------ Overloaded operators ------------
    ten33 operator+(const ten33 &t) const
    {
        return ten33<T>(xx + t.xx, xy + t.xy, xz + t.xz,
                        yx + t.yx, yy + t.yy, yz + t.yz,
                        zx + t.zx, zy + t.zy, zz + t.zz);
    }
    
    void operator+=(const ten33 &t)
    {
        xx += t.xx; xy += t.xy; xz += t.xz;
        yx += t.yx; yy += t.yy; yz += t.yz;
        zx += t.zx; zy += t.zy; zz += t.zz;
    }
    
    ten33 operator-(const ten33 &t) const
    {
        return ten33<T>(xx - t.xx, xy - t.xy, xz - t.xz,
                        yx - t.yx, yy - t.yy, yz - t.yz,
                        zx - t.zx, zy - t.zy, zz - t.zz);
    }
    
    void operator-=(const ten33 &t)
    {
        xx -= t.xx; xy -= t.xy; xz -= t.xz;
        yx -= t.yx; yy -= t.yy; yz -= t.yz;
        zx -= t.zx; zy -= t.zy; zz -= t.zz;
    }
    
    ten33 operator*(const T &value) const
    {
        return ten33<T>(xx * value, xy * value, xz * value,
                        yx * value, yy * value, yz * value,
                        zx * value, zy * value, zz * value);
    }
    
    void operator*=(const T &value)
    {
        xx *= value; xy *= value; xz *= value;
        yx *= value; yy *= value; yz *= value;
        zx *= value; zy *= value; zz *= value;
    }
    
    ten33 operator/(const T &value) const
    {
        T inv = 1 / value;
        return ten33<T>(xx * inv, xy * inv, xz * inv,
                        yx * inv, yy * inv, yz * inv,
                        zx * inv, zy * inv, zz * inv);
    }
    
    void operator/=(const T &value)
    {
        T inv = 1 / value;
        xx *= inv; xy *= inv; xz *= inv;
        yx *= inv; yy *= inv; yz *= inv;
        zx *= inv; zy *= inv; zz *= inv;
    }
    
    bool operator==(const ten33 &t) const
    {
        return ((xx == t.xx) && (xy == t.xy) && (xz == t.xz) &&
                (yx == t.yx) && (yy == t.yy) && (yz == t.yz) &&
                (zx == t.zx) && (zy == t.zy) && (zz == t.zz));
    }
    
    bool operator!=(const ten33 &t) const
    {
        return !(*this == t);
    }
    
    ten33 operator-() const
    {
        return ten22<T>(-xx, -xy, -xz,
                        -yx, -yy, -yz,
                        -zx, -zy, -zz);
    }
    
    friend std::ostream& operator<<(std::ostream& os, ten33& t) {
        return os << "{xx:" << t.xx << " xy:" << t.xy << " xz:" << t.xz << "}" << std::endl
        << "{yx:" << t.yx << " yy:" << t.yy << " yz:" << t.yz << "}" << std::endl
        << "{zx:" << t.zx << " zy:" << t.zy << " zz:" << t.zz << "}";
    }
};

//
// ---------- tensor functions ---------------
//

// ---------- dyadic/outer product -----------
// vu = sum_i sum_j v_i u_j e_i e_j) - note vector u is transposed
template <typename T>
inline ten22<T> dyadic(const vec2<T> &v, const vec2<T> &u)
{
    return ten22<T>(v.x * u.x, v.x * u.y,
                    v.y * u.x, v.y * u.y);
}
template <typename T>
inline ten33<T> dyadic(const vec3<T> &v, const vec3<T> &u)
{
    return ten33<T>(v.x * u.x, v.x * u.y, v.x * u.z,
                    v.y * u.x, v.y * u.y, v.y * u.z,
                    v.z * u.x, v.z * u.y, v.z * u.z);
}

// ---------- tensor dot vector product -------
// t.v = sum_i e_i (sum_j t_ij v_j)
template <typename T>
inline vec2<T> dot(const ten22<T> &t, const vec2<T> &v)
{
    return vec2<T>(t.xx * v.x + t.xy * v.y,
                   t.yx * v.x + t.yy * v.y);
}
template <typename T>
inline vec3<T> dot(const ten33<T> &t, const vec3<T> &v)
{
    return vec3<T>(t.xx * v.x + t.xy * v.y + t.xz * v.z,
                   t.yx * v.x + t.yy * v.y + t.yz * v.z,
                   t.zx * v.x + t.zy * v.y + t.zz * v.z);
}

// ---------- vector dot tensor product -------
// v.t = sum_i e_i (sum_j v_j t_ji) - note vector v is transposed
template <typename T>
inline vec2<T> dot(const vec2<T> &v, const ten22<T> &t)
{
    return vec2<T>(v.x * t.xx + v.y * t.yx,
                   v.x * t.xy + v.y * t.yy);
}
template <typename T>
inline vec3<T> dot(const vec3<T> &v, const ten33<T> &t)
{
    return vec3<T>(v.x * t.xx + v.y * t.yx + v.z * t.zx,
                   v.x * t.xy + v.y * t.yy + v.z * t.zy,
                   v.x * t.xz + v.y * t.yz + v.z * t.zz);
}

// ---------- tensor dot product --------------
// t1.t2 = sum_i sum_l e_i e_l (sum_j t1_ij t2_jl)
template <typename T>
inline ten22<T> dot(const ten22<T> &t1, const ten22<T> &t2)
{
    return ten22<T>(t1.xx * t2.xx + t1.xy * t2.yx, t1.xx * t2.xy + t1.xy * t2.yy,
                    t1.yx * t2.xx + t1.yy * t2.yx, t1.yx * t2.xy + t1.yy * t2.yy);
}
template <typename T>
inline ten33<T> dot(const ten33<T> &t1, const ten33<T> &t2)
{
    return ten33<T>(t1.xx * t2.xx + t1.xy * t2.yx + t1.xz * t2.zx, t1.xx * t2.xy + t1.xy * t2.yy + t1.xz * t2.zy, t1.xx * t2.xz + t1.xy * t2.yz + t1.xz * t2.zz,
                    t1.yx * t2.xx + t1.yy * t2.yx + t1.yz * t2.zx, t1.yx * t2.xy + t1.yy * t2.yy + t1.yz * t2.zy, t1.yx * t2.xz + t1.yy * t2.yz + t1.yz * t2.zz,
                    t1.zx * t2.xx + t1.zy * t2.yx + t1.zz * t2.zx, t1.zx * t2.xy + t1.zy * t2.yy + t1.zz * t2.zy, t1.zx * t2.xz + t1.zy * t2.yz + t1.zz * t2.zz);
}

// ---------- tensor double dot product -------
// t1:t2 = sum_i sum_j t1_ij t2_ji
template <typename T>
inline T ddot(const ten22<T> &t1, const ten22<T> &t2)
{
    return (t1.xx * t2.xx + t1.xy * t2.yx +
            t1.yx * t2.xy + t1.yy * t2.yy);
}
template <typename T>
inline T ddot(const ten33<T> &t1, const ten33<T> &t2)
{
    return (t1.xx * t2.xx + t1.xy * t2.yx + t1.xz * t2.zx +
            t1.yx * t2.xy + t1.yy * t2.yy + t1.yz * t2.zy +
            t1.zx * t2.xz + t1.zy * t2.yz + t1.zz * t2.zz);
}

// ---------- tensor magnitude square ----------
// mag t_ij = 1/2 t:t^T = 1/2 sum_i sum_j t_ij : t_ij
template <typename T>
inline T mag2(const ten22<T> &t)
{
    ten22<T> temp = t;
    temp.transpose();
    return (0.5 * ddot<T>(t, temp));
}

template <typename T>
inline T mag2(const ten33<T> &t)
{
    ten33<T> temp = t;
    temp.transpose();
    return (0.5 * ddot<T>(t, temp));
}

// ---------- tensor magnitude -----------------
// mag t_ij = sqrt(1/2 t:t^T) = sqrt(1/2 sum_i sum_j t_ij : t_ij)
template <typename T>
inline T mag(const ten22<T> &t)
{
    ten22<T> temp = t;
    temp.transpose();
    return std::sqrt(0.5 * ddot<T>(t, temp));
}

template <typename T>
inline T mag(const ten33<T> &t)
{
    ten33<T> temp = t;
    temp.transpose();
    return std::sqrt(0.5 * ddot<T>(t, temp));
}

// ---------- tensor trace ---------------------
// trace t_ij = sum_i t_ii = I:t
template <typename T>
inline T trace(const ten22<T> &t)
{
    return (t.xx + t.yy);
}

template <typename T>
inline T trace(const ten33<T> &t)
{
    return (t.xx + t.yy + t.zz);
}

// ---------- tensor determinant ---------------
// det(t) = sum_i sum_j sum_k e_ijk t_i1 t_j2 t_k3
template <typename T>
inline T det(const ten22<T> &t)
{
    return  (t.xx * t.yy - t.yx * t.xy);
}
template <typename T>
inline T det(const ten33<T> &t)
{
    return  (t.xx * (t.yy * t.zz - t.zy * t.yz)) -
    (t.xy * (t.yx * t.zz - t.zx * t.yz)) +
    (t.xz * (t.yx * t.zy - t.zx * t.yy));
}

// ---------- tensor inverse -------------------
// Exist if and only if det(t) != 0
template <typename T>
inline ten22<T> inv(const ten22<T> &t)
{
    T invdet = 1.0 / det<T>(t);
    if (std::isinf(invdet))
    {
        std::cerr << "ERROR: Determinant of tensor is zero!" << std::endl;
        exit(1);
    }
    return ten22<T>(t.yy, -t.xy,
                    -t.yx,  t.xx)
    * invdet;
}

template <typename T>
inline ten33<T> inv(const ten33<T> &t)
{
    T invdet = 1.0 / det<T>(t);
    if (std::isinf(invdet))
    {
        std::cerr << "ERROR: Determinant of tensor is zero!" << std::endl;
        exit(1);
    }
    return ten33<T>((t.yy * t.zz - t.zy * t.yz), -(t.xy * t.zz - t.zy * t.xz),  (t.xy * t.yz - t.yy * t.xz),
                    -(t.yx * t.zz - t.zx * t.yz),  (t.xx * t.zz - t.zx * t.xz), -(t.xx * t.yz - t.yx * t.xz),
                    (t.yx * t.zy - t.zx * t.yy), -(t.xx * t.zy - t.zx * t.xy),  (t.xx * t.yy - t.yx * t.xy))
    * invdet;
}

#endif /* ten_h */
