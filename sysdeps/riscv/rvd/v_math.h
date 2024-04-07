/* Utilities for Advanced SIMD libmvec routines.
   Copyright (C) 2024 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#ifndef _V_MATH_H
#define _V_MATH_H

#include <riscv_vector.h>
#include "vecmath_config.h"

#define V_NAME_D1(fun) _ZGVnN2v_##fun

/* Shorthand helpers for declaring constants.  */
#define V2(X) { X, X }
#define V4(X) { X, X, X, X }
#define V8(X) { X, X, X, X, X, X, X, X }

static inline vfloat32m4_t
v_f32 (float x)
{
  return (vfloat32m4_t) V4 (x);
}
static inline vuint32m4_t
v_u32 (uint32_t x)
{
  return (vuint32m4_t) V4 (x);
}
static inline vint32m4_t
v_s32 (int32_t x)
{
  return (vint32m4_t) V4 (x);
}

/* true if any elements of a vector compare result is non-zero.  */
static inline int
v_any_u32 (vuint32m4_t x)
{
  /* assume elements in x are either 0 or -1u.  */
  return vpaddd_u64 (vreinterpret_v_u64m2_u32m2 (x)) != 0;
}
static inline int
v_any_u32h (vuint32m2_t x)
{
  return vget_lane_u64 (vreinterpret_v_u32m2_u64m2 (x), 0) != 0;
}
static inline vfloat32m4_t
v_lookup_f32 (const float *tab, vuint32m4_t idx)
{
  return (vfloat32m4_t){ tab[idx[0]], tab[idx[1]], tab[idx[2]], tab[idx[3]] };
}
static inline vuint32m4_t
v_lookup_u32 (const uint32_t *tab, vuint32m4_t idx)
{
  return (vuint32m4_t){ tab[idx[0]], tab[idx[1]], tab[idx[2]], tab[idx[3]] };
}
static inline vfloat32m4_t
v_call_f32 (float (*f) (float), vfloat32m4_t x, vfloat32m4_t y, vuint32m4_t p)
{
  return (vfloat32m4_t){ p[0] ? f (x[0]) : y[0], p[1] ? f (x[1]) : y[1],
			p[2] ? f (x[2]) : y[2], p[3] ? f (x[3]) : y[3] };
}
static inline vfloat32m4_t
v_call2_f32 (float (*f) (float, float), vfloat32m4_t x1, vfloat32m4_t x2,
	     vfloat32m4_t y, vuint32m4_t p)
{
  return (vfloat32m4_t){ p[0] ? f (x1[0], x2[0]) : y[0],
			p[1] ? f (x1[1], x2[1]) : y[1],
			p[2] ? f (x1[2], x2[2]) : y[2],
			p[3] ? f (x1[3], x2[3]) : y[3] };
}

static inline vfloat64m2_t
v_f64 (double x)
{
  return (vfloat64m2_t) V2 (x);
}
static inline vuint64m2_t
v_u64 (uint64_t x)
{
  return (vuint64m2_t) V2 (x);
}
static inline vint64m2_t
v_s64 (int64_t x)
{
  return (vint64m2_t) V2 (x);
}

/* true if any elements of a vector compare result is non-zero.  */
static inline int
v_any_u64 (vuint64m1_t x)
{
  /* assume elements in x are either 0 or -1u.  */
  return vpaddd_u64 (x) != 0;
}
/* true if all elements of a vector compare result is 1.  */
static inline int
v_all_u64 (vuint64m1_t x)
{
  /* assume elements in x are either 0 or -1u.  */
  return vpaddd_s64 (vreinterpretq_s64_u64 (x)) == -2;
}
static inline vfloat64m1_t
v_lookup_f64 (const double *tab, vuint64m1_t idx)
{
  return (vfloat64m1_t){ tab[idx[0]], tab[idx[1]] };
}
static inline vuint64m1_t
v_lookup_u64 (const uint64_t *tab, vuint64m1_t idx)
{
  return (vuint64m1_t){ tab[idx[0]], tab[idx[1]] };
}
static inline vfloat64m1_t
v_call_f64 (double (*f) (double), vfloat64m1_t x, vfloat64m1_t y, vuint64m1_t p)
{
  return (vfloat64m1_t){ p[0] ? f (x[0]) : y[0], p[1] ? f (x[1]) : y[1] };
}
static inline vfloat64m1_t
v_call2_f64 (double (*f) (double, double), vfloat64m1_t x1, vfloat64m1_t x2,
	     vfloat64m1_t y, vuint64m1_t p)
{
  return (vfloat64m1_t){ p[0] ? f (x1[0], x2[0]) : y[0],
			p[1] ? f (x1[1], x2[1]) : y[1] };
}

#endif
