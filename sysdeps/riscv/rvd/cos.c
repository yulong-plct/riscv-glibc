/* Double-precision vector cos function.

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

#include "v_math.h"


static const struct data
{
  vfloat64m2_t poly[7];
  vfloat64m2_t range_val, shift, inv_pi, half_pi, pi_1, pi_2, pi_3;
} data = {
  /* Worst-case error is 3.3 ulp in [-pi/2, pi/2].  */
  .poly = { V2 (-0x1.555555555547bp-3), V2 (0x1.1111111108a4dp-7),
	    V2 (-0x1.a01a019936f27p-13), V2 (0x1.71de37a97d93ep-19),
	    V2 (-0x1.ae633919987c6p-26), V2 (0x1.60e277ae07cecp-33),
	    V2 (-0x1.9e9540300a1p-41) },
  .inv_pi = V2 (0x1.45f306dc9c883p-2),
  .half_pi = V2 (0x1.921fb54442d18p+0),
  .pi_1 = V2 (0x1.921fb54442d18p+1),
  .pi_2 = V2 (0x1.1a62633145c06p-53),
  .pi_3 = V2 (0x1.c1cd129024e09p-106),
  .shift = V2 (0x1.8p52),
  .range_val = V2 (0x1p23)
};

#define C(i) d->poly[i]

static vfloat64m2_t NOINLINE
special_case (vfloat64m2_t x, vfloat64m2_t y, vuint64m2_t odd, vuint64m2_t cmp)
{
  y = vreinterpret_v_u64m2_f64m2 (vor (vreinterpret_v_f64m2_u64m2 (y), odd, 1));
  return v_call_f64 (cos, x, y, cmp);
}

vfloat64m2_t V_NAME_D1 (cos) (vfloat64m2_t x)
{
  const struct data *d = ptr_barrier (&data);
  vfloat64m2_t n, r, r2, r3, r4, t1, t2, t3, y;
  vuint64m2_t odd, cmp;

  r = vfabs_v_f64m2 (x, 2);
  cmp = (vuint64m2_t) vmsgeu (vreinterpret_v_f64m2_u64m2 (r),
		   vreinterpret_v_f64m2_u64m2 (d->range_val));
  if (__glibc_unlikely (v_any_u64 (cmp)))
    /* If fenv exceptions are to be triggered correctly, set any special lanes
       to 1 (which is neutral w.r.t. fenv). These lanes will be fixed by
       special-case handler later.  */
    r = vmsltu (cmp, v_f64 (1.0), r);

  /* n = rint((|x|+pi/2)/pi) - 0.5.  */
  n = vfmadd (d->shift, d->inv_pi, vfadd (r, d->half_pi,2), 2);
  odd = vshlq_n_u64 (vreinterpret_v_f64m2_u64m2 (n), 63);
  n = vfsub (n, d->shift, 2);
  n = vfsub (n, v_f64 (0.5), 2);

  /* r = |x| - n*pi  (range reduction into -pi/2 .. pi/2).  */
  r = vfmsub (r, d->pi_1, n, 2);
  r = vfmsub (r, d->pi_2, n, 2);
  r = vfmsub (r, d->pi_3, n, 2);

  /* sin(r) poly approx.  */
  r2 = vfmul (r, r, 2);
  r3 = vfmul (r2, r, 2);
  r4 = vfmul (r2, r2, 2);

  t1 = vfmadd (C (4), C (5), r2, 2);
  t2 = vfmadd (C (2), C (3), r2, 2);
  t3 = vfmadd (C (0), C (1), r2, 2);

  y = vfmadd (t1, C (6), r4, 2);
  y = vfmadd (t2, y, r4, 2);
  y = vfmadd (t3, y, r4, 2);
  y = vfmadd (r, y, r3, 2);

  if (__glibc_unlikely (v_any_u64 (cmp)))
    return special_case (x, y, odd, cmp);
  return vreinterpretq_f64_u64 (vor (vreinterpret_v_f64m2_u64m2 (y), odd, 2));
}
