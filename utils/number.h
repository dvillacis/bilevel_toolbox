/*
 * dticode/src/number.h
 *
 * Written 2012 by Tuomo Valkonen <tuomov@iki.fi>, University of Graz.
 *
 */


#ifndef _NUMBER_H
#define _NUMBER_H

#include <float.h>
#include <complex.h>
#include <math.h>

#if USE_FLOAT

typedef float FNumber;
typedef float complex FNumber_complex;

#define fn_max fmaxf
#define fn_min fminf
#define fn_abs fabsf
#define fn_sqrt sqrtf
#define fn_cexp cexpf
#define fn_exp expf
#define fn_cos cosf
#define fn_sin sinf
#define fn_real crealf
#define fn_imag cimagf
#define fn_atan2 atan2f
#define fn_log logf
#define fn_ceil ceilf
#if !_OPENACC
#define fn_cfromarg(x) cexpf(I*(x))
#define fn_cabs cabsf
#define fn_cabs2 cabs2f
#define fn_carg cargf
#define fn_conj conjf
#endif

#define FNUMBER(X) X
#define FNUMBER_MAX FLT_MAX
#define FNUMBER_MIN FLT_MIN
#define FNUMBER_EPS (FLT_EPSILON*0.5)
#define FNUMBER_RDX FLT_RADIX

#define FNUMBER_MPI_TYPE MPI_FLOAT

#else /* !USE_FLOAT */

typedef double FNumber;
typedef double complex FNumber_complex;

#define fn_max fmax
#define fn_min fmin
#define fn_abs fabs
#define fn_sqrt sqrt
#define fn_cexp cexp
#define fn_exp exp
#define fn_cos cos
#define fn_sin sin
#define fn_real creal
#define fn_imag cimag
#define fn_atan2 atan2
#define fn_log log
#define fn_ceil ceil
#if !_OPENACC
#define fn_cfromarg(x) cexp(I*(x))
#define fn_cabs cabs
#define fn_cabs2 cabs2
#define fn_carg carg
#define fn_conj conj
#endif

#define FNUMBER(X) X##f
#define FNUMBER_MAX DBL_MAX
#define FNUMBER_MIN DBL_MIN
#define FNUMBER_EPS (DBL_EPSILON*0.5)
#define FNUMBER_RDX FLT_RADIX

#define FNUMBER_MPI_TYPE MPI_DOUBLE

#endif /* USE_FLOAT */

#define FNUMBER_PREC (FNUMBER_RDX*FNUMBER_EPS)

/* FLT_MIN/DBL_MIN seems to agree with slamch('s'), dlamch('s') */
#define FNUMBER_SFMIN FNUMBER_MIN

#if _OPENACC
// PGCC does not support these intrinsics directly in accelerator regions
/*
#define fn_cfromarg(x) (cosf(x)+I*sinf(x))
#define fn_cabs(x) sqrt(crealf(x)*crealf(x)+cimagf(x)*cimagf(x))
#define fn_cabs2(x) (crealf(x)*crealf(x)+cimagf(x)*cimagf(x))
#define fn_carg(x) atan2f(cimagf(x), crealf(x))
*/

static __inline__ FNumber_complex fn_cfromarg(FNumber x)
{
    return fn_cos(x)+I*fn_sin(x);
}

static __inline__ FNumber fn_cabs2(FNumber_complex x)
{
    return fn_real(x)*fn_real(x)+fn_imag(x)*fn_imag(x);
}

static __inline__ FNumber fn_cabs(FNumber_complex x)
{
    return fn_sqrt(fn_cabs2(x));
}

static __inline__ FNumber fn_carg(FNumber_complex x)
{
    return fn_atan2(fn_imag(x), fn_real(x));
}

static __inline__ FNumber fn_conj(FNumber_complex x)
{
    return fn_real(x)-I*fn_imag(x);
}

#endif

static __inline__ FNumber fn_pow2(FNumber x)
{
    return x*x;
}

#endif /* _NUMBER_H */
