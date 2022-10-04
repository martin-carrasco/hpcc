#include <math.h>
#include <time.h>
#include <float.h>
/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* mpifft.c
 */



#ifndef HPL_NO_MPI_DATATYPE         /* Use MPI user-defined data type */
#define HPL_USE_MPI_DATATYPE
#endif
 
#ifndef HPL_COPY_L  /* do not copy L, use MPI user-defined data types */
#define HPL_NO_COPY_L
#endif
 
#ifndef HPL_DETAILED_TIMING         /* Do not enable detailed timings */
#define HPL_NO_DETAILED_TIMING
#endif
 
#ifndef HPL_CALL_VSIPL          /* Call the Fortran 77 BLAS interface */
#ifndef HPL_CALL_CBLAS                       /* there can be only one */
#define HPL_CALL_FBLAS
#endif
#endif

/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */

// hpl_misc.h
#ifdef __STDC__
#define HPL_STDC_HEADERS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HPL_STDC_HEADERS
#include <stdarg.h>
#define STDC_ARGS(p)           p
#else
#include <varargs.h>
#define STDC_ARGS(p)           ()
#endif

#ifdef HPL_CALL_VSIPL
#include <vsip.h>
#endif
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    HPL_rone             1.0
#define    HPL_rtwo             2.0
#define    HPL_rzero            0.0
/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#define    Mabs( a_ )          ( ( (a_) <   0  ) ? -(a_) : (a_) )
#define    Mmin( a_, b_ )      ( ( (a_) < (b_) ) ?  (a_) : (b_) )
#define    Mmax( a_, b_ )      ( ( (a_) > (b_) ) ?  (a_) : (b_) )

#define    Mfloor(a,b) (((a)>0) ? (((a)/(b))) : (-(((-(a))+(b)-1)/(b))))
#define    Mceil(a,b)           ( ( (a)+(b)-1 ) / (b) )
#define    Miceil(a,b) (((a)>0) ? ((((a)+(b)-1)/(b))) : (-((-(a))/(b))))

#define    Mupcase(C)          (((C)>96 && (C)<123) ? (C) & 0xDF : (C))
#define    Mlowcase(C)         (((C)>64 && (C)< 91) ? (C) | 32   : (C))
/*
 * Mptr returns a pointer to a_( i_, j_ ) for readability reasons and
 * also less silly errors ...
 */
#define    Mptr( a_, i_, j_, lda_ ) \
   ( (a_) + (size_t)(i_) + (size_t)(j_)*(size_t)(lda_) )
/*
 * Align pointer
 */
#define    HPL_PTR( ptr_, al_ ) \
                      ( ( ( (size_t)(ptr_)+(al_)-1 ) / (al_) ) * (al_) ) 


enum HPL_ORDER
{  HplRowMajor = 101,  HplColumnMajor  = 102 };
enum HPL_TRANS
{  HplNoTrans  = 111,  HplTrans        = 112,  HplConjTrans    = 113 };
enum HPL_UPLO
{  HplUpper    = 121,  HplLower        = 122 };
enum HPL_DIAG
{  HplNonUnit  = 131,  HplUnit         = 132 };
enum HPL_SIDE
{  HplLeft     = 141,  HplRight        = 142 }; 

#ifdef HPL_CALL_CBLAS
/*
 * ---------------------------------------------------------------------
 * The C interface of the BLAS is available ...
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    CBLAS_INDEX         int
 
#define    CBLAS_ORDER         HPL_ORDER
#define    CblasRowMajor       HplRowMajor
#define    CblasColMajor       HplColMajor
 
#define    CBLAS_TRANSPOSE     HPL_TRANS
#define    CblasNoTrans        HplNoTrans
#define    CblasTrans          HplTrans
#define    CblasConjTrans      HplConjTrans
 
#define    CBLAS_UPLO          HPL_UPLO
#define    CblasUpper          HplUpper
#define    CblasLower          HplLower
 
#define    CBLAS_DIAG          HPL_DIAG
#define    CblasNonUnit        HplNonUnit
#define    CblasUnit           HplUnit
 
#define    CBLAS_SIDE          HPL_SIDE
#define    CblasLeft           HplLeft
#define    CblasRight          HplRight
/*
 * ---------------------------------------------------------------------
 * CBLAS Function prototypes
 * ---------------------------------------------------------------------
 */
CBLAS_INDEX       cblas_idamax
STDC_ARGS(
(  const int,       const double *,  const int ) );
void              cblas_dswap
STDC_ARGS(
(  const int,       double *,        const int,       double *,
   const int ) );
void              cblas_dcopy
STDC_ARGS(
(  const int,       const double *,  const int,       double *,
   const int ) );
void              cblas_daxpy
STDC_ARGS(
(  const int,       const double,    const double *,  const int,
   double *,        const int ) );
void              cblas_dscal
STDC_ARGS(
(  const int,       const double,    double *,        const int ) );

void              cblas_dgemv
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,
   const int,       const int,       const double,    const double *,
   const int,       const double *,  const int,       const double,
   double *,        const int ) );

void              cblas_dger
STDC_ARGS(
(  const enum CBLAS_ORDER,           const int,       const int,
   const double,    const double *,  const int,       const double *,
   const int,       double *,        const int ) );
void              cblas_dtrsv
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_UPLO,
   const enum CBLAS_TRANSPOSE,       const enum CBLAS_DIAG,
   const int,       const double *,  const int,       double *,
   const int ) );

void              cblas_dgemm
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,
   const enum CBLAS_TRANSPOSE,       const int,       const int,
   const int,       const double,    const double *,  const int,
   const double *,  const int,       const double,    double *,
   const int ) );
void              cblas_dtrsm
STDC_ARGS(
(  const enum CBLAS_ORDER,           const enum CBLAS_SIDE,
   const enum CBLAS_UPLO,            const enum CBLAS_TRANSPOSE,
   const enum CBLAS_DIAG,            const int,       const int,
   const double,    const double *,  const int,       double *,
   const int ) );
/*
 * ---------------------------------------------------------------------
 * HPL C BLAS macro definition
 * ---------------------------------------------------------------------
 */
#define    HPL_dswap           cblas_dswap
#define    HPL_dcopy           cblas_dcopy
#define    HPL_daxpy           cblas_daxpy
#define    HPL_dscal           cblas_dscal
#define    HPL_idamax          cblas_idamax

#define    HPL_dgemv           cblas_dgemv
#define    HPL_dtrsv           cblas_dtrsv
#define    HPL_dger            cblas_dger

#define    HPL_dgemm           cblas_dgemm
#define    HPL_dtrsm           cblas_dtrsm

#endif

#ifdef HPL_CALL_FBLAS
/*
 * ---------------------------------------------------------------------
 * Use the Fortran 77 interface of the BLAS ...
 * ---------------------------------------------------------------------
 * Defaults: Add_, F77_INTEGER=int, StringSunStyle
 * ---------------------------------------------------------------------
 */
#ifndef NoChange
#ifndef UpCase
#ifndef Add__
#ifndef Add_

#define Add_

#endif
#endif
#endif
#endif

#ifndef F77_INTEGER
#define    F77_INTEGER         int
#else
#define    HPL_USE_F77_INTEGER_DEF
#endif

#ifndef StringCrayStyle
#ifndef StringStructVal
#ifndef StringStructPtr
#ifndef StringSunStyle

#define StringSunStyle

#endif
#endif
#endif
#endif
/*
 * ---------------------------------------------------------------------
 * Fortran 77 <-> C interface
 * ---------------------------------------------------------------------
 *
 * These macros identifies how Fortran routines will be called.
 *
 * Add_     : the Fortran compiler expects the name of C functions to be
 * in all lower case and to have an underscore postfixed it (Suns, Intel
 * compilers expect this).
 *
 * NoChange : the Fortran compiler expects the name of C functions to be
 * in all lower case (IBM RS6K compilers do this).
 *
 * UpCase   : the Fortran compiler expects the name of C functions to be
 * in all upcase. (Cray compilers expect this).
 *
 * Add__    : the Fortran compiler in use is f2c, a Fortran to C conver-
 * ter.
 */
#ifdef NoChange
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE DGEMM(...)          dgemm(...)
 */
#define    F77dswap               dswap
#define    F77dscal               dscal
#define    F77dcopy               dcopy
#define    F77daxpy               daxpy
#define    F77idamax              idamax

#define    F77dgemv               dgemv
#define    F77dtrsv               dtrsv
#define    F77dger                dger

#define    F77dgemm               dgemm
#define    F77dtrsm               dtrsm

#endif

#ifdef UpCase
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE DGEMM(...)          DGEMM(...)
 */
#ifdef CRAY_BLAS
                                                                                
#define    F77dswap               SSWAP
#define    F77dscal               SSCAL
#define    F77dcopy               SCOPY
#define    F77daxpy               SAXPY
#define    F77idamax              ISAMAX
                                                                                
#define    F77dgemv               SGEMV
#define    F77dtrsv               STRSV
#define    F77dger                SGER
                                                                                
#define    F77dgemm               SGEMM
#define    F77dtrsm               STRSM
                                                                                
#else

#define    F77dswap               DSWAP
#define    F77dscal               DSCAL
#define    F77dcopy               DCOPY
#define    F77daxpy               DAXPY
#define    F77idamax              IDAMAX

#define    F77dgemv               DGEMV
#define    F77dtrsv               DTRSV
#define    F77dger                DGER

#define    F77dgemm               DGEMM
#define    F77dtrsm               DTRSM

#endif

#endif

#ifdef Add_
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine  with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE DGEMM(...)          dgemm_(...)
 */
#define    F77dswap               dswap_
#define    F77dscal               dscal_
#define    F77dcopy               dcopy_
#define    F77daxpy               daxpy_
#define    F77idamax              idamax_

#define    F77dgemv               dgemv_
#define    F77dtrsv               dtrsv_
#define    F77dger                dger_

#define    F77dgemm               dgemm_
#define    F77dtrsm               dtrsm_

#endif

#ifdef Add__
/*
 * These defines  set  up  the  naming scheme required to have a FORTRAN
 * routine called by a C routine  with the following  FORTRAN to C inter-
 * face:
 *
 *          FORTRAN DECLARATION            C CALL
 *          SUBROUTINE DGEMM(...)          dgemm_(...)
 */
#define    F77dswap               dswap_
#define    F77dscal               dscal_
#define    F77dcopy               dcopy_
#define    F77daxpy               daxpy_
#define    F77idamax              idamax_
 
#define    F77dgemv               dgemv_
#define    F77dtrsv               dtrsv_
#define    F77dger                dger_
 
#define    F77dgemm               dgemm_
#define    F77dtrsm               dtrsm_
 
#endif
/*
 * ---------------------------------------------------------------------
 * Typedef definitions and conversion utilities
 * ---------------------------------------------------------------------
 */
#ifdef StringCrayStyle

#include <fortran.h>
                      /* Type of character argument in a FORTRAN call */
#define    F77_CHAR            _fcd
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(_fcdtocp(c) ))
#define    HPL_C2F_CHAR(c)     (_cptofcd(&(c), 1))

#define    F77_CHAR_DECL       F77_CHAR          /* input CHARACTER*1 */

#endif
/* ------------------------------------------------------------------ */
#ifdef StringStructVal
                      /* Type of character argument in a FORTRAN call */
typedef struct { char *cp; F77_INTEGER len; } F77_CHAR;
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(c.cp))

#define    F77_CHAR_DECL       F77_CHAR          /* input CHARACTER*1 */

#endif
/* ------------------------------------------------------------------ */
#ifdef StringStructPtr
                      /* Type of character argument in a FORTRAN call */
typedef struct { char *cp; F77_INTEGER len; } F77_CHAR;
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(c->cp))

#define    F77_CHAR_DECL       F77_CHAR *        /* input CHARACTER*1 */

#endif
/* ------------------------------------------------------------------ */
#ifdef StringSunStyle
                      /* Type of character argument in a FORTRAN call */
#define    F77_CHAR            char *
                                    /* Character conversion utilities */
#define    HPL_F2C_CHAR(c)     (*(c))
#define    HPL_C2F_CHAR(c)     (&(c))

#define    F77_CHAR_DECL       F77_CHAR          /* input CHARACTER*1 */
#define    F77_1_CHAR          , F77_INTEGER
#define    F77_2_CHAR          F77_1_CHAR F77_1_CHAR
#define    F77_3_CHAR          F77_2_CHAR F77_1_CHAR
#define    F77_4_CHAR          F77_3_CHAR F77_1_CHAR

#endif
/* ------------------------------------------------------------------ */

#ifndef F77_1_CHAR
#define    F77_1_CHAR
#define    F77_2_CHAR
#define    F77_3_CHAR
#define    F77_4_CHAR
#endif

#define    F77_INT_DECL        const F77_INTEGER *   /* input integer */
#define    F77_SIN_DECL        const double *         /* input scalar */
#define    F77_VIN_DECL        const double *         /* input vector */
#define    F77_VINOUT_DECL     double *        /* input/output matrix */
#define    F77_MIN_DECL        const double *         /* input matrix */
#define    F77_MINOUT_DECL     double *        /* input/output matrix */
 
#ifdef CRAY_PVP_ENV                      /* Type of FORTRAN functions */
#define    F77_VOID_FUN        extern fortran void      /* subroutine */
#define    F77_INT_FUN         extern fortran int /* integer function */
#else
#define    F77_VOID_FUN        extern void              /* subroutine */
#define    F77_INT_FUN         extern int         /* integer function */
#endif
/*
 * ---------------------------------------------------------------------
 * Fortran 77 BLAS function prototypes
 * ---------------------------------------------------------------------
 */
F77_VOID_FUN    F77dswap
STDC_ARGS(
(  F77_INT_DECL,    F77_VINOUT_DECL, F77_INT_DECL,    F77_VINOUT_DECL,
   F77_INT_DECL ) );
F77_VOID_FUN    F77dscal
STDC_ARGS(
(  F77_INT_DECL,    F77_SIN_DECL,    F77_VINOUT_DECL, F77_INT_DECL ) );
F77_VOID_FUN    F77dcopy
STDC_ARGS(
(  F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL,    F77_VINOUT_DECL,
   F77_INT_DECL ) );
F77_VOID_FUN    F77daxpy
STDC_ARGS(
(  F77_INT_DECL,    F77_SIN_DECL,    F77_VIN_DECL,    F77_INT_DECL,
   F77_VINOUT_DECL, F77_INT_DECL ) );
F77_INT_FUN     F77idamax
STDC_ARGS(
(  F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL ) );

F77_VOID_FUN    F77dgemv
STDC_ARGS(
(  F77_CHAR_DECL,   F77_INT_DECL,    F77_INT_DECL,    F77_SIN_DECL,
   F77_MIN_DECL,    F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL,
   F77_SIN_DECL,    F77_VINOUT_DECL, F77_INT_DECL     F77_1_CHAR ) );
F77_VOID_FUN    F77dger
STDC_ARGS(
(  F77_INT_DECL,    F77_INT_DECL,    F77_SIN_DECL,    F77_VIN_DECL,
   F77_INT_DECL,    F77_VIN_DECL,    F77_INT_DECL,    F77_MINOUT_DECL,
   F77_INT_DECL ) );
F77_VOID_FUN    F77dtrsv
STDC_ARGS(
(  F77_CHAR_DECL,   F77_CHAR_DECL,   F77_CHAR_DECL,   F77_INT_DECL,
   F77_MIN_DECL,    F77_INT_DECL,    F77_VINOUT_DECL, F77_INT_DECL
   F77_3_CHAR ) );

F77_VOID_FUN    F77dgemm
STDC_ARGS(
(  F77_CHAR_DECL,   F77_CHAR_DECL,   F77_INT_DECL,    F77_INT_DECL,
   F77_INT_DECL,    F77_SIN_DECL,    F77_MIN_DECL,    F77_INT_DECL,
   F77_MIN_DECL,    F77_INT_DECL,    F77_SIN_DECL,    F77_MINOUT_DECL,
   F77_INT_DECL     F77_2_CHAR ) );
F77_VOID_FUN    F77dtrsm
STDC_ARGS(
(  F77_CHAR_DECL,   F77_CHAR_DECL,   F77_CHAR_DECL,   F77_CHAR_DECL,
   F77_INT_DECL,    F77_INT_DECL,    F77_SIN_DECL,    F77_MIN_DECL,
   F77_INT_DECL,    F77_MINOUT_DECL, F77_INT_DECL     F77_4_CHAR ) );

#endif
/*
 * ---------------------------------------------------------------------
 * HPL BLAS Function prototypes
 * ---------------------------------------------------------------------
 */
#ifndef HPL_CALL_CBLAS

int                              HPL_idamax
STDC_ARGS( (
   const int,
   const double *,
   const int
) );
void                             HPL_daxpy
STDC_ARGS( (
   const int,
   const double,
   const double *,
   const int,
   double *,
   const int
) );
void                             HPL_dcopy
STDC_ARGS( (
   const int,
   const double *,
   const int,
   double *,
   const int
) );
void                             HPL_dscal
STDC_ARGS( (
   const int,
   const double,
   double *,
   const int
) );
void                             HPL_dswap
STDC_ARGS( (
   const int,
   double *,
   const int,
   double *,
   const int
) );
void                             HPL_dgemv
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_TRANS,
   const int,
   const int,
   const double,
   const double *,
   const int,
   const double *,
   const int,
   const double,
   double *,
   const int
) );
void                             HPL_dger
STDC_ARGS( (
   const enum HPL_ORDER,
   const int,
   const int,
   const double,
   const double *,
   const int,
   double *,
   const int,
   double *,
   const int
) );
void                             HPL_dtrsv
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_UPLO,
   const enum HPL_TRANS,
   const enum HPL_DIAG,
   const int,
   const double *,
   const int,
   double *,
   const int
) );
void                             HPL_dgemm
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_TRANS,
   const enum HPL_TRANS,
   const int,
   const int,
   const int,
   const double,
   const double *,
   const int,
   const double *,
   const int,
   const double,
   double *,
   const int
) );
void                             HPL_dtrsm
STDC_ARGS( (
   const enum HPL_ORDER,
   const enum HPL_SIDE,
   const enum HPL_UPLO,
   const enum HPL_TRANS,
   const enum HPL_DIAG,
   const int,
   const int,
   const double,
   const double *,
   const int,
   double *,
   const int
) );

#endif

typedef enum
{ HPL_NORM_A = 800, HPL_NORM_1 = 801, HPL_NORM_I = 802 } HPL_T_NORM;

typedef enum
{
   HPL_MACH_EPS   = 900,                /* relative machine precision */
   HPL_MACH_SFMIN = 901, /* safe minimum st 1/sfmin does not overflow */
   HPL_MACH_BASE  = 902,                /* base = base of the machine */
   HPL_MACH_PREC  = 903,                          /* prec  = eps*base */
   HPL_MACH_MLEN  = 904,   /* number of (base) digits in the mantissa */
   HPL_MACH_RND   = 905,        /* 1.0 if rounding occurs in addition */
   HPL_MACH_EMIN  = 906,   /* min exponent before (gradual) underflow */
   HPL_MACH_RMIN  = 907,        /* underflow threshold base**(emin-1) */
   HPL_MACH_EMAX  = 908,          /* largest exponent before overflow */
   HPL_MACH_RMAX  = 909  /* overflow threshold - (base**emax)*(1-eps) */
 
} HPL_T_MACH;
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_fprintf
STDC_ARGS( (
   FILE *,
   const char *,
   ...
) );
void                             HPL_warn
STDC_ARGS( (
   FILE *,
   int,
   const char *,
   const char *,
   ...
) );
void                             HPL_abort
STDC_ARGS( (
   int,
   const char *,
   const char *,
   ...
) );
void                             HPL_dlacpy
STDC_ARGS( (
   const int,
   const int,
   const double *,
   const int,
   double *,
   const int
) );
void                             HPL_dlatcpy
STDC_ARGS( (
   const int,
   const int,
   const double *,
   const int,
   double *,
   const int
) );
void                             HPL_dlaprnt
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const int,
   const int,
   const char *
) );
double                           HPL_dlange
STDC_ARGS( (
   const HPL_T_NORM,
   const int,
   const int,
   const double *,
   const int
) );
double                           HPL_dlamch
STDC_ARGS( (
   const HPL_T_MACH
) );

typedef enum
{
   HPL_LEFT_LOOKING  = 301,           /* Left looking lu fact variant */
   HPL_CROUT         = 302,                  /* Crout lu fact variant */
   HPL_RIGHT_LOOKING = 303           /* Right looking lu fact variant */
} HPL_T_FACT;
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void              HPL_dgesv
STDC_ARGS(
(  const int,       const int,       const int,       const HPL_T_FACT,
   const HPL_T_FACT,                 const int,       double *,
   const int,       int * ) );
void              HPL_ipid
STDC_ARGS(
(  const int,       double *,        int *,           int *,
   int *,           int *,           int *,           int *,
   const int,       const int,       const int,       const int,
   const int ) );

#include "mpi.h"

typedef enum { HPL_INT       = 100, HPL_DOUBLE       = 101 } HPL_T_TYPE;
 
typedef enum
{
   HPL_ROW_MAJOR     = 201,
   HPL_COLUMN_MAJOR  = 202
} HPL_T_ORDER;

typedef struct HPL_S_grid
{
   MPI_Comm        all_comm;                     /* grid communicator */
   MPI_Comm        row_comm;                      /* row communicator */
   MPI_Comm        col_comm;                   /* column communicator */
   HPL_T_ORDER     order;        /* ordering of the procs in the grid */
   int             iam;                        /* my rank in the grid */
   int             myrow;                /* my row number in the grid */
   int             mycol;             /* my column number in the grid */
   int             nprow;          /* the total # of rows in the grid */
   int             npcol;       /* the total # of columns in the grid */
   int             nprocs;        /* the total # of procs in the grid */
   int             row_ip2;          /* largest power of two <= nprow */
   int             row_hdim;     /* row_ip2 procs hypercube dimension */
   int             row_ip2m1;      /* largest power of two <= nprow-1 */
   int             row_mask;        /* row_ip2m1 procs hypercube mask */
   int             col_ip2;          /* largest power of two <= npcol */
   int             col_hdim;     /* col_ip2 procs hypercube dimension */
   int             col_ip2m1;      /* largest power of two <= npcol-1 */
   int             col_mask;        /* col_ip2m1 procs hypercube mask */
} HPL_T_grid;

/*
 * ---------------------------------------------------------------------
 * Data Structures
 * ---------------------------------------------------------------------
 */
typedef void (*HPL_T_OP)
(  const int,       const void *,    void *,          const HPL_T_TYPE );
/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#define    HPL_2_MPI_TYPE( typ ) \
                           ( ( typ == HPL_INT ? MPI_INT : MPI_DOUBLE ) )
/*
 * The following macros perform common modulo operations;  All functions
 * except MPosMod assume arguments are < d (i.e., arguments are themsel-
 * ves within modulo range).
 */
                                                /* increment with mod */
#define    MModInc(I, d)       if(++(I) == (d)) (I) = 0
                                                /* decrement with mod */
#define    MModDec(I, d)       if(--(I) == -1) (I) = (d)-1
                                                   /* positive modulo */
#define    MPosMod(I, d)       ( (I) - ((I)/(d))*(d) )
                                                   /* add two numbers */
#define    MModAdd(I1, I2, d) \
           ( ( (I1) + (I2) < (d) ) ? (I1) + (I2) : (I1) + (I2) - (d) )
                                                        /* add 1 to # */
#define    MModAdd1(I, d) ( ((I) != (d)-1) ? (I) + 1 : 0 )
                                              /* subtract two numbers */
#define    MModSub(I1, I2, d) \
           ( ( (I1) < (I2) ) ? (d) + (I1) - (I2) : (I1) - (I2) )
                                                      /* sub 1 from # */
#define    MModSub1(I, d) ( ((I)!=0) ? (I)-1 : (d)-1 )
/*
 * ---------------------------------------------------------------------
 * grid function prototypes
 * ---------------------------------------------------------------------
 */
int                              HPL_grid_init
STDC_ARGS( (
   MPI_Comm,
   const HPL_T_ORDER,
   const int,
   const int,
   HPL_T_grid *
) );
int                              HPL_grid_exit
STDC_ARGS( (
   HPL_T_grid *
) );

int                              HPL_grid_info
STDC_ARGS( (
   const HPL_T_grid *,
   int *,
   int *,
   int *,
   int *
) );
int                              HPL_pnum
STDC_ARGS( (
   const HPL_T_grid *,
   const int,
   const int
) );

int                              HPL_barrier
STDC_ARGS( (
   MPI_Comm
) );
int                              HPL_broadcast
STDC_ARGS( (
   void *,
   const int,
   const HPL_T_TYPE,
   const int,
   MPI_Comm
) );
int                              HPL_reduce
STDC_ARGS( (
   void *,
   const int,
   const HPL_T_TYPE,
   const HPL_T_OP ,
   const int,
   MPI_Comm
) );
int                              HPL_all_reduce
STDC_ARGS( (
   void *,
   const int,
   const HPL_T_TYPE,
   const HPL_T_OP ,
   MPI_Comm
) );

void                             HPL_max
STDC_ARGS( (
   const int,
   const void *,
   void *,
   const HPL_T_TYPE
) );
void                             HPL_min
STDC_ARGS( (
   const int,
   const void *,
   void *,
   const HPL_T_TYPE
) );
void                             HPL_sum
STDC_ARGS( (
   const int,
   const void *,
   void *,
   const HPL_T_TYPE
) );

#define    Mindxg2p( ig_, inb_, nb_, proc_, src_, nprocs_ )            \
           {                                                           \
              if( ( (ig_) >= (inb_) ) && ( (src_) >= 0 ) &&            \
                  ( (nprocs_) > 1 ) )                                  \
              {                                                        \
                 proc_  = (src_) + 1 + ( (ig_)-(inb_) ) / (nb_);       \
                 proc_ -= ( proc_ / (nprocs_) ) * (nprocs_);           \
              }                                                        \
              else                                                     \
              {                                                        \
                 proc_ = (src_);                                       \
              }                                                        \
           }

#define    Mindxg2l( il_, ig_, inb_, nb_, proc_, src_, nprocs_ )       \
           {                                                           \
              if( ( (ig_) < (inb_) ) || ( (src_) == -1 ) ||            \
                  ( (nprocs_) == 1 ) ) { il_ = (ig_); }                \
              else                                                     \
              {                                                        \
                 int i__, j__;                                         \
                 j__ = ( i__ = ( (ig_)-(inb_) ) / (nb_) ) / (nprocs_); \
                 il_ = (nb_)*( j__ - i__ ) +                           \
                       ( (i__ + 1 - ( j__ + 1 ) * (nprocs_) ) ?        \
                         (ig_) - (inb_) : (ig_) );                     \
              }                                                        \
           }

#define    Mindxg2lp( il_, proc_, ig_, inb_, nb_, src_, nprocs_ )      \
           {                                                           \
              if( ( (ig_) < (inb_) ) || ( (src_) == -1 ) ||            \
                  ( (nprocs_) == 1 ) )                                 \
              { il_ = (ig_); proc_ = (src_); }                         \
              else                                                     \
              {                                                        \
                 int i__, j__;                                         \
                 j__ = ( i__ = ( (ig_)-(inb_) ) / (nb_) ) / (nprocs_); \
                 il_ = (nb_)*(j__-i__) +                               \
                       ( ( i__ + 1 - ( j__ + 1 ) * (nprocs_) ) ?       \
                         (ig_) - (inb_) : (ig_) );                     \
                 proc_  = (src_) + 1 + i__;                            \
                 proc_ -= ( proc_ / (nprocs_) ) * (nprocs_);           \
              }                                                        \
           }
/*
 * Mindxl2g computes the global index ig_ corresponding to the local
 * index il_ in process proc_.
 */
#define    Mindxl2g( ig_, il_, inb_, nb_, proc_, src_, nprocs_ )       \
           {                                                           \
              if( ( (src_) >= 0 ) && ( (nprocs_) > 1 ) )               \
              {                                                        \
                 if( (proc_) == (src_) )                               \
                 {                                                     \
                    if( (il_) < (inb_) ) ig_ = (il_);                  \
                    else                 ig_ = (il_) +                 \
                       (nb_)*((nprocs_)-1)*(((il_)-(inb_))/(nb_) + 1); \
                 }                                                     \
                 else if( (proc_) < (src_) )                           \
                 {                                                     \
                    ig_ = (il_) + (inb_) +                             \
                          (nb_)*(  ((nprocs_)-1)*((il_)/(nb_)) +       \
                                   (proc_)-(src_)-1+(nprocs_) );       \
                 }                                                     \
                 else                                                  \
                 {                                                     \
                    ig_ =  (il_) + (inb_) +                            \
                           (nb_)*( ((nprocs_)-1)*((il_)/(nb_)) +       \
                           (proc_)-(src_)-1 );                         \
                 }                                                     \
              }                                                        \
              else                                                     \
              {                                                        \
                 ig_ = (il_);                                          \
              }                                                        \
           }
/*
 * MnumrocI computes the # of local indexes  np_ residing in the process
 * of coordinate  proc_  corresponding to the interval of global indexes
 * i_:i_+n_-1  assuming  that the global index 0 resides in  the process
 * src_,  and that the indexes are distributed from src_ using the para-
 * meters inb_, nb_ and nprocs_.
 */
#define    MnumrocI( np_, n_, i_, inb_, nb_, proc_, src_, nprocs_ )    \
           {                                                           \
              if( ( (src_) >= 0 ) && ( (nprocs_) > 1 ) )               \
              {                                                        \
                 int inb__, mydist__, n__, nblk__, quot__, src__;      \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 )                  \
                 {                                                     \
                    nblk__ = (-inb__) / (nb_) + 1;                     \
                    src__  = (src_) + nblk__;                          \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_);        \
                    inb__ += nblk__*(nb_);                             \
                    if( ( n__ = (n_) - inb__ ) <= 0 )                  \
                    {                                                  \
                       if( (proc_) == src__ ) np_ = (n_);              \
                       else                   np_ = 0;                 \
                    }                                                  \
                    else                                               \
                    {                                                  \
                       if( ( mydist__ = (proc_) - src__ ) < 0 )        \
                          mydist__ += (nprocs_);                       \
                       nblk__    = n__ / (nb_) + 1;                    \
                       mydist__ -= nblk__ -                            \
                          (quot__ = (nblk__ / (nprocs_))) * (nprocs_); \
                       if( mydist__ < 0 )                              \
                       {                                               \
                          if( (proc_) != src__ )                       \
                             np_ = (nb_) + (nb_) * quot__;             \
                          else                                         \
                             np_ = inb__ + (nb_) * quot__;             \
                       }                                               \
                       else if( mydist__ > 0 )                         \
                       {                                               \
                          np_ = (nb_) * quot__;                        \
                       }                                               \
                       else                                            \
                       {                                               \
                          if( (proc_) != src__ )                       \
                             np_ = n__ +(nb_)+(nb_)*(quot__ - nblk__); \
                          else                                         \
                             np_ = (n_)+      (nb_)*(quot__ - nblk__); \
                       }                                               \
                    }                                                  \
                 }                                                     \
                 else                                                  \
                 {                                                     \
                    if( ( n__ = (n_) - inb__ ) <= 0 )                  \
                    {                                                  \
                       if( (proc_) == (src_) ) np_ = (n_);             \
                       else                    np_ = 0;                \
                    }                                                  \
                    else                                               \
                    {                                                  \
                       if( ( mydist__ = (proc_) - (src_) ) < 0 )       \
                          mydist__ += (nprocs_);                       \
                       nblk__    = n__ / (nb_) + 1;                    \
                       mydist__ -= nblk__ -                            \
                          ( quot__ = (nblk__ / (nprocs_)) )*(nprocs_); \
                       if( mydist__ < 0 )                              \
                       {                                               \
                          if( (proc_) != (src_) )                      \
                             np_ = (nb_) + (nb_) * quot__;             \
                          else                                         \
                             np_ = inb__ + (nb_) * quot__;             \
                       }                                               \
                       else if( mydist__ > 0 )                         \
                       {                                               \
                          np_ = (nb_) * quot__;                        \
                       }                                               \
                       else                                            \
                       {                                               \
                          if( (proc_) != (src_) )                      \
                             np_ = n__ +(nb_)+(nb_)*(quot__ - nblk__); \
                          else                                         \
                             np_ = (n_)+      (nb_)*(quot__ - nblk__); \
                       }                                               \
                    }                                                  \
                 }                                                     \
              }                                                        \
              else                                                     \
              {                                                        \
                 np_ = (n_);                                           \
              }                                                        \
           }

#define    Mnumroc( np_, n_, inb_, nb_, proc_, src_, nprocs_ )         \
           MnumrocI( np_, n_, 0, inb_, nb_, proc_, src_, nprocs_ )
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_indxg2lp
STDC_ARGS( (
   int *,
   int *,
   const int,
   const int,
   const int,
   const int,
   const int
) );
int                              HPL_indxg2l
STDC_ARGS( (
   const int,
   const int,
   const int,
   const int,
   const int
) );
int                              HPL_indxg2p
STDC_ARGS( (
   const int,
   const int,
   const int,
   const int,
   const int
) );
int                              HPL_indxl2g
STDC_ARGS( (
   const int,
   const int,
   const int,
   const int,
   const int,
   const int
) );
void                             HPL_infog2l
STDC_ARGS( (
   int,
   int,
   const int,
   const int,
   const int,
   const int,
   const int,
   const int,
   const int,
   const int,
   const int,
   const int,
   int *,
   int *,
   int *,
   int *
) );
int                              HPL_numroc
STDC_ARGS( (
   const int,
   const int,
   const int,
   const int,
   const int,
   const int
) );
int                              HPL_numrocI
STDC_ARGS( (
   const int,
   const int,
   const int,
   const int,
   const int,
   const int,
   const int
) );

void                             HPL_dlaswp00N
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const int *
) );
void                             HPL_dlaswp10N
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const int *
) );
void                             HPL_dlaswp01N
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp01T
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp02N
STDC_ARGS( (
   const int,
   const int,
   const double *,
   const int,
   double *,
   double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp03N
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const double *,
   const double *,
   const int
) );
void                             HPL_dlaswp03T
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const double *,
   const double *,
   const int
) );
void                             HPL_dlaswp04N
STDC_ARGS( (
   const int,
   const int,
   const int,
   double *,
   const int,
   double *,
   const int,
   const double *,
   const double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp04T
STDC_ARGS( (
   const int,
   const int,
   const int,
   double *,
   const int,
   double *,
   const int,
   const double *,
   const double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp05N
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp05T
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const double *,
   const int,
   const int *,
   const int *
) );
void                             HPL_dlaswp06N
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   double *,
   const int,
   const int *
) );
void                             HPL_dlaswp06T
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   double *,
   const int,
   const int *
) );

void                             HPL_pabort
STDC_ARGS( (
   int,
   const char *,
   const char *,
   ...
) );
void                             HPL_pwarn
STDC_ARGS( (
   FILE *,
   int,
   const char *,
   const char *,
   ...
) );
void                             HPL_pdlaprnt
STDC_ARGS( (
   const HPL_T_grid *,
   const int,
   const int,
   const int,
   double *,
   const int,
   const int,
   const int,
   const char *
) );
double                           HPL_pdlamch
STDC_ARGS( (
   MPI_Comm,
   const HPL_T_MACH
) );
double                           HPL_pdlange
STDC_ARGS( (
   const HPL_T_grid *,
   const HPL_T_NORM,
   const int,
   const int,
   const int,
   const double *,
   const int
) );

typedef struct HPL_S_panel
{
   struct HPL_S_grid   * grid;             /* ptr to the process grid */
   struct HPL_S_palg   * algo;          /* ptr to the algo parameters */
   struct HPL_S_pmat   * pmat;         /* ptr to the local array info */
   double              * A;              /* ptr to trailing part of A */
   double              * WORK;                          /* work space */
   double              * L2;                              /* ptr to L */
   double              * L1;       /* ptr to jb x jb upper block of A */
   double              * DPIV;    /* ptr to replicated jb pivot array */
   double              * DINFO;      /* ptr to replicated scalar info */
   double              * U;                               /* ptr to U */
   int                 * IWORK;     /* integer workspace for swapping */
   void                * * * buffers[2];   /* buffers for panel bcast */
   int                 counts [2];          /* counts for panel bcast */
   MPI_Datatype        dtypes [2];      /* data types for panel bcast */
   MPI_Request         request[1];        /* requests for panel bcast */
   MPI_Status          status [1];          /* status for panel bcast */
   int                 nb;            /* distribution blocking factor */
   int                 jb;                             /* panel width */
   int                 m;   /* global # of rows of trailing part of A */
   int                 n;   /* global # of cols of trailing part of A */
   int                 ia;  /* global row index of trailing part of A */
   int                 ja;  /* global col index of trailing part of A */
   int                 mp;   /* local # of rows of trailing part of A */
   int                 nq;   /* local # of cols of trailing part of A */
   int                 ii;   /* local row index of trailing part of A */
   int                 jj;   /* local col index of trailing part of A */
   int                 lda;           /* local leading dim of array A */
   int                 prow;  /* proc. row owning 1st row of trail. A */
   int                 pcol;  /* proc. col owning 1st col of trail. A */
   int                 msgid;           /* message id for panel bcast */
   int                 ldl2;         /* local leading dim of array L2 */
   int                 len;      /* length of the buffer to broadcast */
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * Ablock;                           /* A block */
   vsip_block_d        * L1block;                         /* L1 block */
   vsip_block_d        * L2block;                         /* L2 block */
   vsip_block_d        * Ublock;                           /* U block */
#endif
} HPL_T_panel;

typedef enum
{
   HPL_1RING         = 401,                        /* Increasing ring */
   HPL_1RING_M       = 402,             /* Increasing ring (modified) */
   HPL_2RING         = 403,                      /* Increasing 2-ring */
   HPL_2RING_M       = 404,           /* Increasing 2-ring (modified) */
   HPL_BLONG         = 405,                         /* long broadcast */
   HPL_BLONG_M       = 406               /* long broadcast (modified) */
} HPL_T_TOP;
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    HPL_FAILURE            0
#define    HPL_SUCCESS            1
#define    HPL_KEEP_TESTING       2
/*
 * ---------------------------------------------------------------------
 * comm function prototypes
 * ---------------------------------------------------------------------
 */
int                              HPL_send
STDC_ARGS( (
   double *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_recv
STDC_ARGS( (
   double *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_sdrv
STDC_ARGS( (
   double *,
   int,
   int,
   double *,
   int,
   int,
   int,
   MPI_Comm
) );
int                              HPL_binit
STDC_ARGS( (
   HPL_T_panel *
) );
int                              HPL_bcast
STDC_ARGS( (
   HPL_T_panel *,
   int *
) );
int                              HPL_bwait
STDC_ARGS( (
   HPL_T_panel *
) );
int                              HPL_packL
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int
) );
void                             HPL_copyL
STDC_ARGS( (
   HPL_T_panel *
) );
 
int HPL_binit_1ring STDC_ARGS( ( HPL_T_panel *        ) );
int HPL_bcast_1ring STDC_ARGS( ( HPL_T_panel *, int * ) );
int HPL_bwait_1ring STDC_ARGS( ( HPL_T_panel *        ) );
 
int HPL_binit_1rinM STDC_ARGS( ( HPL_T_panel *        ) );
int HPL_bcast_1rinM STDC_ARGS( ( HPL_T_panel *, int * ) );
int HPL_bwait_1rinM STDC_ARGS( ( HPL_T_panel *        ) );
 
int HPL_binit_2ring STDC_ARGS( ( HPL_T_panel *        ) );
int HPL_bcast_2ring STDC_ARGS( ( HPL_T_panel *, int * ) );
int HPL_bwait_2ring STDC_ARGS( ( HPL_T_panel *        ) );
 
int HPL_binit_2rinM STDC_ARGS( ( HPL_T_panel *        ) );
int HPL_bcast_2rinM STDC_ARGS( ( HPL_T_panel *, int * ) );
int HPL_bwait_2rinM STDC_ARGS( ( HPL_T_panel *        ) );
 
int HPL_binit_blong STDC_ARGS( ( HPL_T_panel *        ) );
int HPL_bcast_blong STDC_ARGS( ( HPL_T_panel *, int * ) );
int HPL_bwait_blong STDC_ARGS( ( HPL_T_panel *        ) );
 
int HPL_binit_blonM STDC_ARGS( ( HPL_T_panel *        ) );
int HPL_bcast_blonM STDC_ARGS( ( HPL_T_panel *, int * ) );
int HPL_bwait_blonM STDC_ARGS( ( HPL_T_panel *        ) );


typedef void (*HPL_T_PFA_FUN)
(  HPL_T_panel *,   const int,       const int,       const int,
   double * );
typedef void (*HPL_T_RFA_FUN)
(  HPL_T_panel *,   const int,       const int,       const int,
   double * );
typedef void (*HPL_T_UPD_FUN)
(  HPL_T_panel *,   int *,           HPL_T_panel *,   const int ); 
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_dlocmax
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );

void                             HPL_dlocswpN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   double *
) );
void                             HPL_dlocswpT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   double *
) );
void                             HPL_pdmxswp
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );

void                             HPL_pdpancrN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdpancrT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdpanllN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdpanllT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdpanrlN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdpanrlT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );

void                             HPL_pdrpancrN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdrpancrT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdrpanllN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdrpanllT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdrpanrlN
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );
void                             HPL_pdrpanrlT
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int,
   const int,
   double *
) );

void                             HPL_pdfact
STDC_ARGS( (
   HPL_T_panel *
) );

typedef enum
{
   HPL_SWAP00        = 451,                      /* Use HPL_pdlaswp00 */
   HPL_SWAP01        = 452,                      /* Use HPL_pdlaswp01 */
   HPL_SW_MIX        = 453, /* Use HPL_pdlaswp00_ for small number of */
                            /* columns, and HPL_pdlaswp01_ otherwise. */
   HPL_NO_SWP        = 499
} HPL_T_SWAP;

typedef struct HPL_S_palg
{
   HPL_T_TOP           btopo;               /* row broadcast topology */
   int                 depth;                     /* look-ahead depth */
   int                 nbdiv;            /* recursive division factor */
   int                 nbmin;         /* recursion stopping criterium */
   HPL_T_FACT          pfact;                   /* panel fact variant */
   HPL_T_FACT          rfact;               /* recursive fact variant */
   HPL_T_PFA_FUN       pffun;              /* panel fact function ptr */
   HPL_T_RFA_FUN       rffun;          /* recursive fact function ptr */
   HPL_T_UPD_FUN       upfun;                      /* update function */
   HPL_T_SWAP          fswap;                   /* Swapping algorithm */
   int                 fsthr;                   /* Swapping threshold */
   int                 equil;                        /* Equilibration */
   int                 align;              /* data alignment constant */
} HPL_T_palg;

typedef struct HPL_S_pmat
{
#ifdef HPL_CALL_VSIPL
   vsip_block_d        * block;
#endif
   double              * A;            /* pointer to local piece of A */
   double              * X;             /* pointer to solution vector */
   int                 n;                      /* global problem size */
   int                 nb;                         /* blocking factor */
   int                 ld;                 /* local leading dimension */
   int                 mp;                    /* local number of rows */
   int                 nq;                 /* local number of columns */
   int                 info;                    /* computational flag */
} HPL_T_pmat;
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    MSGID_BEGIN_PFACT   1001              /* message id ranges */
#define    MSGID_END_PFACT     2000
#define    MSGID_BEGIN_FACT    2001
#define    MSGID_END_FACT      3000
#define    MSGID_BEGIN_PTRSV   3001
#define    MSGID_END_PTRSV     4000
 
#define    MSGID_BEGIN_COLL    9001
#define    MSGID_END_COLL     10000
/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#define    MNxtMgid( id_, beg_, end_ ) \
                             (( (id_)+1 > (end_) ?  (beg_) : (id_)+1 ))
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_pipid
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   int *
) );
void                             HPL_plindx0
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_pdlaswp00N
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );
void                             HPL_pdlaswp00T
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );

void                             HPL_perm
STDC_ARGS( (
   const int,
   int *,
   int *,
   int *
) );
void                             HPL_logsort
STDC_ARGS( (
   const int,
   const int,
   int *,
   int *,
   int *
) );
void                             HPL_plindx10
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int *,
   int *,
   int *,
   int *
) );
void                             HPL_plindx1
STDC_ARGS( (
   HPL_T_panel *,
   const int,
   const int *,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_spreadN
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const enum HPL_SIDE,
   const int,
   double *,
   const int,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_spreadT
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const enum HPL_SIDE,
   const int,
   double *,
   const int,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_equil
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const enum HPL_TRANS,
   const int,
   double *,
   const int,
   int *,
   const int *,
   const int *,
   int *
) );
void                             HPL_rollN
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int,
   double *,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_rollT
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int,
   double *,
   const int,
   const int *,
   const int *,
   const int *
) );
void                             HPL_pdlaswp01N
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );
void                             HPL_pdlaswp01T
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );

void                             HPL_pdupdateNN
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );
void                             HPL_pdupdateNT
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );
void                             HPL_pdupdateTN
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );
void                             HPL_pdupdateTT
STDC_ARGS( (
   HPL_T_panel *,
   int *,
   HPL_T_panel *,
   const int
) );

void                             HPL_pdgesv0
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pmat *
) );
void                             HPL_pdgesvK1
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pmat *
) );
void                             HPL_pdgesvK2
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pmat *
) );
void                             HPL_pdgesv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   HPL_T_pmat *
) );
 
void                             HPL_pdtrsv
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_pmat *
) );

void                             HPL_pdpanel_new
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   const int,
   HPL_T_pmat *,
   const int,
   const int,
   const int,
   HPL_T_panel * *
) );
void                             HPL_pdpanel_init
STDC_ARGS( (
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   const int,
   HPL_T_pmat *,
   const int,
   const int,
   const int,
   HPL_T_panel *
) );
int                              HPL_pdpanel_disp
STDC_ARGS( (
   HPL_T_panel * *
) );
int                              HPL_pdpanel_free
STDC_ARGS( (
   HPL_T_panel *
) );

#define    HPL_NTIMER              64
#define    HPL_TIMER_STARTFLAG    5.0
#define    HPL_TIMER_ERROR       -1.0
/*
 * ---------------------------------------------------------------------
 * type definitions
 * ---------------------------------------------------------------------
 */
typedef enum
{  HPL_WALL_TIME = 101, HPL_CPU_TIME  = 102 } HPL_T_TIME;
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
double          HPL_timer_cputime    STDC_ARGS(     ( void      ) );
double          HPL_timer_walltime   STDC_ARGS(     ( void      ) );

void            HPL_timer            STDC_ARGS(     ( const int ) );
void            HPL_timer_boot       STDC_ARGS(     ( void      ) );
void            HPL_timer_enable     STDC_ARGS(     ( void      ) );
void            HPL_timer_disable    STDC_ARGS(     ( void      ) );
double          HPL_timer_inquire
STDC_ARGS(
(  const HPL_T_TIME,                 const int ) );


#define    HPL_MULT0         1284865837
#define    HPL_MULT1         1481765933
#define    HPL_IADD0         1
#define    HPL_IADD1         0
#define    HPL_DIVFAC        2147483648.0
#define    HPL_POW16         65536.0
#define    HPL_HALF          0.5
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_dmatgen
STDC_ARGS( (
   const int,
   const int,
   double *,
   const int,
   const int
) );
void                             HPL_lmul
STDC_ARGS( (
   int *,
   int *,
   int *
) );
void                             HPL_ladd
STDC_ARGS( (
   int *,
   int *,
   int *
) );
void                             HPL_xjumpm
STDC_ARGS( (
   const int,
   int *,
   int *,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_setran
STDC_ARGS( (
   const int,
   int *
) );
void                             HPL_jumpit
STDC_ARGS( (
   int *,
   int *,
   int *,
   int *
) );
double                           HPL_rand STDC_ARGS( ( void ) );


void            HPL_dinfo
STDC_ARGS(
(  FILE * *,        int *,           int *,           int *,
   HPL_T_FACT *,    int *,           int *,           int *, 
   int *,           int *,           HPL_T_FACT *,    int *,
   double *,        double * ) );
void            HPL_dtest
STDC_ARGS(
(  FILE *,          const int,       const int,       const int,
   HPL_T_FACT,      HPL_T_FACT,      const int,       const double,
   const double,    int *,           int *,           int * ) );

#define    HPL_NPTIMER             64
#define    HPL_PTIMER_STARTFLAG   5.0
#define    HPL_PTIMER_ERROR      -1.0
/*
 * ---------------------------------------------------------------------
 * type definitions
 * ---------------------------------------------------------------------
 */
typedef enum
{  HPL_WALL_PTIME = 101, HPL_CPU_PTIME  = 102 } HPL_T_PTIME;

typedef enum
{ HPL_AMAX_PTIME  = 201, HPL_AMIN_PTIME = 202, HPL_SUM_PTIME  = 203 }
HPL_T_PTIME_OP;
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
double          HPL_ptimer_cputime   STDC_ARGS(     ( void      ) );
double          HPL_ptimer_walltime  STDC_ARGS(     ( void      ) );

void            HPL_ptimer           STDC_ARGS(     ( const int ) );
void            HPL_ptimer_boot      STDC_ARGS(     ( void      ) );
void            HPL_ptimer_combine
STDC_ARGS(
(  MPI_Comm comm,   const HPL_T_PTIME_OP,             const HPL_T_PTIME,
   const int,       const int,       double * ) );
void            HPL_ptimer_disable   STDC_ARGS(     ( void      ) );
void            HPL_ptimer_enable    STDC_ARGS(     ( void      ) );
double          HPL_ptimer_inquire
STDC_ARGS(
(  const HPL_T_PTIME,                const int ) );

void                             HPL_pdmatgen
STDC_ARGS( (
   const HPL_T_grid *,
   const int,
   const int,
   const int,
   double *,
   const int,
   const int
) );

typedef struct HPL_S_test
{
   double              epsil;                      /* epsilon machine */
   double              thrsh;                            /* threshold */
   FILE *              outfp;       /* output stream (only in proc 0) */
   int                 kfail;                    /* # of tests failed */
   int                 kpass;                    /* # of tests passed */
   int                 kskip;                   /* # of tests skipped */
   int                 ktest;                /* total number of tests */
} HPL_T_test;

typedef struct {
  double Gflops, time, eps, RnormI, Anorm1, AnormI, Xnorm1, XnormI, BnormI;
  int N, NB, nprow, npcol, depth, nbdiv, nbmin;
  char cpfact, crfact, ctop, order;
} HPL_RuntimeData;

/*
 * ---------------------------------------------------------------------
 * #define macro constants for testing only
 * ---------------------------------------------------------------------
 */
#define    HPL_LINE_MAX         256
#define    HPL_MAX_PARAM         20
#define    HPL_ISEED            100
/*
 * ---------------------------------------------------------------------
 * global timers for timing analysis only
 * ---------------------------------------------------------------------
 */
#ifdef HPL_DETAILED_TIMING
#define    HPL_TIMING_BEG        11 /* timer 0 reserved, used by main */
#define    HPL_TIMING_N           6 /* number of timers defined below */
#define    HPL_TIMING_RPFACT     11 /* starting from here, contiguous */
#define    HPL_TIMING_PFACT      12
#define    HPL_TIMING_MXSWP      13
#define    HPL_TIMING_UPDATE     14
#define    HPL_TIMING_LASWP      15
#define    HPL_TIMING_PTRSV      16
#endif
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void                             HPL_pdinfo
STDC_ARGS( (
   HPL_T_test *,
   int *,
   int *,
   int *,
   int *,
   HPL_T_ORDER *,
   int *,
   int *,
   int *,
   int *,
   HPL_T_FACT *,
   int *,
   int *,
   int *,
   int *,
   int *,
   HPL_T_FACT *,
   int *,
   HPL_T_TOP *,
   int *,
   int *,
   HPL_T_SWAP *,
   int *,
   int *,
   int *,
   int *,
   int *
) );
void                             HPL_pdtest
STDC_ARGS( (
   HPL_T_test *,
   HPL_T_grid *,
   HPL_T_palg *,
   const int,
   const int,
   HPL_RuntimeData *
) );

/*
 * End of hpl.h
 */

#define HPCC_VERSION_MAJOR   1
#define HPCC_VERSION_MINOR   5
#define HPCC_VERSION_MICRO   0
#define HPCC_VERSION_RELEASE 'f'





#define MPIFFT_TIMING_COUNT 8

/* Define 64-bit types and corresponding format strings for printf() and constants */
#ifdef LONG_IS_64BITS
typedef unsigned long u64Int;
typedef long s64Int;
#define FSTR64 "%ld"
#define FSTRU64 "%lu"
#define ZERO64B 0L
#else
typedef unsigned long long u64Int;
typedef long long s64Int;
#define FSTR64 "%lld"
#define FSTRU64 "%llu"
#define ZERO64B 0LL
#endif

typedef struct {
  double GBs, time, residual;
  int n, nb, nprow, npcol;
} PTRANS_RuntimeData;

/* parameters of execution */
typedef struct {
  /* HPL section */
   HPL_T_test                 test;
   int                        nval  [HPL_MAX_PARAM],
                              nbval [HPL_MAX_PARAM],
                              pval  [HPL_MAX_PARAM],
                              qval  [HPL_MAX_PARAM],
                              nbmval[HPL_MAX_PARAM],
                              ndvval[HPL_MAX_PARAM],
                              ndhval[HPL_MAX_PARAM];
   HPL_T_ORDER                porder;
   HPL_T_FACT                 pfaval[HPL_MAX_PARAM],
                              rfaval[HPL_MAX_PARAM];
   HPL_T_TOP                  topval[HPL_MAX_PARAM];
   HPL_T_FACT                 rpfa;
   HPL_T_SWAP                 fswap;
   int ns, nbs, npqs, npfs, nbms, ndvs, nrfs, ntps, ndhs, tswap, L1notran, Unotran, equil, align;

  /* HPCC section */
  char inFname[256 + 1], outFname[256 + 1];
  int PTRANSns, PTRANSnval[2 * HPL_MAX_PARAM];
  int PTRANSnbs, PTRANSnbval[2 * HPL_MAX_PARAM];
  int PTRANSnpqs, PTRANSpval[2 * HPL_MAX_PARAM], PTRANSqval[2 * HPL_MAX_PARAM];
  double MPIRandomAccess_LCG_GUPs, MPIRandomAccess_GUPs, Star_LCG_GUPs, Single_LCG_GUPs, StarGUPs, SingleGUPs,
    MPIRandomAccess_ErrorsFraction, MPIRandomAccess_time, MPIRandomAccess_CheckTime,
    MPIRandomAccess_TimeBound,
    MPIRandomAccess_LCG_ErrorsFraction, MPIRandomAccess_LCG_time, MPIRandomAccess_LCG_CheckTime,
    MPIRandomAccess_LCG_TimeBound,
    StarStreamCopyGBs, StarStreamScaleGBs,
    StarStreamAddGBs, StarStreamTriadGBs, SingleStreamCopyGBs, SingleStreamScaleGBs,
    SingleStreamAddGBs, SingleStreamTriadGBs, StarDGEMMGflops, SingleDGEMMGflops;
  double StarFFTGflops, SingleFFTGflops, MPIFFTGflops, MPIFFT_maxErr;
  double MaxPingPongLatency, RandomlyOrderedRingLatency, MinPingPongBandwidth,
    NaturallyOrderedRingBandwidth, RandomlyOrderedRingBandwidth,
    MinPingPongLatency, AvgPingPongLatency, MaxPingPongBandwidth, AvgPingPongBandwidth,
    NaturallyOrderedRingLatency;
  int DGEMM_N;
  int StreamThreads, StreamVectorSize;
  int FFT_N;
  int MPIFFT_Procs;
  int MPIRandomAccess_LCG_Algorithm, MPIRandomAccess_Algorithm;

  HPL_RuntimeData HPLrdata;
  PTRANS_RuntimeData PTRANSrdata;

  int Failure; /* over all failure of the benchmark */

  double MPIFFTtimingsForward[MPIFFT_TIMING_COUNT], MPIFFTtimingsBackward[MPIFFT_TIMING_COUNT];

  size_t HPLMaxProcMem;
  int HPLMaxProc, HPLMinProc;
  int RunHPL, RunStarDGEMM, RunSingleDGEMM,
    RunPTRANS, RunStarStream, RunSingleStream,
    RunMPIRandomAccess_LCG, RunStarRandomAccess_LCG, RunSingleRandomAccess_LCG,
    RunMPIRandomAccess, RunStarRandomAccess, RunSingleRandomAccess,
    RunStarFFT, RunSingleFFT, RunMPIFFT,
    RunLatencyBandwidth;

  int FFTEnblk, FFTEnp, FFTEl2size;
  s64Int RandomAccess_LCG_N, RandomAccess_N, MPIRandomAccess_LCG_ExeUpdates, MPIRandomAccess_ExeUpdates,
    MPIRandomAccess_LCG_N, MPIRandomAccess_N, MPIRandomAccess_LCG_Errors, MPIRandomAccess_Errors, MPIFFT_N;
} HPCC_Params;
/*
This is what needs to be done to add a new benchmark:
-  Add the benchmark code to the directory structure (headers, makefiles)
-  Add benchmark output data to the HPCC_Params structure.
-  Initialize the HPCC_Params structure data in HPCC_Init().
-  Add a call to the benchmark function in main().
-  Make sure that all the processes fill out the structure with the same data.
-  Print the output of the benchmark in HPCC_Finalize().
-  For tests that have "Star" and "Single" variants (DGEMM, RandomAccess, STREAM) the function
that performs the test returns a value (0 or 1) that indicates runtime failure and also returns
benchamark failure (due to wrong optimization that causes numerical error) by setting
params->Failure.
*/

/*
int HPCC_external_init(int argc, char *argv[], void *extdata);
int HPCC_external_finalize(int argc, char *argv[], void *extdata);

extern int HPCC_Init(HPCC_Params *params);
extern int HPCC_Finalize(HPCC_Params *params);
*/
int HPCC_LocalVectorSize(HPCC_Params *params, int vecCnt, size_t size, int pow2);


extern int
HPCC_Defaults(HPL_T_test *TEST, int *NS, int *N,
              int *NBS, int *NB,
              HPL_T_ORDER *PMAPPIN,
              int *NPQS, int *P, int *Q,
              int *NPFS, HPL_T_FACT *PF,
              int *NBMS, int *NBM,
              int *NDVS, int *NDV,
              int *NRFS, HPL_T_FACT *RF,
              int *NTPS, HPL_T_TOP *TP,
              int *NDHS, int *DH,
              HPL_T_SWAP *FSWAP, int *TSWAP, int *L1NOTRAN, int *UNOTRAN, int *EQUIL, int *ALIGN, MPI_Comm comm);

extern int HPL_main(int ARGC, char **ARGV, HPL_RuntimeData *rdata, int *failure);
extern float HPL_slamch (const HPL_T_MACH);
extern double HPCC_dweps();
extern float HPCC_sweps();
extern int HPCC_StarDGEMM(HPCC_Params *params);
extern int HPCC_SingleDGEMM(HPCC_Params *params);
extern int PTRANS(HPCC_Params *params);
extern int HPCC_MPIRandomAccess_LCG(HPCC_Params *params);
extern int HPCC_SingleRandomAccess_LCG(HPCC_Params *params);
extern int HPCC_StarRandomAccess_LCG(HPCC_Params *params);
extern int HPCC_MPIRandomAccess(HPCC_Params *params);
extern int HPCC_SingleRandomAccess(HPCC_Params *params);
extern int HPCC_StarRandomAccess(HPCC_Params *params);
extern int HPCC_SingleStream(HPCC_Params *params);
extern int HPCC_StarStream(HPCC_Params *params);
extern int HPCC_StarFFT(HPCC_Params *params);
extern int HPCC_SingleFFT(HPCC_Params *params);
extern int HPCC_MPIFFT(HPCC_Params *params);

extern int HPCC_TestFFT(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure);
extern int HPCC_TestDGEMM(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure);
extern int MaxMem(int nprocs, int imrow, int imcol, int nmat, int *mval, int *nval, int nbmat,
  int *mbval, int *nbval, int ngrids, int *npval, int *nqval, long *maxMem);
extern int HPCC_Stream(HPCC_Params *params, int doIO, MPI_Comm comm, int world_rank,
  double *copyGBs, double *scaleGBs, double *addGBs, double *triadGBs,
  int *failure);
extern void main_bench_lat_bw(HPCC_Params *params);

extern int pdtrans(char *trans, int *m, int *n, int * mb, int *nb, double *a, int *lda,
  double *beta, double *c__, int *ldc, int *imrow, int *imcol, double *work, int *iwork);
extern FILE* pdtransinfo(int *nmat, int *mval, int *nval, int *ldval,
  int *nbmat, int *mbval, int *nbval, int *ldnbval, int *ngrids, int *npval, int *nqval,
  int *ldpqval, int *iaseed, int *imrow, int *imcol, float *thresh, int *iam, int *nprocs,
  double *eps, char *infname, int *fcl, char *outfname);
extern int pdmatgen(int *ictxt, char *aform, char *diag, int *m, int *n, int *mb, int *nb, double*a,
  int *lda, int *iarow, int *iacol, int *iseed, int *iroff, int *irnum, int *icoff, int *icnum,
  int * myrow, int *mycol, int *nprow, int *npcol, double alpha);
extern int pdmatcmp(int *ictxt, int *m_, int *n_, double *a, int *lda_, double *aCopy, int *ldc_,
  double *error);
extern int pxerbla(int *ictxt, char *srname, int *info);
extern int slcombine_(int *ictxt, char *scope, char *op, char * timetype, int *n, int *ibeg,
  double *times);
extern int icopy_(int *n, int *sx, int *incx, int * sy, int *incy);
extern int numroc_(int *, int *, int *, int *, int *);
extern int slboot_(void);
extern int sltimer_(int *i__);
extern int ilcm_(int *, int *);
extern int iceil_(int *, int *);
extern double pdrand();
extern int setran_(int *, int *, int *);
extern int jumpit_(int *, int *, int *, int *);
extern int xjumpm_(int *, int *, int *, int *, int *, int *, int *);


/* ---------------------------------------------------------------------- */

#define DPRN(i,v) do{printf(__FILE__ "(%d)@%d:" #v "=%g\n",__LINE__,i,(double)(v));fflush(stdout);}while(0)

#define BEGIN_IO(r,fn,f) if(0==r){f=fopen(fn,"a");if(!f)fprintf(f=stderr,"Problem with appending to file '%s'\n",fn)
#define END_IO(r,f) fflush(f); if (f!=stdout && f!=stderr) fclose(f);} f=(FILE*)(NULL)

#ifdef HPCC_MEMALLCTR
extern int HPCC_alloc_init(size_t total_size);
extern int HPCC_alloc_finalize();
extern void *HPCC_malloc(size_t size);
extern void HPCC_free(void *ptr);
#define HPCC_fftw_malloc HPCC_malloc
#define HPCC_fftw_free HPCC_free
#define HPCC_XMALLOC(t,s) ((t*)HPCC_malloc(sizeof(t)*(s)))
#else
#define HPCC_malloc malloc
#define HPCC_free free
#define HPCC_fftw_malloc fftw_malloc
#define HPCC_fftw_free fftw_free
#define HPCC_XMALLOC(t,s) XMALLOC(t,s)
#endif

#define XMALLOC(t,s) ((t*)malloc(sizeof(t)*(s)))
#define XCALLOC(t,s) ((t*)calloc((s),sizeof(t)))


#define FFTE_NDA2 65536
#define FFTE_NDA3 4096
#define FFTE_NDA4 256

/* Parameters that affect performance */

/*
  Blocking parameter. Suggested values:
   8 for Pentium III and Athlon
  16 for Pentium4, Athlon XP, Opteron, Itanium and Itanium2
*/
#ifndef FFTE_NBLK
#define FFTE_NBLK 16
#endif

/*
  Padding parameter to avoid cache conflicts.
  Suggested values:
  2 for Pentium III
  4 for Athlon, Athlon XP, Opteron, Itanium
  8 for Pentium4 and Itanium2
*/
#ifndef FFTE_NP
#define FFTE_NP 8
#endif

/* Size of Level 2 cache */
#ifndef FFTE_L2SIZE
#define FFTE_L2SIZE 1048576
#endif

#ifdef LONG_IS_64BITS
typedef unsigned long u64Int_t;
typedef long s64Int_t;
#else
typedef unsigned long long u64Int_t;
typedef long long s64Int_t;
#endif


#ifdef USING_FFTW

#include <fftw.h>

#else

typedef double fftw_real;
typedef struct {
     fftw_real re, im;
} fftw_complex_orig;
typedef fftw_real HPCC_Complex[2];
typedef HPCC_Complex fftw_complex;

typedef enum {
     FFTW_FORWARD = -1, FFTW_BACKWARD = 1
} fftw_direction;
#endif

struct hpcc_fftw_plan_struct {
  fftw_complex *w1, *w2, *ww1, *ww2, *ww3, *ww4, *c, *d;
  int n, c_size, d_size;
  int flags;
  fftw_direction dir;
};
typedef struct hpcc_fftw_plan_struct *hpcc_fftw_plan;

extern hpcc_fftw_plan HPCC_fftw_create_plan(int n, fftw_direction dir, int flags);
extern void HPCC_fftw_destroy_plan(hpcc_fftw_plan plan);
extern void HPCC_fftw_one(hpcc_fftw_plan plan, fftw_complex *in, fftw_complex *out);

extern int HPCC_ipow(int x, int p);

extern int HPCC_zfft1d(int n, fftw_complex *a, fftw_complex *b, int iopt, hpcc_fftw_plan p);

int HPCC_fft235(fftw_complex *a, fftw_complex *b, fftw_complex *w, int n, const int *ip);
int HPCC_settbl(fftw_complex *w, int n);

int HPCC_factor235(int n, int *ip);
int HPCC_factor235_8(s64Int_t n, int *ip);

extern int HPCC_bcnrand(u64Int_t n, u64Int_t a, void *x);

#define ARR2D(a, i, j, lda) a[(i)+(j)*(lda)]
#define PTR2D(a, i, j, lda) (a+(i)+(j)*(lda))
#define ARR3D(a, i, j, k, lda1, lda2) a[(i)+(lda1)*((j)+(k)*(lda2))]
#define PTR3D(a, i, j, k, lda1, lda2) (a+(i)+(lda1)*((j)+(k)*(lda2)))
#define ARR4D(a, i, j, k, l, lda1, lda2, lda3) a[(i)+(lda1)*((j)+(lda2)*((k)+(lda3)*(l)))]
#define c_mul3v(v, v1, v2) c_re(v) = c_re(v1) * c_re(v2) - c_im(v1) * c_im(v2); c_im(v) = c_re(v1) * c_im(v2) + c_im(v1) * c_re(v2)
#define c_assgn(d, s) c_re(d)=c_re(s);c_im(d)=c_im(s)
#define V3MIN(r, e, v) r = (e); V2MIN(r, v)
#define V2MIN(r, v) r = (v) < r ? (v) : r
#define EMAX(d, v, e) d=(e); d=d>(v)?d:(v)

#define    Mmax( a_, b_ )      ( ( (a_) > (b_) ) ?  (a_) : (b_) )


#ifndef USING_FFTW

typedef struct hpcc_fftw_plan_struct *fftw_plan;

#define c_re(c)  ((c)[0])
#define c_im(c)  ((c)[1])

#define fftw_malloc malloc
#define fftw_free free
/* flags for the planner */
#define  FFTW_ESTIMATE (0)
#define  FFTW_MEASURE  (1)

#define FFTW_OUT_OF_PLACE (0)
#define FFTW_IN_PLACE (8)
#define FFTW_USE_WISDOM (16)

#define fftw_create_plan HPCC_fftw_create_plan
#define fftw_destroy_plan HPCC_fftw_destroy_plan
#define fftw_one HPCC_fftw_one

#endif



int HPCC_ipow(int x, int p);
/*
extern int HPCC_zfft1d(int n, fftw_complex *a, fftw_complex *b, int iopt, hpcc_fftw_plan p);
extern int HPCC_fft235(fftw_complex *a, fftw_complex *b, fftw_complex *w, int n, const int *ip);
extern int HPCC_settbl(fftw_complex *w, int n);
*/

// int HPCC_factor235(int n, int *ip);
// int HPCC_factor235_8(s64Int_t n, int *ip);

int HPCC_bcnrand(u64Int_t n, u64Int_t a, void *x);

#define ARR2D(a, i, j, lda) a[(i)+(j)*(lda)]
#define PTR2D(a, i, j, lda) (a+(i)+(j)*(lda))
#define ARR3D(a, i, j, k, lda1, lda2) a[(i)+(lda1)*((j)+(k)*(lda2))]
#define PTR3D(a, i, j, k, lda1, lda2) (a+(i)+(lda1)*((j)+(k)*(lda2)))
#define ARR4D(a, i, j, k, l, lda1, lda2, lda3) a[(i)+(lda1)*((j)+(lda2)*((k)+(lda3)*(l)))]
#define c_mul3v(v, v1, v2) c_re(v) = c_re(v1) * c_re(v2) - c_im(v1) * c_im(v2); c_im(v) = c_re(v1) * c_im(v2) + c_im(v1) * c_re(v2)
#define c_assgn(d, s) c_re(d)=c_re(s);c_im(d)=c_im(s)
#define V3MIN(r, e, v) r = (e); V2MIN(r, v)
#define V2MIN(r, v) r = (v) < r ? (v) : r
#define EMAX(d, v, e) d=(e); d=d>(v)?d:(v)

#define    Mmax( a_, b_ )      ( ( (a_) > (b_) ) ?  (a_) : (b_) )

#ifdef USING_FFTW
#include <fftw_mpi.h>
#else
#include <mpi.h>


typedef struct hpcc_fftw_mpi_plan_struct *fftw_mpi_plan;
#define fftw_mpi_create_plan  HPCC_fftw_mpi_create_plan
#define fftw_mpi_destroy_plan HPCC_fftw_mpi_destroy_plan
#define fftw_mpi HPCC_fftw_mpi
#define fftw_mpi_local_sizes HPCC_fftw_mpi_local_sizes
#endif

struct hpcc_fftw_mpi_plan_struct {
  MPI_Comm comm;
  MPI_Datatype cmplx;
  fftw_complex *wx, *wy, *wz, *c, *work;
  s64Int_t n;
  int flags, c_size;
  fftw_direction dir;
  double *timings;
};
typedef struct hpcc_fftw_mpi_plan_struct *hpcc_fftw_mpi_plan;

hpcc_fftw_mpi_plan
HPCC_fftw_mpi_create_plan(MPI_Comm comm, s64Int_t n, fftw_direction dir, int flags);
void HPCC_fftw_mpi_destroy_plan(hpcc_fftw_mpi_plan plan);
extern void HPCC_fftw_mpi(hpcc_fftw_mpi_plan p, int n_fields, fftw_complex *local_data,
                     fftw_complex *work);
extern void HPCC_fftw_mpi_local_sizes(hpcc_fftw_mpi_plan p, s64Int_t *local_n,
              s64Int_t *local_start, s64Int_t *local_n_after_transform,
              s64Int_t *local_start_after_transform, s64Int_t *total_local_size);

int
HPCC_pzfft1d(s64Int_t n, fftw_complex *a, fftw_complex *b, fftw_complex *w, int me, int npu, int iopt,
  hpcc_fftw_mpi_plan p);


double *HPCC_fft_timings_forward, *HPCC_fft_timings_backward;

/** HPL_dlamch.c **/
#ifndef DBL_DIGITS
#define DBL_DIGITS 53
#endif

/*
 * ---------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------
 */
static void     HPL_dlamc1
STDC_ARGS(
(  int *,           int *,           int *,           int * ) );
static void     HPL_dlamc2
STDC_ARGS(
(  int *,           int *,           int *,           double *,
   int *,           double *,        int *,           double * ) );
static double   HPL_dlamc3
STDC_ARGS(
(  const double,    const double ) );
static void     HPL_dlamc4
STDC_ARGS(
(  int *,           const double,    const int ) );
static void     HPL_dlamc5
STDC_ARGS(
(  const int,       const int,       const int,       const int,
   int *,           double * ) );
static double   HPL_dipow
STDC_ARGS(
(  const double,    const int ) );

#ifdef HPL_STDC_HEADERS
double HPL_dlamch
(
   const HPL_T_MACH                 CMACH
)
#else
double HPL_dlamch
( CMACH )
   const HPL_T_MACH                 CMACH;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_dlamch determines  machine-specific  arithmetic constants such as
 * the relative machine precision  (eps),  the safe minimum (sfmin) such
 * that 1 / sfmin does not overflow, the base of the machine (base), the
 * precision (prec), the  number of (base) digits  in the  mantissa (t),
 * whether rounding occurs in addition (rnd=1.0 and 0.0 otherwise),  the
 * minimum exponent before  (gradual)  underflow (emin),  the  underflow
 * threshold (rmin) base**(emin-1), the largest exponent before overflow
 * (emax), the overflow threshold (rmax) (base**emax)*(1-eps).
 *
 * Notes
 * =====
 * 
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function dlamch.f  (version 2.0 -- 1992), that  was  itself
 * based on the function ENVRON  by Malcolm and incorporated suggestions
 * by Gentleman and Marovich. See                                       
 *  
 * Malcolm M. A.,  Algorithms  to  reveal  properties  of floating-point
 * arithmetic.,  Comms. of the ACM, 15, 949-951 (1972).                 
 *  
 * Gentleman W. M. and Marovich S. B.,  More  on algorithms  that reveal
 * properties of  floating point arithmetic units.,  Comms. of  the ACM,
 * 17, 276-277 (1974).
 * 
 * Arguments
 * =========
 *
 * CMACH   (local input)                 const HPL_T_MACH
 *         Specifies the value to be returned by HPL_dlamch             
 *            = HPL_MACH_EPS,   HPL_dlamch := eps (default)             
 *            = HPL_MACH_SFMIN, HPL_dlamch := sfmin                     
 *            = HPL_MACH_BASE,  HPL_dlamch := base                      
 *            = HPL_MACH_PREC,  HPL_dlamch := eps*base                  
 *            = HPL_MACH_MLEN,  HPL_dlamch := t                         
 *            = HPL_MACH_RND,   HPL_dlamch := rnd                       
 *            = HPL_MACH_EMIN,  HPL_dlamch := emin                      
 *            = HPL_MACH_RMIN,  HPL_dlamch := rmin                      
 *            = HPL_MACH_EMAX,  HPL_dlamch := emax                      
 *            = HPL_MACH_RMAX,  HPL_dlamch := rmax                      
 *          
 *         where                                                        
 *          
 *            eps   = relative machine precision,                       
 *            sfmin = safe minimum,                                     
 *            base  = base of the machine,                              
 *            prec  = eps*base,                                         
 *            t     = number of digits in the mantissa,                 
 *            rnd   = 1.0 if rounding occurs in addition,               
 *            emin  = minimum exponent before underflow,                
 *            rmin  = underflow threshold,                              
 *            emax  = largest exponent before overflow,                 
 *            rmax  = overflow threshold.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   /*static*/ double              eps, sfmin, base, t, rnd, emin, rmin, emax,
                              rmax, prec;
   double                     small;
   /*static*/ int                 first=0/*1*/;
   int                        beta=0, imax=0, imin=0, it=0, lrnd=0;
/* ..
 * .. Executable Statements ..
 */
   eps = DBL_EPSILON / FLT_RADIX;
   base = FLT_RADIX;
   prec = DBL_EPSILON;
   t = DBL_DIGITS;
   rnd = FLT_ROUNDS < 2 ? HPL_rone : HPL_rzero;
   emin = DBL_MIN_EXP;
   rmin = DBL_MIN;
   emax = DBL_MAX_EXP;
   rmax = DBL_MAX;

   sfmin = rmin;
   small = HPL_rone / rmax;
   if (small >= sfmin)
     sfmin = small * ( HPL_rone + eps );

   if( first != 0 )
   {
      first = 0;
      HPL_dlamc2( &beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax );
      base  = (double)(beta);  t     = (double)(it);
      if( lrnd != 0 )
      { rnd = HPL_rone;  eps = HPL_dipow( base, 1 - it ) / HPL_rtwo; }
      else
      { rnd = HPL_rzero; eps = HPL_dipow( base, 1 - it );            }
      prec  = eps * base;  emin  = (double)(imin); emax  = (double)(imax);
      sfmin = rmin;        small = HPL_rone / rmax;
/*
 * Use  SMALL  plus a bit,  to avoid the possibility of rounding causing
 * overflow when computing  1/sfmin.
 */
      if( small >= sfmin ) sfmin = small * ( HPL_rone + eps );
   }

   if( CMACH == HPL_MACH_EPS   ) return( eps   );
   if( CMACH == HPL_MACH_SFMIN ) return( sfmin );
   if( CMACH == HPL_MACH_BASE  ) return( base  );
   if( CMACH == HPL_MACH_PREC  ) return( prec  );
   if( CMACH == HPL_MACH_MLEN  ) return( t     );
   if( CMACH == HPL_MACH_RND   ) return( rnd   );
   if( CMACH == HPL_MACH_EMIN  ) return( emin  );
   if( CMACH == HPL_MACH_RMIN  ) return( rmin  );
   if( CMACH == HPL_MACH_EMAX  ) return( emax  );
   if( CMACH == HPL_MACH_RMAX  ) return( rmax  );

   return( eps );
/*
 * End of HPL_dlamch
 */
}

#ifdef HPL_STDC_HEADERS
static void HPL_dlamc1
(
   int                        * BETA,
   int                        * T,
   int                        * RND,
   int                        * IEEE1
)
#else
static void HPL_dlamc1
( BETA, T, RND, IEEE1 )
/*
 * .. Scalar Arguments ..
 */
   int                        * BETA, * IEEE1, * RND, * T;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_dlamc1  determines  the machine parameters given by BETA, T, RND,
 * and IEEE1.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function dlamc1.f  (version 2.0 -- 1992), that  was  itself
 * based on the function ENVRON  by Malcolm and incorporated suggestions
 * by Gentleman and Marovich. See
 *
 * Malcolm M. A.,  Algorithms  to  reveal  properties  of floating-point
 * arithmetic.,  Comms. of the ACM, 15, 949-951 (1972).
 *
 * Gentleman W. M. and Marovich S. B.,  More  on algorithms  that reveal
 * properties of  floating point arithmetic units.,  Comms. of  the ACM,
 * 17, 276-277 (1974).
 *
 * Arguments
 * =========
 *
 * BETA    (local output)              int *
 *         The base of the machine.
 *
 * T       (local output)              int *
 *         The number of ( BETA ) digits in the mantissa.
 *
 * RND     (local output)              int *
 *         Specifies whether proper rounding (RND=1) or chopping (RND=0)
 *         occurs in addition.  This may not be a  reliable guide to the
 *         way in which the machine performs its arithmetic.
 *
 * IEEE1   (local output)              int *
 *         Specifies  whether  rounding  appears  to be done in the IEEE
 *         `round to nearest' style (IEEE1=1), (IEEE1=0) otherwise.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   double                     a, b, c, f, one, qtr, savec, t1, t2;
   static int                 first=1, lbeta, lieee1, lrnd, lt;
/* ..
 * .. Executable Statements ..
 */
   if( first != 0 )
   {
      first = 0; one = HPL_rone;
/*
 * lbeta, lieee1, lt and lrnd are the local values of BETA, IEEE1, T and
 * RND. Throughout this routine we use the function HPL_dlamc3 to ensure
 * that relevant values are stored and not held in registers, or are not
 * affected by optimizers.
 *
 * Compute  a = 2.0**m  with the  smallest  positive integer m such that
 * fl( a + 1.0 ) == a.
 */
      a = HPL_rone; c = HPL_rone;
      do
      { a *= HPL_rtwo; c = HPL_dlamc3( a, one ); c = HPL_dlamc3( c, -a ); }
      while( c == HPL_rone );
/*
 * Now compute b = 2.0**m with the smallest positive integer m such that
 * fl( a + b ) > a.
 */
      b = HPL_rone; c = HPL_dlamc3( a, b );
      while( c == a ) { b *= HPL_rtwo; c = HPL_dlamc3( a, b ); }
/*
 * Now compute the base.  a and c  are  neighbouring floating point num-
 * bers in the interval ( BETA**T, BETA**( T + 1 ) ) and so their diffe-
 * rence is BETA.  Adding 0.25 to c is to ensure that it is truncated to
 * BETA and not (BETA-1).
 */
      qtr = one / 4.0; savec = c;
      c   = HPL_dlamc3( c, -a ); lbeta = (int)(c+qtr);
/*
 * Now  determine  whether  rounding or chopping occurs, by adding a bit
 * less than BETA/2 and a bit more than BETA/2 to a.
 */
      b = (double)(lbeta);
      f = HPL_dlamc3( b / HPL_rtwo, -b / 100.0 ); c = HPL_dlamc3( f, a );
      if( c == a ) { lrnd = 1; } else { lrnd = 0; }
      f = HPL_dlamc3( b / HPL_rtwo,  b / 100.0 ); c = HPL_dlamc3( f, a );
      if( ( lrnd != 0 ) && ( c == a ) ) lrnd = 0;
/*
 * Try  and decide whether rounding is done in the  IEEE  round to nea-
 * rest style.  b/2 is half a unit in the last place of the two numbers
 * a  and savec. Furthermore, a is even, i.e. has last bit zero, and sa-
 * vec is odd.  Thus adding b/2 to a should not change a, but adding b/2
 * to savec should change savec.
 */
      t1 = HPL_dlamc3( b / HPL_rtwo, a );
      t2 = HPL_dlamc3( b / HPL_rtwo, savec );
      if ( ( t1 == a ) && ( t2 > savec ) && ( lrnd != 0 ) ) lieee1 = 1;
      else                                                  lieee1 = 0;
/*
 * Now find the mantissa, T. It should be the integer part of log to the
 * base BETA of a, however it is safer to determine T by powering. So we
 * find T as the smallest positive integer for which fl( beta**t + 1.0 )
 * is equal to 1.0.
 */
      lt = 0; a = HPL_rone; c = HPL_rone;

      do
      {
         lt++; a *= (double)(lbeta);
         c = HPL_dlamc3( a, one ); c = HPL_dlamc3( c,  -a );
      } while( c == HPL_rone );
   }

   *BETA  = lbeta; *T = lt; *RND = lrnd; *IEEE1 = lieee1;
} 

#ifdef HPL_STDC_HEADERS
static void HPL_dlamc2
(
   int                        * BETA, 
   int                        * T,
   int                        * RND,
   double                     * EPS,
   int                        * EMIN,
   double                     * RMIN,
   int                        * EMAX,
   double                     * RMAX
)
#else
static void HPL_dlamc2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
/*
 * .. Scalar Arguments ..
 */
   int                        * BETA, * EMAX, * EMIN, * RND, * T;
   double                     * EPS, * RMAX, * RMIN;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_dlamc2  determines the machine  parameters specified in its argu-
 * ment list.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function  dlamc2.f (version 2.0 -- 1992), that  was  itself
 * based on a function PARANOIA  by  W. Kahan of the University of Cali-
 * fornia at Berkeley for the computation of the  relative machine epsi-
 * lon eps.
 *
 * Arguments
 * =========
 *
 * BETA    (local output)              int *
 *         The base of the machine.
 *
 * T       (local output)              int *
 *         The number of ( BETA ) digits in the mantissa.
 *
 * RND     (local output)              int *
 *         Specifies whether proper rounding (RND=1) or chopping (RND=0)
 *         occurs in addition. This may not be a reliable  guide to  the
 *         way in which the machine performs its arithmetic.
 *
 * EPS     (local output)              double *
 *         The smallest positive number such that fl( 1.0 - EPS ) < 1.0,
 *         where fl denotes the computed value.
 *
 * EMIN    (local output)              int *
 *         The minimum exponent before (gradual) underflow occurs.
 *
 * RMIN    (local output)              double *
 *         The smallest  normalized  number  for  the  machine, given by
 *         BASE**( EMIN - 1 ), where  BASE  is the floating  point value
 *         of BETA.
 *
 * EMAX    (local output)              int *
 *         The maximum exponent before overflow occurs.
 *
 * RMAX    (local output)              double *
 *         The  largest  positive  number  for  the  machine,  given  by
 *         BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating  point
 *         value of BETA.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   static double              leps, lrmax, lrmin;
   double                     a, b, c, half, one, rbase, sixth, small,
                              third, two, zero;
   static int                 first=1, iwarn=0, lbeta=0, lemax, lemin,
                              lt=0;
   int                        gnmin=0, gpmin=0, i, ieee, lieee1=0,
                              lrnd=0, ngnmin=0, ngpmin=0;
/* ..
 * .. Executable Statements ..
 */
   if( first != 0 )
   {
      first = 0; zero = HPL_rzero; one = HPL_rone; two = HPL_rtwo;
/*
 * lbeta, lt, lrnd, leps, lemin and lrmin are the local values of  BETA,
 * T, RND, EPS, EMIN and RMIN.
 *
 * Throughout this routine we use the function HPL_dlamc3 to ensure that
 * relevant values are stored and not held in registers,  or are not af-
 * fected by optimizers.
 *
 * HPL_dlamc1 returns the parameters  lbeta, lt, lrnd and lieee1.
 */
      HPL_dlamc1( &lbeta, &lt, &lrnd, &lieee1 );
/*
 * Start to find eps.
 */
      b = (double)(lbeta); a = HPL_dipow( b, -lt ); leps = a;
/*
 * Try some tricks to see whether or not this is the correct  EPS.
 */
      b     = two / 3.0; 
      half  = one / HPL_rtwo;
      sixth = HPL_dlamc3( b, -half );
      third = HPL_dlamc3( sixth, sixth );
      b     = HPL_dlamc3( third, -half );
      b     = HPL_dlamc3( b, sixth );
      b     = Mabs( b ); if( b < leps ) b = leps;

      leps = HPL_rone;

      while( ( leps > b ) && ( b > zero ) )
      {
         leps = b;
         c = HPL_dlamc3( half * leps,
                         HPL_dipow( two, 5 ) * HPL_dipow( leps, 2 ) );
         c = HPL_dlamc3( half, -c ); b = HPL_dlamc3( half, c );
         c = HPL_dlamc3( half, -b ); b = HPL_dlamc3( half, c );
      }
      if( a < leps ) leps = a;
/*
 * Computation of EPS complete.
 *
 * Now find  EMIN.  Let a = + or - 1, and + or - (1 + BASE**(-3)).  Keep
 * dividing a by BETA until (gradual) underflow occurs. This is detected
 * when we cannot recover the previous a.
 */
      rbase = one / (double)(lbeta); small = one;
      for( i = 0; i < 3; i++ ) small = HPL_dlamc3( small * rbase, zero );
      a = HPL_dlamc3( one, small );
      HPL_dlamc4( &ngpmin, one, lbeta ); HPL_dlamc4( &ngnmin, -one, lbeta );
      HPL_dlamc4( &gpmin,    a, lbeta ); HPL_dlamc4( &gnmin,    -a, lbeta );

      ieee = 0;

      if( ( ngpmin == ngnmin ) && ( gpmin == gnmin ) )
      {
         if( ngpmin == gpmin )
         {
/*
 * Non twos-complement machines, no gradual underflow; e.g.,  VAX )
 */
            lemin = ngpmin;
         }
         else if( ( gpmin-ngpmin ) == 3 )
         {
/*
 * Non twos-complement machines with gradual underflow; e.g., IEEE stan-
 * dard followers
 */
            lemin = ngpmin - 1 + lt; ieee = 1;
         }
         else
         {
/*
 * A guess; no known machine
 */
            lemin = Mmin( ngpmin, gpmin );
            iwarn = 1;
         }
      }
      else if( ( ngpmin == gpmin ) && ( ngnmin == gnmin ) )
      {
         if( Mabs( ngpmin-ngnmin ) == 1 )
         {
/*
 * Twos-complement machines, no gradual underflow; e.g., CYBER 205
 */
            lemin = Mmax( ngpmin, ngnmin );
         }
         else
         {
/*
 * A guess; no known machine
 */
            lemin = Mmin( ngpmin, ngnmin );
            iwarn = 1;
         }
      }
      else if( ( Mabs( ngpmin-ngnmin ) == 1 ) && ( gpmin == gnmin ) )
      {
         if( ( gpmin - Mmin( ngpmin, ngnmin ) ) == 3 )
         {
/*
 * Twos-complement machines with gradual underflow; no known machine
 */
            lemin = Mmax( ngpmin, ngnmin ) - 1 + lt;
         }
         else
         {
/*
 * A guess; no known machine
 */
            lemin = Mmin( ngpmin, ngnmin );
            iwarn = 1;
         }
      }
      else
      {
/*
 * A guess; no known machine
 */
         lemin = Mmin( ngpmin, ngnmin ); lemin = Mmin( lemin, gpmin );
         lemin = Mmin( lemin, gnmin ); iwarn = 1;
      }
/*
 * Comment out this if block if EMIN is ok
 */
      if( iwarn != 0 )
      {
         first = 1;
         HPL_fprintf( stderr, "\n %s %8d\n%s\n%s\n%s\n",
"WARNING. The value EMIN may be incorrect:- EMIN =", lemin,
"If, after inspection, the value EMIN looks acceptable, please comment ",
"out the  if  block  as marked within the code of routine  HPL_dlamc2, ",
"otherwise supply EMIN explicitly." );
      }
/*
 * Assume IEEE arithmetic if we found denormalised  numbers above, or if
 * arithmetic seems to round in the  IEEE style,  determined  in routine
 * HPL_dlamc1.  A true  IEEE  machine should have both things true; how-
 * ever, faulty machines may have one or the other.
 */
      if( ( ieee != 0 ) || ( lieee1 != 0 ) ) ieee = 1;
      else                                   ieee = 0;
/*
 * Compute  RMIN by successive division by  BETA. We could compute  RMIN
 * as BASE**( EMIN - 1 ), but some machines underflow during this compu-
 * tation.
 */
      lrmin = HPL_rone;
      for( i = 0; i < 1 - lemin; i++ )
         lrmin = HPL_dlamc3( lrmin*rbase, zero );
/*
 * Finally, call HPL_dlamc5 to compute emax and rmax.
 */
      HPL_dlamc5( lbeta, lt, lemin, ieee, &lemax, &lrmax );
   }
   *BETA = lbeta; *T    = lt;    *RND  = lrnd;  *EPS  = leps;
   *EMIN = lemin; *RMIN = lrmin; *EMAX = lemax; *RMAX = lrmax;
} 

#ifdef HPL_STDC_HEADERS
static double HPL_dlamc3( const double A, const double B )
#else
static double HPL_dlamc3( A, B )
/*
 * .. Scalar Arguments ..
 */
   const double               A, B;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_dlamc3  is intended to force a and b  to be stored prior to doing
 * the addition of  a  and  b,  for  use  in situations where optimizers
 * might hold one of these in a register.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function dlamc3.f (version 2.0 -- 1992).
 *
 * Arguments
 * =========
 *
 * A, B    (local input)               double
 *         The values a and b.
 *
 * ---------------------------------------------------------------------
 */
/* ..
 * .. Executable Statements ..
 */
   return( A + B );
} 

#ifdef HPL_STDC_HEADERS
static void HPL_dlamc4
(
   int                        * EMIN,
   const double               START,
   const int                  BASE
)
#else
static void HPL_dlamc4( EMIN, START, BASE )
/*
 * .. Scalar Arguments ..
 */
   int                        * EMIN;
   const int                  BASE;
   const double               START;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_dlamc4 is a service function for HPL_dlamc2.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function dlamc4.f (version 2.0 -- 1992).
 *
 * Arguments
 * =========
 *
 * EMIN    (local output)              int *
 *         The minimum exponent before  (gradual) underflow, computed by
 *         setting A = START and dividing  by  BASE until the previous A
 *         can not be recovered.
 *
 * START   (local input)               double
 *         The starting point for determining EMIN.
 *
 * BASE    (local input)               int
 *         The base of the machine.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   double                     a, b1, b2, c1, c2, d1, d2, one, rbase, zero;
   int                        i;
/* ..
 * .. Executable Statements ..
 */
   a     = START; one = HPL_rone; rbase = one / (double)(BASE);
   zero  = HPL_rzero;
   *EMIN = 1; b1 = HPL_dlamc3( a * rbase, zero ); c1 = c2 = d1 = d2 = a;

   do
   {
      (*EMIN)--; a = b1;
      b1 = HPL_dlamc3( a /  BASE,  zero );
      c1 = HPL_dlamc3( b1 *  BASE, zero );
      d1 = zero; for( i = 0; i < BASE; i++ ) d1 = d1 + b1;
      b2 = HPL_dlamc3( a * rbase,  zero );
      c2 = HPL_dlamc3( b2 / rbase, zero );
      d2 = zero; for( i = 0; i < BASE; i++ ) d2 = d2 + b2;
   } while( ( c1 == a ) && ( c2 == a ) &&  ( d1 == a ) && ( d2 == a ) );
} 

#ifdef HPL_STDC_HEADERS
static void HPL_dlamc5
(
   const int                  BETA,
   const int                  P, 
   const int                  EMIN,
   const int                  IEEE,
   int                        * EMAX,
   double                     * RMAX
)
#else
static void HPL_dlamc5( BETA, P, EMIN, IEEE, EMAX, RMAX )
/*
 * .. Scalar Arguments ..
 */
   const int                  BETA, EMIN, IEEE, P; 
   int                        * EMAX;
   double                     * RMAX;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_dlamc5  attempts  to compute RMAX, the largest machine  floating-
 * point number, without overflow.  It assumes that EMAX + abs(EMIN) sum
 * approximately to a power of 2.  It will fail  on machines where  this
 * assumption does not hold, for example, the  Cyber 205 (EMIN = -28625,
 * EMAX = 28718).  It will also fail if  the value supplied for  EMIN is
 * too large (i.e. too close to zero), probably with overflow.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function dlamc5.f (version 2.0 -- 1992).
 *
 * Arguments
 * =========
 *
 * BETA    (local input)               int
 *         The base of floating-point arithmetic.
 *
 * P       (local input)               int
 *         The number of base BETA digits in the mantissa of a floating-
 *         point value.
 *
 * EMIN    (local input)               int
 *         The minimum exponent before (gradual) underflow.
 *
 * IEEE    (local input)               int
 *         A logical flag specifying whether or not  the arithmetic sys-
 *         tem is thought to comply with the IEEE standard.
 *
 * EMAX    (local output)              int *
 *         The largest exponent before overflow.
 *
 * RMAX    (local output)              double *
 *         The largest machine floating-point number.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   double                     oldy=HPL_rzero, recbas, y, z;
   int                        exbits=1, expsum, i, lexp=1, nbits, try_,
                              uexp;
/* ..
 * .. Executable Statements ..
 */
/*
 * First compute  lexp  and  uexp, two powers of 2 that bound abs(EMIN).
 * We then assume that  EMAX + abs( EMIN ) will sum approximately to the
 * bound that  is closest to abs( EMIN ). (EMAX  is the  exponent of the
 * required number RMAX).
 */
l_10:
   try_ = (int)( (unsigned int)(lexp) << 1 );
   if( try_ <= ( -EMIN ) ) { lexp = try_; exbits++; goto l_10; }

   if( lexp == -EMIN ) { uexp = lexp; } else { uexp = try_; exbits++; }
/*
 * Now -lexp is less than or equal to EMIN, and -uexp is greater than or
 * equal to EMIN. exbits is the number of bits needed to store the expo-
 * nent.
 */
   if( ( uexp+EMIN ) > ( -lexp-EMIN ) )
   { expsum = (int)( (unsigned int)(lexp) << 1 ); }
   else
   { expsum = (int)( (unsigned int)(uexp) << 1 ); }
/*
 * expsum is the exponent range, approximately equal to EMAX - EMIN + 1.
 */
   *EMAX = expsum + EMIN - 1;
/*
 * nbits  is  the total number of bits needed to store a  floating-point
 * number.
 */
   nbits = 1 + exbits + P;

   if( ( nbits % 2 == 1 ) && ( BETA == 2 ) )
   {
/*
 * Either there are an odd number of bits used to store a floating-point
 * number, which is unlikely, or some bits are not used in the represen-
 * tation of numbers,  which is possible,  (e.g. Cray machines)  or  the
 * mantissa has an implicit bit, (e.g. IEEE machines, Dec Vax machines),
 * which is perhaps the most likely. We have to assume the last alterna-
 * tive.  If this is true,  then we need to reduce  EMAX  by one because
 * there must be some way of representing zero  in an  implicit-bit sys-
 * tem. On machines like Cray we are reducing EMAX by one unnecessarily.
 */
      (*EMAX)--;
   }

   if( IEEE != 0 )
   {
/*
 * Assume we are on an IEEE  machine which reserves one exponent for in-
 * finity and NaN.
 */
      (*EMAX)--;
   }
/*
 * Now create RMAX, the largest machine number, which should be equal to
 * (1.0 - BETA**(-P)) * BETA**EMAX . First compute 1.0-BETA**(-P), being
 * careful that the result is less than 1.0.
 */
   recbas = HPL_rone / (double)(BETA);
   z      = (double)(BETA) - HPL_rone;
   y      = HPL_rzero;

   for( i = 0; i < P; i++ )
   { z *= recbas; if( y < HPL_rone ) oldy = y; y = HPL_dlamc3( y, z ); }

   if( y >= HPL_rone ) y = oldy;
/*
 * Now multiply by BETA**EMAX to get RMAX.
 */
   for( i = 0; i < *EMAX; i++ ) y = HPL_dlamc3( y * BETA, HPL_rzero );

   *RMAX = y;
/*
 * End of HPL_dlamch
 */
} 

#ifdef HPL_STDC_HEADERS
static double HPL_dipow
(
   const double               X,
   const int                  N
)
#else
static double HPL_dipow( X, N )
/*
 * .. Scalar Arguments ..
 */
   const int                  N;
   const double               X;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_dipow computes the integer n-th power of a real scalar x.
 *
 * Arguments
 * =========
 *
 * X       (local input)               const double
 *         The real scalar x.
 *
 * N       (local input)               const int
 *         The integer power to raise x to.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   double                     r, y=HPL_rone;
   int                        k, n;
/* ..
 * .. Executable Statements ..
 */
   if( X == HPL_rzero ) return( HPL_rzero );
   if( N < 0 ) { n = -N; r = HPL_rone / X; } else { n = N; r = X; }
   for( k = 0; k < n; k++ ) y *= r; 

   return( y );
}

/*
 * End of HPL_dlamch.c
 */

/** HPL_fprintf.c **/
#ifdef HPL_STDC_HEADERS
void HPL_fprintf
(
   FILE *                           STREAM,
   const char *                     FORM,
   ...                              
)
#else
void HPL_fprintf( va_alist )
va_dcl
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_fprintf is a wrapper around fprintf flushing the output stream.
 * 
 *
 * Arguments
 * =========
 *
 * STREAM  (local input)                 FILE *
 *         On entry, STREAM specifies the output stream.
 *
 * FORM    (local input)                 const char *
 *         On entry, FORM specifies the format, i.e., how the subsequent
 *         arguments are converted for output.
 *
 *         (local input)                 ...
 *         On entry,  ...  is the list of arguments to be printed within
 *         the format string.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   va_list                    argptr;
   char                       cline[256];
#ifndef HPL_STDC_HEADERS
   FILE                       * STREAM;
   char                       * FORM;
#endif
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_STDC_HEADERS
   va_start( argptr, FORM );
#else
   va_start( argptr );
   STREAM = va_arg( argptr, FILE * );
   FORM   = va_arg( argptr, char * );
#endif
   (void) vsprintf( cline, FORM, argptr );
   va_end( argptr ); 

   (void) fprintf( STREAM, "%s", cline );
   (void) fflush( STREAM );
/*
 * End of HPL_fprintf
 */
} 


/*
 * Begin of bcnrand.c
 */
typedef u64Int_t Big[2];

/* r = a * b */
static void
ddmuldd(u64Int_t a, u64Int_t b, Big r) {
  u64Int_t a0, a1, b0, b1, hb, acc, acc1;

  /* 'hb' should be 0xFFFFFFFF (first 32-bits set to one) */
  hb = 65535L;
  hb = (hb << 16) | hb;

  /* split 'a' and 'b' into two 32-bit quantities */
  a0 = a & hb;
  a1 = (a >> 32) & hb;
  b0 = b & hb;
  b1 = (b >> 32) & hb;

  acc = a0 * b0;
  r[0] = acc & hb;
  acc >>= 32;
  acc += a1 * b0;
  acc1 = acc >> 32;
  acc &= hb;
  acc += a0 * b1;
  r[0] += (acc & hb) << 32;
  acc >>= 32;
  acc += acc1;
  acc += a1 * b1;
  r[1] = acc;
}

/* r = a - b */
static void
ddsub(Big a, Big b, Big r) {
  u64Int_t mx = 0;
  mx = ~mx;

  r[1] = a[1] - b[1];
  if (a[0] >= b[0])
    r[0] = a[0] - b[0];
  else {
    r[1] -= 1;
    r[0] = mx - b[0] + 1 + a[0];
  }
}

/* q = d / v; reminder Ur */
static void
dddiv(Big d, u64Int_t v, Big q, u64Int_t *Ur) {
  u64Int_t r1, r0, v1, v0, msb = 1, mx = 0, one = 1;
  int i;

  msb <<= 63;
  mx = ~mx;
  q[0] = q[1] = 0;

  if (v <= d[1]) {
    q[1] = d[1] / v;
    r1 = d[1] % v;
  } else {
    r1 = d[1];
    q[1] = 0;
  }
  r0 = d[0];

  while (r1) {
    v1 = 0;
    v0 = v;

    for (i = 0; v1 <= r1; i++) {
      v1 <<= 1;
      if (msb & v0) v1 |= 1;
      v0 <<= 1;
    }
    do {
    i--;
    v0 >>= 1;
    v0 &= mx;
    if (1 & v1) v0 |= msb;
    v1 >>= 1;
    } while (v1 == r1 && v0 > r0); /* make sure (v1,v0) is not too big */

    q[0] += one << i;
    r1 -= v1;
    if (r0 >= v0)
      r0 -= v0;
    else {
      r0 += mx - v0 + 1;
      r1 -= 1;
    }

  }

  q[0] += r0 / v;
  r0 %= v;

  if (Ur) *Ur = r0;
}

/*
!   expm2 = 2^p mod am.  p2 is a table with powers of 2, i.e., p2(i) = 2^i.
!   This routine uses a left-to-right binary exponentiation scheme.
*/
static u64Int_t
expm2(u64Int_t p, u64Int_t am) {
  u64Int_t p2, p1, pt1, r;
  Big ddm, dd1, dd2;
  int i;

  for (p2 = i = 1; i < 54; i++) {
    p2 <<= 1;
    if (p2 > p) break;
  }

  p1 = p;
  pt1 = p2 >> 1;
  r = 1;
  ddm[0] = am;
  ddm[1] = 0;

  while (1) {
    if (p1 >= pt1) {
      /* r = mod(2.0 * r, am) */
      ddmuldd( 2, r, dd1 );
      if (dd1[0] > am) {
	ddsub(dd1, ddm, dd2);
	dd1[0] = dd2[0];
	dd1[1] = dd2[1];
      }
      r = dd1[0];
      p1 = p1 - pt1;
    }

    pt1 /= 2;
    if (pt1 >= 1) {
      /* r = mod(r * r, am) */
      ddmuldd( r, r, dd1 );
      dddiv( dd1, am, dd2, &r );
      continue;
    }
    break;
  }

  return r;
}

/*
  Let minA = 3^33 + 100
  If `a' is smaller than `minA' then `a' is incremented by `minA' this value.
  In this way, you can seed the generator with small integers and the requirements
  will be fullfilled internally.
 */
int
HPCC_bcnrand(u64Int_t n, u64Int_t a, void *vx) {
  u64Int_t d1, d2, d3, t53, p3i, ui, minA;
  s64Int_t sd1, sp3i;
  Big dd1, dd2;
  int i;
  double rcp, two64, v, *x = (double *)vx;

  /* minA = 3.d0 ** 33 + 100.d0 */
  minA = 20709114;
  minA <<= 28;
  minA += 106609639;

  /* make sure `a' is big enough */
  if (a < minA) a += minA;

  t53 = 1;
  t53 <<= 53;

  d1 = 1;
  for (i = 0; i < 53; i++) {
    d1 *= 3;
    if (d1 > a) break;
  }

  /* two64 = 2 ** 64 */
  two64 = 2.0;
  for (i = 0; i < 6; i++)
    two64 *= two64;

  p3i = d1 / 3;
  sp3i = (s64Int_t)p3i;
  rcp = 1.0 / p3i;

  /*
!   Calculate starting element.  This code performs the following:
!   d1 = [int[p3i/2] * 2^(a-p3i)] mod p3i.
  */
  /* d1 = (p3i/2 * (2 ** (a-p3i))) % p3i */
  d2 = expm2( a - p3i, p3i );
  d3 = p3i / 2;
  ddmuldd( d2, d3, dd1 );
  dddiv( dd1, p3i, dd2, &d1 );

  x[0] = d1 * rcp;
  for (ui = 1; ui < n; ui++) {
    /* dd1 = d1 * t53 */
    dd1[1] = (d1 >> 11);
    dd1[0] = (d1 << 53);

    /* Approximate `dd1/p3i' (the result should be off by 1) */
    v = ((two64 * (double)dd1[1]) + (double)dd1[0]) * rcp;

    ddmuldd( (u64Int_t)v, p3i, dd2 );

    /* The value of `dd1-dd2' should between `-p3i' and 'p3i',
     hence upper halves of `dd1' and `dd2' can be ignored */
    sd1 = (s64Int_t)(dd1[0] - dd2[0]);

    /* Check the above approximation */
    if (sd1 < 0) sd1 += sp3i;
    if (sd1 > sp3i) sd1 -= sp3i;

    /* d1 = (d1 * t53) % p3i */
    d1 = (u64Int_t)sd1;
    x[ui] = d1 * rcp;
  }

  return 0;
}
/*
 * End of bcnrand.c
 */

/*
 * Begin of wrapmpifftw.c
 */

#define    Mmax3( a_, b_, c_ )      ( (a_) > (b_) ?  ((a_) > (c_) ? (a_) : (c_)) : ((b_) > (c_) ? (b_) : (c_)) )

static int
GetNXYZ(s64Int_t n, int npu) {
  int ip[3], lnx[3], lny[3], lnz[3], lnpu[3];
  int i, nx, ny, nz, nxyz;

  HPCC_factor235( npu, lnpu );
  HPCC_factor235_8( n, ip );

  for (i = 0; i < 3; ++i) {
    EMAX( lnz[i], lnpu[i], (ip[i]+1)/3 );
    EMAX( lnx[i], lnpu[i], (ip[i]-lnz[i]+1)/2 );
    lny[i] = ip[i] - lnx[i] - lnz[i];
  }

  nx = HPCC_ipow( 2, lnx[0] ) * HPCC_ipow( 3, lnx[1] ) * HPCC_ipow( 5, lnx[2] );
  ny = HPCC_ipow( 2, lny[0] ) * HPCC_ipow( 3, lny[1] ) * HPCC_ipow( 5, lny[2] );
  nz = HPCC_ipow( 2, lnz[0] ) * HPCC_ipow( 3, lnz[1] ) * HPCC_ipow( 5, lnz[2] );

  nxyz = Mmax3( nx, ny, nz );

  return nxyz;
}

hpcc_fftw_mpi_plan
HPCC_fftw_mpi_create_plan(MPI_Comm comm, s64Int_t n, fftw_direction dir, int flags) {
  hpcc_fftw_mpi_plan p;
  fftw_complex *a = NULL, *b = NULL;
  int nxyz;
  int rank, size;

  MPI_Comm_size( comm, &size );
  MPI_Comm_rank( comm, &rank );

  p = (hpcc_fftw_mpi_plan)fftw_malloc( sizeof *p );
  if (! p) return p;

  nxyz = GetNXYZ( n, size );

  p->wx = (fftw_complex *)HPCC_fftw_malloc( (nxyz/2 + FFTE_NP) * (sizeof *p->wx) );
  p->wy = (fftw_complex *)HPCC_fftw_malloc( (nxyz/2 + FFTE_NP) * (sizeof *p->wy) );
  p->wz = (fftw_complex *)HPCC_fftw_malloc( (nxyz/2 + FFTE_NP) * (sizeof *p->wz) );
  p->work = (fftw_complex *)HPCC_fftw_malloc( n / size * 3 / 2 * (sizeof *p->work) );

  p->c_size = (nxyz+FFTE_NP) * (FFTE_NBLK + 1) + FFTE_NP;
#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single
    {
      int i;
      i = omp_get_num_threads();
      p->c = (fftw_complex *)HPCC_fftw_malloc( p->c_size * (sizeof *p->c) * i );
    }
  }
#else
  p->c = (fftw_complex *)HPCC_fftw_malloc( p->c_size * (sizeof *p->c) );
#endif

  if (! p->wx || ! p->wy || ! p->wz || ! p->work || ! p->c) {
    if (p->c) HPCC_fftw_free( p->c );
    if (p->work) HPCC_fftw_free( p->work );
    if (p->wz) HPCC_fftw_free( p->wz );
    if (p->wy) HPCC_fftw_free( p->wy );
    if (p->wx) HPCC_fftw_free( p->wx );
    fftw_free( p );
    return NULL;
  }

  p->n = n;
  p->comm = comm;
  p->dir = dir;
  p->flags = flags;

  MPI_Type_contiguous( 2, MPI_DOUBLE, &p->cmplx );
  MPI_Type_commit( &p->cmplx );

  if (FFTW_FORWARD == p->dir)
    p->timings = HPCC_fft_timings_forward;
  else
    p->timings = HPCC_fft_timings_backward;

  HPCC_pzfft1d( n, a, b, p->work, rank, size, 0, p );

  return p;
}

void
HPCC_fftw_mpi_destroy_plan(hpcc_fftw_mpi_plan p) {
  if (!p) return;

  MPI_Type_free( &p->cmplx );

  HPCC_fftw_free( p->work );
  HPCC_fftw_free( p->c );
  HPCC_fftw_free( p->wz );
  HPCC_fftw_free( p->wy );
  HPCC_fftw_free( p->wx );
  fftw_free( p );
}

void
HPCC_fftw_mpi(hpcc_fftw_mpi_plan p, int n_fields, fftw_complex *local_data, fftw_complex *work){
  int rank, size;
  s64Int_t n;
  int i, ln;

  MPI_Comm_size( p->comm, &size );
  MPI_Comm_rank( p->comm, &rank );

  n = p->n;

  if (FFTW_FORWARD == p->dir)
    HPCC_pzfft1d( n, local_data, work, p->work, rank, size, -1, p );
  else
    HPCC_pzfft1d( n, local_data, work, p->work, rank, size, +1, p );

  ln = n / size;
  for (i = 0; i < ln; ++i) {
    c_assgn( local_data[i], work[i] );
  }
}

void
HPCC_fftw_mpi_local_sizes(hpcc_fftw_mpi_plan p, s64Int_t *local_n, s64Int_t *local_start,
  s64Int_t *local_n_after_transform, s64Int_t *local_start_after_transform, s64Int_t *total_local_size) {
  int rank, size;
  s64Int_t n;
  MPI_Comm_size( p->comm, &size );
  MPI_Comm_rank( p->comm, &rank );
  n = p->n;
  if (local_n) *local_n = n / size;
  if (local_start) *local_start = n / size * rank;
  if (local_n_after_transform) *local_n_after_transform = n / size;
  if (local_start_after_transform) *local_start_after_transform = n / size * rank;
  if (total_local_size) *total_local_size = n / size;
}

/*
 * End of wrapmpifftw.c
 */

/*
 * Begin of pzfft1d.c
 */

static void
ztrans(fftw_complex *a, fftw_complex *b, int n1, int n2) {
  int ii, jj, i, j;
  int tmin1, tmin2;
  int lda, ldb;

  lda = n1;
  ldb = n2;

#ifdef _OPENMP
#pragma omp parallel for private(i,j,jj,tmin1,tmin2)
#endif
  for (ii = 0; ii < n1; ii += FFTE_NBLK)
    for (jj = 0; jj < n2; jj += FFTE_NBLK) {

      V3MIN( tmin1, ii + FFTE_NBLK, n1 );
      for (i = ii; i < tmin1; ++i) {

        V3MIN( tmin2, jj + FFTE_NBLK, n2 );
        for (j = jj; j < tmin2; ++j) {
          c_assgn( ARR2D( b, j, i, ldb ), ARR2D( a, i, j, lda ) );
        }
      }
    }
}	/* ztrans */

static void
pztrans(fftw_complex *a, fftw_complex *b, s64Int_t nn, hpcc_fftw_mpi_plan p, int npu) {
  int i, nn2;

  nn2 = nn / npu;

  if (1 == npu)
    for (i = 0; i < nn2; ++i) {
      c_assgn( b[i], a[i] );
    }
  else
    MPI_Alltoall( a, nn2, p->cmplx, b, nn2, p->cmplx, p->comm );
}	/* pztrans */

static void
pzfft1d0(fftw_complex *a, fftw_complex *a2, fftw_complex *apxyz, fftw_complex *axyzp, fftw_complex *b,
  fftw_complex *bxyzp, fftw_complex *bzyx, fftw_complex *cy, fftw_complex *cz, fftw_complex *d,
  fftw_complex *wx, fftw_complex *wy, fftw_complex *wz, fftw_complex *ww, fftw_complex *www,
  int nx, int ny, int nz, hpcc_fftw_mpi_plan p, int npu, const int *lnx, const int *lny, const int *lnz) {

  int i, j, k, l, ii, jj, kk;
  int tmin1, tmin2, tmin3;
  int nnx, nnz;
  s64Int_t nn;
  int ldcz, lda2_1, lda2_2, ldaxyzp1, ldaxyzp2, ldaxyzp3, ldbxyzp1, ldbxyzp2, ldbxyzp3, ldww, ldcy;
  int ldwww1, ldwww2, ldwww3, ldapxyz1, ldapxyz2, ldapxyz3, ldbzyx1, ldbzyx2, lda1, lda2;
  fftw_complex ztmp1, ztmp2, ztmp3;

  ldcz = nz + FFTE_NP;
  lda2_1 = nx / npu;
  lda2_2 = ny;
  ldaxyzp1 = nx / npu;
  ldaxyzp2 = ny;
  ldaxyzp3 = nz / npu;
  ldbxyzp1 = nx / npu;
  ldbxyzp2 = ny;
  ldbxyzp3 = nz / npu;
  ldww = ny;
  ldcy = ny + FFTE_NP;
  ldwww1 = npu;
  ldwww2 = nx / npu;
  ldwww3 = ny;
  ldapxyz1 = npu;
  ldapxyz2 = nx / npu;
  ldapxyz3 = ny;
  ldbzyx1 = nz / npu;
  ldbzyx2 = ny;
  lda1 = nx;
  lda2 = ny;

  nnx = nx / npu;
  nnz = nz / npu;
  nn = (s64Int_t)nx * ny * nz / npu;

#ifdef _OPENMP
#pragma omp for private(i,k,l,ii,kk,tmin1,tmin2)
#endif
  for (j = 0; j < ny; ++j) {
    for (ii = 0; ii < nnx; ii += FFTE_NBLK) {
      for (kk = 0; kk < nz; kk += FFTE_NBLK) {

        V3MIN( tmin1, ii + FFTE_NBLK, nnx );
        for (i = ii; i < tmin1; ++i) {

          V3MIN( tmin2, kk + FFTE_NBLK, nz );
          for (k = kk; k < tmin2; ++k) {
            c_assgn( ARR2D( cz, k, i-ii, ldcz ), ARR3D( a2, i, j, k, lda2_1, lda2_2 ) );
          }
        }
      }

      V3MIN( tmin2, ii + FFTE_NBLK, nnx );
      for (i = ii; i < tmin2; ++i)
        HPCC_fft235( PTR2D( cz, 0, i-ii, ldcz ), d, wz, nz, lnz );

      for (l = 0; l < npu; ++l) {
        for (k = 0; k < nnz; ++k) {

          /* reusing tmin2 from above */
          for (i = ii; i < tmin2; ++i) {
            c_assgn( ARR4D( axyzp, i, j, k, l, ldaxyzp1, ldaxyzp2, ldaxyzp3 ),
                     ARR2D( cz, l + k*npu, i-ii, ldcz ) );
          }
        }
      }
    }
  }

#ifdef _OPENMP
#pragma omp single
  {
#endif
  p->timings[3] = MPI_Wtime();

  pztrans( a, b, nn, p, npu );

  p->timings[4] = MPI_Wtime();
#ifdef _OPENMP
  }
#endif

#ifdef _OPENMP
#pragma omp for private(i,j,l,ii,jj,kk,tmin1,tmin2)
#endif
  for (k = 0; k < nnz; ++k) {
    for (l = 0; l < npu; ++l) {
      for (ii = 0; ii < nnx; ii += FFTE_NBLK) {
        for (jj = 0; jj < ny; jj += FFTE_NBLK) {

          V3MIN( tmin1, ii + FFTE_NBLK, nnx );
          for (i = ii; i < tmin1; ++i) {

            V3MIN( tmin2, jj + FFTE_NBLK, ny );
            for (j = jj; j < tmin2; ++j) {
              c_assgn( ztmp1, ARR4D( bxyzp, i, j, k, l, ldbxyzp1, ldbxyzp2, ldbxyzp3 ) );
              c_assgn( ztmp2, ARR2D( ww, j, k, ldww ) );
              c_mul3v(ztmp3, ztmp1, ztmp2);
              c_assgn( ARR2D( cy, j, i-ii, ldcy ), ztmp3 );
            }
          }
        }

        V3MIN( tmin1, ii + FFTE_NBLK, nnx );
        for (i = ii; i < tmin1; ++i)
          HPCC_fft235( PTR2D( cy, 0, i-ii, ldcy ), d, wy, ny, lny );

        for (j = 0; j < ny; ++j) {
        V3MIN( tmin1, ii + FFTE_NBLK, nnx );
          for (i = ii; i < tmin1; ++i) {
            c_assgn( ztmp1, ARR2D( cy, j, i-ii, ldcy ) );
            c_assgn( ztmp2, ARR4D( www, l, i, j, k, ldwww1, ldwww2, ldwww3 ) );
            c_mul3v(ztmp3, ztmp1, ztmp2);
            c_assgn( ARR4D( apxyz, l, i, j, k, ldapxyz1, ldapxyz2, ldapxyz3 ), ztmp3 );
          }
        }
      }
    }

    for (j = 0; j < ny; ++j)
      HPCC_fft235( PTR3D( a, 0, j, k, lda1, lda2 ), d, wx, nx, lnx );
  }

#ifdef _OPENMP
#pragma omp for private(i,j,k,jj,kk,tmin1,tmin2,tmin3)
#endif
  for (ii = 0; ii < nx; ii += FFTE_NBLK) {
    for (jj = 0; jj < ny; jj += FFTE_NBLK) {
      for (kk = 0; kk < nnz; kk += FFTE_NBLK) {

        V3MIN( tmin1, ii + FFTE_NBLK, nx );
        for (i = ii; i < tmin1; ++i) {

          V3MIN( tmin2, jj + FFTE_NBLK, ny );
          for (j = jj; j < tmin2; ++j) {

            V3MIN( tmin3, kk + FFTE_NBLK, nnz );
            for (k = kk; k < tmin3; ++k) {
              c_assgn( ARR3D( bzyx, k, j, i, ldbzyx1, ldbzyx2 ), ARR3D( a, i, j, k, lda1, lda2 ) );
            }
          }
        }
      }
    }
  }
}	/* pzfft1d0 */

static void
psettbl2(fftw_complex *w, int ny, int nz, int me, int npu) {
  int j, k;
  int ldw;
  double pi2, px;
  int tmin1;

  ldw = ny;

  pi2 = 8.0 * atan(1.0);
  px = -pi2 / ((double)ny * nz);

  tmin1 = nz / npu;
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
  for (k = 0; k < tmin1; ++k)
    for (j = 0; j < ny; ++j) {
      c_re( ARR2D( w, j, k, ldw ) ) = cos(px * j * (me + (double)k * npu));
      c_im( ARR2D( w, j, k, ldw ) ) = sin(px * j * (me + (double)k * npu));
    }
}	/* psettbl2 */

static void
psettbl3(fftw_complex *w, int nx, int ny, int nz, int me, int npu) {
  int i, j, k;
  int ldw1, ldw2;
  int tmin1;
  double pi2, px;

  ldw1 = nx;
  ldw2 = ny;

  pi2 = 8.0 * atan(1.0);
  px = -pi2 / ((double)nx * ny * nz);

  tmin1 = nz / npu;
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
  for (k = 0; k < tmin1; ++k)
    for (j = 0; j < ny; ++j)
      for (i = 0; i < nx; ++i) {
        c_re( ARR3D( w, i, j, k, ldw1, ldw2 ) ) = cos( px * i * (me + (double)k * npu + (double)j * nz));
        c_im( ARR3D( w, i, j, k, ldw1, ldw2 ) ) = sin( px * i * (me + (double)k * npu + (double)j * nz));
      }
}	/* psettbl3 */

int
HPCC_pzfft1d(s64Int_t n, fftw_complex *a, fftw_complex *b, fftw_complex *w, int me, int npu, int iopt,
  hpcc_fftw_mpi_plan p) {

  int ip[3], lnx[3], lny[3], lnz[3], lnpu[3];
  s64Int_t nn;
  int i, inn, nn2, nd, nx, ny, nz;
  fftw_complex *wx, *wy, *wz, *c;
  double dn;

  p->timings[0] = MPI_Wtime();

  wx = p->wx;
  wy = p->wy;
  wz = p->wz;
  c = p->c;

  nn = n / npu; inn = (int)nn;
  nn2 = nn / npu;

  HPCC_factor235( npu, lnpu );
  HPCC_factor235_8( n, ip );

  for (i = 0; i < 3; ++i) {
    EMAX( lnz[i], lnpu[i], (ip[i]+1)/3 );
    EMAX( lnx[i], lnpu[i], (ip[i]-lnz[i]+1)/2 );
    lny[i] = ip[i] - lnx[i] - lnz[i];
  }

  nx = HPCC_ipow( 2, lnx[0] ) * HPCC_ipow( 3, lnx[1] ) * HPCC_ipow( 5, lnx[2] );
  ny = HPCC_ipow( 2, lny[0] ) * HPCC_ipow( 3, lny[1] ) * HPCC_ipow( 5, lny[2] );
  nz = HPCC_ipow( 2, lnz[0] ) * HPCC_ipow( 3, lnz[1] ) * HPCC_ipow( 5, lnz[2] );

  if (0 == iopt) {
    HPCC_settbl( wx, nx );
    HPCC_settbl( wy, ny );
    HPCC_settbl( wz, nz );
    psettbl2( w, ny, nz, me, npu );
    psettbl3( w + ny * (nz / npu), nx, ny, nz, me, npu );
    return 0;
  }

  if (1 == iopt || 2 == iopt) {
    for (i = 0; i < inn; ++i) {
      c_im( a[i] ) = -c_im( a[i] );
    }
  }

  p->timings[1] = MPI_Wtime();

  if (-1 == iopt || 1 == iopt || -2 == iopt) {
    ztrans( a, b, npu, nn2 );
    pztrans( b, a, nn, p, npu );
  }

  p->timings[2] = MPI_Wtime();

  nd = ((ny > nz ? ny : nz) + FFTE_NP) * FFTE_NBLK + FFTE_NP;

#ifdef _OPENMP
#pragma omp parallel private(c,i)
   {
    i = omp_get_thread_num();
    c = p->c + i*p->c_size;
#endif

  pzfft1d0( a, a, a, a, b, b, b, c, c, c + nd, wx, wy, wz, w, w + ny*(nz/npu), nx, ny, nz, p, npu, lnx, lny, lnz );

#ifdef _OPENMP
   }
#endif

  p->timings[5] = MPI_Wtime();

  if (-1 == iopt || 1 == iopt || 2 == iopt) {
    pztrans( b, a, nn, p, npu );
    ztrans( a, b, nn2, npu );
  }

  p->timings[6] = MPI_Wtime();

  if (1 == iopt || 2 == iopt) {
    dn = 1.0 / n;
    for (i = 0; i < inn; ++i) {
      c_re( b[i] ) *= dn;
      c_im( b[i] ) *= -dn;
    }
  }

  p->timings[7] = MPI_Wtime();

  return 0;
}	/* HPCC_pzfft1d */

/*
 * End of pzfft1d.c
 */

/*
 * Begin of fft235.c
 */
static void
fft2(fftw_complex *a, fftw_complex *b, int m) {
  int i, lda, ldb;
  double x0, x1, y0, y1;

  lda = m;
  ldb = m;

  for (i = 0; i < m; ++i) {
    x0 = c_re( ARR2D( a, i, 0, lda ) );
    y0 = c_im( ARR2D( a, i, 0, lda ) );
    x1 = c_re( ARR2D( a, i, 1, lda ) );
    y1 = c_im( ARR2D( a, i, 1, lda ) );
    c_re( ARR2D( b, i, 0, ldb ) ) = x0 + x1;
    c_im( ARR2D( b, i, 0, ldb ) ) = y0 + y1;
    c_re( ARR2D( b, i, 1, ldb ) ) = x0 - x1;
    c_im( ARR2D( b, i, 1, ldb ) ) = y0 - y1;
  }
}

static void
fft4a(fftw_complex *a, fftw_complex *b, fftw_complex *w, int l) {
  int j, lda, ldb;
  double wr1, wr2, wr3, wi1, wi2, wi3;
  double x0, x1, x2, x3, y0, y1, y2, y3;

  lda = l;
  ldb = 4;

  for (j = 0; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2 = wr1*wr1 - wi1*wi1;
    wi2 = wr1*wi1 + wr1*wi1;
    wr3 = wr1*wr2 - wi1*wi2;
    wi3 = wr1*wi2 + wi1*wr2;

    x0 = c_re( ARR2D( a, j, 0, lda ) ) + c_re( ARR2D( a, j, 2, lda ) );
    y0 = c_im( ARR2D( a, j, 0, lda ) ) + c_im( ARR2D( a, j, 2, lda ) );
    x1 = c_re( ARR2D( a, j, 0, lda ) ) - c_re( ARR2D( a, j, 2, lda ) );
    y1 = c_im( ARR2D( a, j, 0, lda ) ) - c_im( ARR2D( a, j, 2, lda ) );

    x2 = c_re( ARR2D( a, j, 1, lda ) ) + c_re( ARR2D( a, j, 3, lda ) );
    y2 = c_im( ARR2D( a, j, 1, lda ) ) + c_im( ARR2D( a, j, 3, lda ) );
    x3 = c_im( ARR2D( a, j, 1, lda ) ) - c_im( ARR2D( a, j, 3, lda ) );
    y3 = c_re( ARR2D( a, j, 3, lda ) ) - c_re( ARR2D( a, j, 1, lda ) );

    c_re( ARR2D( b, 0, j, ldb ) ) = x0 + x2;
    c_im( ARR2D( b, 0, j, ldb ) ) = y0 + y2;
    c_re( ARR2D( b, 2, j, ldb ) ) = wr2 * (x0-x2) - wi2 * (y0-y2);
    c_im( ARR2D( b, 2, j, ldb ) ) = wr2 * (y0-y2) + wi2 * (x0-x2);
    c_re( ARR2D( b, 1, j, ldb ) ) = wr1 * (x1+x3) - wi1 * (y1+y3);
    c_im( ARR2D( b, 1, j, ldb ) ) = wr1 * (y1+y3) + wi1 * (x1+x3);
    c_re( ARR2D( b, 3, j, ldb ) ) = wr3 * (x1-x3) - wi3 * (y1-y3);
    c_im( ARR2D( b, 3, j, ldb ) ) = wr3 * (y1-y3) + wi3 * (x1-x3);
  }
}

static void
fft4b(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  int i, j, lda1, lda2, ldb1, ldb2;
  double x0, x1, x2, x3, y0, y1, y2, y3;
  double wr1, wr2, wr3, wi1, wi2, wi3;

  lda1 = m;
  lda2 = l;
  ldb1 = m;
  ldb2 = 4;

  for (i = 0; i < m; ++i) {
    x0 = c_re( ARR3D( a, i, 0, 0, lda1, lda2 ) ) + c_re( ARR3D( a, i, 0, 2, lda1, lda2 ) );
    y0 = c_im( ARR3D( a, i, 0, 0, lda1, lda2 ) ) + c_im( ARR3D( a, i, 0, 2, lda1, lda2 ) );
    x1 = c_re( ARR3D( a, i, 0, 0, lda1, lda2 ) ) - c_re( ARR3D( a, i, 0, 2, lda1, lda2 ) );
    y1 = c_im( ARR3D( a, i, 0, 0, lda1, lda2 ) ) - c_im( ARR3D( a, i, 0, 2, lda1, lda2 ) );

    x2 = c_re( ARR3D( a, i, 0, 1, lda1, lda2 ) ) + c_re( ARR3D( a, i, 0, 3, lda1, lda2 ) );
    y2 = c_im( ARR3D( a, i, 0, 1, lda1, lda2 ) ) + c_im( ARR3D( a, i, 0, 3, lda1, lda2 ) );
    x3 = c_im( ARR3D( a, i, 0, 1, lda1, lda2 ) ) - c_im( ARR3D( a, i, 0, 3, lda1, lda2 ) );
    y3 = c_re( ARR3D( a, i, 0, 3, lda1, lda2 ) ) - c_re( ARR3D( a, i, 0, 1, lda1, lda2 ) );

    c_re( ARR3D( b, i, 0, 0, ldb1, ldb2 ) ) = x0 + x2;
    c_im( ARR3D( b, i, 0, 0, ldb1, ldb2 ) ) = y0 + y2;
    c_re( ARR3D( b, i, 2, 0, ldb1, ldb2 ) ) = x0 - x2;
    c_im( ARR3D( b, i, 2, 0, ldb1, ldb2 ) ) = y0 - y2;

    c_re( ARR3D( b, i, 1, 0, ldb1, ldb2 ) ) = x1 + x3;
    c_im( ARR3D( b, i, 1, 0, ldb1, ldb2 ) ) = y1 + y3;
    c_re( ARR3D( b, i, 3, 0, ldb1, ldb2 ) ) = x1 - x3;
    c_im( ARR3D( b, i, 3, 0, ldb1, ldb2 ) ) = y1 - y3;
  }

  for (j = 1; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2 = wr1*wr1 - wi1*wi1;
    wi2 = wr1*wi1 + wr1*wi1;
    wr3 = wr1*wr2 - wi1*wi2;
    wi3 = wr1*wi2 + wi1*wr2;

    for (i = 0; i < m; ++i) {
      x0 = c_re( ARR3D( a, i, j, 0, lda1, lda2 ) ) + c_re( ARR3D( a, i, j, 2, lda1, lda2 ) );
      y0 = c_im( ARR3D( a, i, j, 0, lda1, lda2 ) ) + c_im( ARR3D( a, i, j, 2, lda1, lda2 ) );
      x1 = c_re( ARR3D( a, i, j, 0, lda1, lda2 ) ) - c_re( ARR3D( a, i, j, 2, lda1, lda2 ) );
      y1 = c_im( ARR3D( a, i, j, 0, lda1, lda2 ) ) - c_im( ARR3D( a, i, j, 2, lda1, lda2 ) );

      x2 = c_re( ARR3D( a, i, j, 1, lda1, lda2 ) ) + c_re( ARR3D( a, i, j, 3, lda1, lda2 ) );
      y2 = c_im( ARR3D( a, i, j, 1, lda1, lda2 ) ) + c_im( ARR3D( a, i, j, 3, lda1, lda2 ) );
      x3 = c_im( ARR3D( a, i, j, 1, lda1, lda2 ) ) - c_im( ARR3D( a, i, j, 3, lda1, lda2 ) );
      y3 = c_re( ARR3D( a, i, j, 3, lda1, lda2 ) ) - c_re( ARR3D( a, i, j, 1, lda1, lda2 ) );

      c_re( ARR3D( b, i, 0, j, ldb1, ldb2 ) ) = x0 + x2;
      c_im( ARR3D( b, i, 0, j, ldb1, ldb2 ) ) = y0 + y2;
      c_re( ARR3D( b, i, 2, j, ldb1, ldb2 ) ) = wr2 * (x0-x2) - wi2 * (y0-y2);
      c_im( ARR3D( b, i, 2, j, ldb1, ldb2 ) ) = wr2 * (y0-y2) + wi2 * (x0-x2);
      c_re( ARR3D( b, i, 1, j, ldb1, ldb2 ) ) = wr1 * (x1+x3) - wi1 * (y1+y3);
      c_im( ARR3D( b, i, 1, j, ldb1, ldb2 ) ) = wr1 * (y1+y3) + wi1 * (x1+x3);
      c_re( ARR3D( b, i, 3, j, ldb1, ldb2 ) ) = wr3 * (x1-x3) - wi3 * (y1-y3);
      c_im( ARR3D( b, i, 3, j, ldb1, ldb2 ) ) = wr3 * (y1-y3) + wi3 * (x1-x3);
    }
  }
}

static void
fft8a(fftw_complex *a, fftw_complex *b, fftw_complex *w, int l) {
  int j, lda, ldb;
  double x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7;
  double wr1, wr2, wr3, wr4, wr5, wr6, wr7, wi1, wi2, wi3, wi4, wi5, wi6, wi7;
  double u0, u1, u2, u3, v0, v1, v2, v3;
  double c81 = 0.70710678118654752;

  lda = l;
  ldb = 8;

  for (j = 0; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2 = wr1*wr1 - wi1*wi1;
    wi2 = wr1*wi1 + wr1*wi1;
    wr3 = wr1*wr2 - wi1*wi2;
    wi3 = wr1*wi2 + wi1*wr2;
    wr4 = wr2*wr2 - wi2*wi2;
    wi4 = wr2*wi2 + wr2*wi2;
    wr5 = wr2*wr3 - wi2*wi3;
    wi5 = wr2*wi3 + wi2*wr3;
    wr6 = wr3*wr3 - wi3*wi3;
    wi6 = wr3*wi3 + wr3*wi3;
    wr7 = wr3*wr4 - wi3*wi4;
    wi7 = wr3*wi4 + wi3*wr4;

    x0 = c_re( ARR2D( a, j, 0, lda ) ) + c_re( ARR2D( a, j, 4, lda ) );
    y0 = c_im( ARR2D( a, j, 0, lda ) ) + c_im( ARR2D( a, j, 4, lda ) );
    x1 = c_re( ARR2D( a, j, 0, lda ) ) - c_re( ARR2D( a, j, 4, lda ) );
    y1 = c_im( ARR2D( a, j, 0, lda ) ) - c_im( ARR2D( a, j, 4, lda ) );

    x2 = c_re( ARR2D( a, j, 2, lda ) ) + c_re( ARR2D( a, j, 6, lda ) );
    y2 = c_im( ARR2D( a, j, 2, lda ) ) + c_im( ARR2D( a, j, 6, lda ) );
    x3 = c_im( ARR2D( a, j, 2, lda ) ) - c_im( ARR2D( a, j, 6, lda ) );
    y3 = c_re( ARR2D( a, j, 6, lda ) ) - c_re( ARR2D( a, j, 2, lda ) );

    u0 = x0 + x2;
    v0 = y0 + y2;
    u1 = x0 - x2;
    v1 = y0 - y2;

    x4 = c_re( ARR2D( a, j, 1, lda ) ) + c_re( ARR2D( a, j, 5, lda ) );
    y4 = c_im( ARR2D( a, j, 1, lda ) ) + c_im( ARR2D( a, j, 5, lda ) );
    x5 = c_re( ARR2D( a, j, 1, lda ) ) - c_re( ARR2D( a, j, 5, lda ) );
    y5 = c_im( ARR2D( a, j, 1, lda ) ) - c_im( ARR2D( a, j, 5, lda ) );

    x6 = c_re( ARR2D( a, j, 3, lda ) ) + c_re( ARR2D( a, j, 7, lda ) );
    y6 = c_im( ARR2D( a, j, 3, lda ) ) + c_im( ARR2D( a, j, 7, lda ) );
    x7 = c_re( ARR2D( a, j, 3, lda ) ) - c_re( ARR2D( a, j, 7, lda ) );
    y7 = c_im( ARR2D( a, j, 3, lda ) ) - c_im( ARR2D( a, j, 7, lda ) );

    u2 = x4 + x6;
    v2 = y4 + y6;
    u3 = y4 - y6;
    v3 = x6 - x4;

    c_re( ARR2D( b, 0, j, ldb ) ) = u0 + u2;
    c_im( ARR2D( b, 0, j, ldb ) ) = v0 + v2;
    c_re( ARR2D( b, 4, j, ldb ) ) = wr4 * (u0-u2) - wi4 * (v0-v2);
    c_im( ARR2D( b, 4, j, ldb ) ) = wr4 * (v0-v2) + wi4 * (u0-u2);
    c_re( ARR2D( b, 2, j, ldb ) ) = wr2 * (u1+u3) - wi2 * (v1+v3);
    c_im( ARR2D( b, 2, j, ldb ) ) = wr2 * (v1+v3) + wi2 * (u1+u3);
    c_re( ARR2D( b, 6, j, ldb ) ) = wr6 * (u1-u3) - wi6 * (v1-v3);
    c_im( ARR2D( b, 6, j, ldb ) ) = wr6 * (v1-v3) + wi6 * (u1-u3);

    u0 = x1 + c81 * (x5 - x7);
    v0 = y1 + c81 * (y5 - y7);
    u1 = x1 - c81 * (x5 - x7);
    v1 = y1 - c81 * (y5 - y7);
    u2 = x3 + c81 * (y5 + y7);
    v2 = y3 - c81 * (x5 + x7);
    u3 = x3 - c81 * (y5 + y7);
    v3 = y3 + c81 * (x5 + x7);

    c_re( ARR2D( b, 1, j, ldb ) ) = wr1 * (u0+u2) - wi1 * (v0+v2);
    c_im( ARR2D( b, 1, j, ldb ) ) = wr1 * (v0+v2) + wi1 * (u0+u2);
    c_re( ARR2D( b, 5, j, ldb ) ) = wr5 * (u1+u3) - wi5 * (v1+v3);
    c_im( ARR2D( b, 5, j, ldb ) ) = wr5 * (v1+v3) + wi5 * (u1+u3);
    c_re( ARR2D( b, 3, j, ldb ) ) = wr3 * (u1-u3) - wi3 * (v1-v3);
    c_im( ARR2D( b, 3, j, ldb ) ) = wr3 * (v1-v3) + wi3 * (u1-u3);
    c_re( ARR2D( b, 7, j, ldb ) ) = wr7 * (u0-u2) - wi7 * (v0-v2);
    c_im( ARR2D( b, 7, j, ldb ) ) = wr7 * (v0-v2) + wi7 * (u0-u2);
  }
}

static void
fft8b(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  int i, j, lda1, lda2, ldb1, ldb2;
  double x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6, y7;
  double wr1, wr2, wr3, wr4, wr5, wr6, wr7, wi1, wi2, wi3, wi4, wi5, wi6, wi7;
  double u0, u1, u2, u3, v0, v1, v2, v3;
  double c81 = 0.70710678118654752;

  lda1 = m;
  lda2 = l;
  ldb1 = m;
  ldb2 = 8;

  for (i = 0; i < m; ++i) {
    x0 = c_re( ARR3D( a, i, 0, 0, lda1, lda2 ) ) + c_re( ARR3D( a, i, 0, 4, lda1, lda2 ) );
    y0 = c_im( ARR3D( a, i, 0, 0, lda1, lda2 ) ) + c_im( ARR3D( a, i, 0, 4, lda1, lda2 ) );
    x1 = c_re( ARR3D( a, i, 0, 0, lda1, lda2 ) ) - c_re( ARR3D( a, i, 0, 4, lda1, lda2 ) );
    y1 = c_im( ARR3D( a, i, 0, 0, lda1, lda2 ) ) - c_im( ARR3D( a, i, 0, 4, lda1, lda2 ) );

    x2 = c_re( ARR3D( a, i, 0, 2, lda1, lda2 ) ) + c_re( ARR3D( a, i, 0, 6, lda1, lda2 ) );
    y2 = c_im( ARR3D( a, i, 0, 2, lda1, lda2 ) ) + c_im( ARR3D( a, i, 0, 6, lda1, lda2 ) );
    x3 = c_im( ARR3D( a, i, 0, 2, lda1, lda2 ) ) - c_im( ARR3D( a, i, 0, 6, lda1, lda2 ) );
    y3 = c_re( ARR3D( a, i, 0, 6, lda1, lda2 ) ) - c_re( ARR3D( a, i, 0, 2, lda1, lda2 ) );

    u0 = x0 + x2;
    v0 = y0 + y2;
    u1 = x0 - x2;
    v1 = y0 - y2;

    x4 = c_re( ARR3D( a, i, 0, 1, lda1, lda2 ) ) + c_re( ARR3D( a, i, 0, 5, lda1, lda2 ) );
    y4 = c_im( ARR3D( a, i, 0, 1, lda1, lda2 ) ) + c_im( ARR3D( a, i, 0, 5, lda1, lda2 ) );
    x5 = c_re( ARR3D( a, i, 0, 1, lda1, lda2 ) ) - c_re( ARR3D( a, i, 0, 5, lda1, lda2 ) );
    y5 = c_im( ARR3D( a, i, 0, 1, lda1, lda2 ) ) - c_im( ARR3D( a, i, 0, 5, lda1, lda2 ) );

    x6 = c_re( ARR3D( a, i, 0, 3, lda1, lda2 ) ) + c_re( ARR3D( a, i, 0, 7, lda1, lda2 ) );
    y6 = c_im( ARR3D( a, i, 0, 3, lda1, lda2 ) ) + c_im( ARR3D( a, i, 0, 7, lda1, lda2 ) );
    x7 = c_re( ARR3D( a, i, 0, 3, lda1, lda2 ) ) - c_re( ARR3D( a, i, 0, 7, lda1, lda2 ) );
    y7 = c_im( ARR3D( a, i, 0, 3, lda1, lda2 ) ) - c_im( ARR3D( a, i, 0, 7, lda1, lda2 ) );

    u2 = x4 + x6;
    v2 = y4 + y6;
    u3 = y4 - y6;
    v3 = x6 - x4;

    c_re( ARR3D( b, i, 0, 0, ldb1, ldb2 ) ) = u0 + u2;
    c_im( ARR3D( b, i, 0, 0, ldb1, ldb2 ) ) = v0 + v2;
    c_re( ARR3D( b, i, 4, 0, ldb1, ldb2 ) ) = u0 - u2;
    c_im( ARR3D( b, i, 4, 0, ldb1, ldb2 ) ) = v0 - v2;

    c_re( ARR3D( b, i, 2, 0, ldb1, ldb2 ) ) = u1 + u3;
    c_im( ARR3D( b, i, 2, 0, ldb1, ldb2 ) ) = v1 + v3;
    c_re( ARR3D( b, i, 6, 0, ldb1, ldb2 ) ) = u1 - u3;
    c_im( ARR3D( b, i, 6, 0, ldb1, ldb2 ) ) = v1 - v3;

    u0 = x1 + c81 * (x5 - x7);
    v0 = y1 + c81 * (y5 - y7);
    u1 = x1 - c81 * (x5 - x7);
    v1 = y1 - c81 * (y5 - y7);
    u2 = x3 + c81 * (y5 + y7);
    v2 = y3 - c81 * (x5 + x7);
    u3 = x3 - c81 * (y5 + y7);
    v3 = y3 + c81 * (x5 + x7);

    c_re( ARR3D( b, i, 1, 0, ldb1, ldb2 ) ) = u0 + u2;
    c_im( ARR3D( b, i, 1, 0, ldb1, ldb2 ) ) = v0 + v2;
    c_re( ARR3D( b, i, 5, 0, ldb1, ldb2 ) ) = u1 + u3;
    c_im( ARR3D( b, i, 5, 0, ldb1, ldb2 ) ) = v1 + v3;

    c_re( ARR3D( b, i, 3, 0, ldb1, ldb2 ) ) = u1 - u3;
    c_im( ARR3D( b, i, 3, 0, ldb1, ldb2 ) ) = v1 - v3;
    c_re( ARR3D( b, i, 7, 0, ldb1, ldb2 ) ) = u0 - u2;
    c_im( ARR3D( b, i, 7, 0, ldb1, ldb2 ) ) = v0 - v2;
  }

  for (j = 1; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2 = wr1*wr1 - wi1*wi1;
    wi2 = wr1*wi1 + wr1*wi1;
    wr3 = wr1*wr2 - wi1*wi2;
    wi3 = wr1*wi2 + wi1*wr2;
    wr4 = wr2*wr2 - wi2*wi2;
    wi4 = wr2*wi2 + wr2*wi2;
    wr5 = wr2*wr3 - wi2*wi3;
    wi5 = wr2*wi3 + wi2*wr3;
    wr6 = wr3*wr3 - wi3*wi3;
    wi6 = wr3*wi3 + wr3*wi3;
    wr7 = wr3*wr4 - wi3*wi4;
    wi7 = wr3*wi4 + wi3*wr4;

    for (i = 0; i < m; ++i) {
      x0 = c_re( ARR3D( a, i, j, 0, lda1, lda2 ) ) + c_re( ARR3D( a, i, j, 4, lda1, lda2 ) );
      y0 = c_im( ARR3D( a, i, j, 0, lda1, lda2 ) ) + c_im( ARR3D( a, i, j, 4, lda1, lda2 ) );
      x1 = c_re( ARR3D( a, i, j, 0, lda1, lda2 ) ) - c_re( ARR3D( a, i, j, 4, lda1, lda2 ) );
      y1 = c_im( ARR3D( a, i, j, 0, lda1, lda2 ) ) - c_im( ARR3D( a, i, j, 4, lda1, lda2 ) );

      x2 = c_re( ARR3D( a, i, j, 2, lda1, lda2 ) ) + c_re( ARR3D( a, i, j, 6, lda1, lda2 ) );
      y2 = c_im( ARR3D( a, i, j, 2, lda1, lda2 ) ) + c_im( ARR3D( a, i, j, 6, lda1, lda2 ) );
      x3 = c_im( ARR3D( a, i, j, 2, lda1, lda2 ) ) - c_im( ARR3D( a, i, j, 6, lda1, lda2 ) );
      y3 = c_re( ARR3D( a, i, j, 6, lda1, lda2 ) ) - c_re( ARR3D( a, i, j, 2, lda1, lda2 ) );

      u0 = x0 + x2;
      v0 = y0 + y2;
      u1 = x0 - x2;
      v1 = y0 - y2;

      x4 = c_re( ARR3D( a, i, j, 1, lda1, lda2 ) ) + c_re( ARR3D( a, i, j, 5, lda1, lda2 ) );
      y4 = c_im( ARR3D( a, i, j, 1, lda1, lda2 ) ) + c_im( ARR3D( a, i, j, 5, lda1, lda2 ) );
      x5 = c_re( ARR3D( a, i, j, 1, lda1, lda2 ) ) - c_re( ARR3D( a, i, j, 5, lda1, lda2 ) );
      y5 = c_im( ARR3D( a, i, j, 1, lda1, lda2 ) ) - c_im( ARR3D( a, i, j, 5, lda1, lda2 ) );

      x6 = c_re( ARR3D( a, i, j, 3, lda1, lda2 ) ) + c_re( ARR3D( a, i, j, 7, lda1, lda2 ) );
      y6 = c_im( ARR3D( a, i, j, 3, lda1, lda2 ) ) + c_im( ARR3D( a, i, j, 7, lda1, lda2 ) );
      x7 = c_re( ARR3D( a, i, j, 3, lda1, lda2 ) ) - c_re( ARR3D( a, i, j, 7, lda1, lda2 ) );
      y7 = c_im( ARR3D( a, i, j, 3, lda1, lda2 ) ) - c_im( ARR3D( a, i, j, 7, lda1, lda2 ) );

      u2 = x4 + x6;
      v2 = y4 + y6;
      u3 = y4 - y6;
      v3 = x6 - x4;

      c_re( ARR3D( b, i, 0, j, ldb1, ldb2 ) ) = u0 + u2;
      c_im( ARR3D( b, i, 0, j, ldb1, ldb2 ) ) = v0 + v2;
      c_re( ARR3D( b, i, 4, j, ldb1, ldb2 ) ) = wr4 * (u0-u2) - wi4 * (v0-v2);
      c_im( ARR3D( b, i, 4, j, ldb1, ldb2 ) ) = wr4 * (v0-v2) + wi4 * (u0-u2);
      c_re( ARR3D( b, i, 2, j, ldb1, ldb2 ) ) = wr2 * (u1+u3) - wi2 * (v1+v3);
      c_im( ARR3D( b, i, 2, j, ldb1, ldb2 ) ) = wr2 * (v1+v3) + wi2 * (u1+u3);
      c_re( ARR3D( b, i, 6, j, ldb1, ldb2 ) ) = wr6 * (u1-u3) - wi6 * (v1-v3);
      c_im( ARR3D( b, i, 6, j, ldb1, ldb2 ) ) = wr6 * (v1-v3) + wi6 * (u1-u3);

      u0 = x1 + c81 * (x5 - x7);
      v0 = y1 + c81 * (y5 - y7);
      u1 = x1 - c81 * (x5 - x7);
      v1 = y1 - c81 * (y5 - y7);
      u2 = x3 + c81 * (y5 + y7);
      v2 = y3 - c81 * (x5 + x7);
      u3 = x3 - c81 * (y5 + y7);
      v3 = y3 + c81 * (x5 + x7);

      c_re( ARR3D( b, i, 1, j, ldb1, ldb2 ) ) = wr1 * (u0+u2) - wi1 * (v0+v2);
      c_im( ARR3D( b, i, 1, j, ldb1, ldb2 ) ) = wr1 * (v0+v2) + wi1 * (u0+u2);
      c_re( ARR3D( b, i, 5, j, ldb1, ldb2 ) ) = wr5 * (u1+u3) - wi5 * (v1+v3);
      c_im( ARR3D( b, i, 5, j, ldb1, ldb2 ) ) = wr5 * (v1+v3) + wi5 * (u1+u3);
      c_re( ARR3D( b, i, 3, j, ldb1, ldb2 ) ) = wr3 * (u1-u3) - wi3 * (v1-v3);
      c_im( ARR3D( b, i, 3, j, ldb1, ldb2 ) ) = wr3 * (v1-v3) + wi3 * (u1-u3);
      c_re( ARR3D( b, i, 7, j, ldb1, ldb2 ) ) = wr7 * (u0-u2) - wi7 * (v0-v2);
      c_im( ARR3D( b, i, 7, j, ldb1, ldb2 ) ) = wr7 * (v0-v2) + wi7 * (u0-u2);
    }
  }
}

static void
fft3a(fftw_complex *a, fftw_complex *b, fftw_complex *w, int l) {
  int j;
  double x0, x1, x2;
  double y0, y1, y2;
  double wr1, wr2;
  double wi1, wi2;
  double c31 = 0.86602540378443865, c32 = 0.5;

  for (j = 0; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2=wr1*wr1-wi1*wi1;
    wi2=wr1*wi1+wr1*wi1;
    x0 = c_re( ARR2D( a, j, 1, l ) ) + c_re( ARR2D( a, j, 2, l ) );
    y0 = c_im( ARR2D( a, j, 1, l ) ) + c_im( ARR2D( a, j, 2, l ) );
    x1 = c_re( ARR2D( a, j, 0, l ) ) - c32 * x0;
    y1 = c_im( ARR2D( a, j, 0, l ) ) - c32 * y0;
    x2 = c31 * ( c_im( ARR2D( a, j, 1, l ) ) - c_im( ARR2D( a, j, 2, l ) ));
    y2 = c31 * ( c_re( ARR2D( a, j, 2, l ) ) - c_re( ARR2D( a, j, 1, l ) ));
    c_re( ARR2D( b, 0, j, 3 ) ) = c_re( ARR2D( a, j, 0, l ) ) + x0;
    c_im( ARR2D( b, 0, j, 3 ) ) = c_im( ARR2D( a, j, 0, l ) ) + y0;
    c_re( ARR2D( b, 1, j, 3 ) ) = wr1*(x1+x2)-wi1*(y1+y2);
    c_im( ARR2D( b, 1, j, 3 ) ) = wr1*(y1+y2)+wi1*(x1+x2);
    c_re( ARR2D( b, 2, j, 3 ) ) = wr2*(x1-x2)-wi2*(y1-y2);
    c_im( ARR2D( b, 2, j, 3 ) ) = wr2*(y1-y2)+wi2*(x1-x2);
  }
}

static void
fft3b(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  int i, j;
  double x0, x1, x2;
  double y0, y1, y2;
  double wr1, wr2;
  double wi1, wi2;
  double c31 = 0.86602540378443865, c32 = 0.5;

  for (i = 0; i < m; ++i) {
    x0 = c_re( ARR3D( a, i, 0, 1, m, l ) ) + c_re( ARR3D( a, i, 0, 2, m, l ) );
    y0 = c_im( ARR3D( a, i, 0, 1, m, l ) ) + c_im( ARR3D( a, i, 0, 2, m, l ) );
    x1 = c_re( ARR3D( a, i, 0, 0, m, l ) ) - c32 * x0;
    y1 = c_im( ARR3D( a, i, 0, 0, m, l ) ) - c32 * y0;
    x2 = c31 * (c_im( ARR3D( a, i, 0, 1, m, l ) ) - c_im( ARR3D( a, i, 0, 2, m, l ) ));
    y2 = c31 * (c_re( ARR3D( a, i, 0, 2, m, l ) ) - c_re( ARR3D( a, i, 0, 1, m, l ) ));
    c_re( ARR3D( b, i, 0, 0, m, 3 ) ) = c_re( ARR3D( a, i, 0, 0, m, l ) ) + x0;
    c_im( ARR3D( b, i, 0, 0, m, 3 ) ) = c_im( ARR3D( a, i, 0, 0, m, l ) ) + y0;
    c_re( ARR3D( b, i, 1, 0, m, 3 ) ) = x1 + x2;
    c_im( ARR3D( b, i, 1, 0, m, 3 ) ) = y1 + y2;
    c_re( ARR3D( b, i, 2, 0, m, 3 ) ) = x1 - x2;
    c_im( ARR3D( b, i, 2, 0, m, 3 ) ) = y1 - y2;
  }

  for (j = 1; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2=wr1*wr1-wi1*wi1;
    wi2=wr1*wi1+wr1*wi1;
    for (i = 0; i < m; ++i) {
      x0 = c_re( ARR3D( a, i, j, 1, m, l ) ) + c_re( ARR3D( a, i, j, 2, m, l ) );
      y0 = c_im( ARR3D( a, i, j, 1, m, l ) ) + c_im( ARR3D( a, i, j, 2, m, l ) );
      x1 = c_re( ARR3D( a, i, j, 0, m, l ) ) - c32 * x0;
      y1 = c_im( ARR3D( a, i, j, 0, m, l ) ) - c32 * y0;
      x2 = c31 * (c_im( ARR3D( a, i, j, 1, m, l ) ) - c_im( ARR3D( a, i, j, 2, m, l ) ));
      y2 = c31 * (c_re( ARR3D( a, i, j, 2, m, l ) ) - c_re( ARR3D( a, i, j, 1, m, l ) ));
      c_re( ARR3D( b, i, 0, j, m, 3 ) ) = c_re( ARR3D( a, i, j, 0, m, l ) ) + x0;
      c_im( ARR3D( b, i, 0, j, m, 3 ) ) = c_im( ARR3D( a, i, j, 0, m, l ) ) + y0;
      c_re( ARR3D( b, i, 1, j, m, 3 ) ) = wr1*(x1+x2)-wi1*(y1+y2);
      c_im( ARR3D( b, i, 1, j, m, 3 ) ) = wr1*(y1+y2)+wi1*(x1+x2);
      c_re( ARR3D( b, i, 2, j, m, 3 ) ) = wr2*(x1-x2)-wi2*(y1-y2);
      c_im( ARR3D( b, i, 2, j, m, 3 ) ) = wr2*(y1-y2)+wi2*(x1-x2);
    }
  }
}

static void
fft5a(fftw_complex *a, fftw_complex *b, fftw_complex *w, int l) {
  int j;
  double wr1, wr2, wr3, wr4;
  double wi1, wi2, wi3, wi4;
  double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
  double y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10;
  double c51 = 0.95105651629515357, c52 = 0.61803398874989485;
  double c53 = 0.55901699437494742, c54 = 0.25;

  for (j = 0; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2=wr1*wr1-wi1*wi1;
    wi2=wr1*wi1+wr1*wi1;
    wr3=wr1*wr2-wi1*wi2;
    wi3=wr1*wi2+wi1*wr2;
    wr4=wr2*wr2-wi2*wi2;
    wi4=wr2*wi2+wr2*wi2;
    x0 = c_re( ARR2D( a, j, 1, l ) ) + c_re( ARR2D( a, j, 4, l ) );
    y0 = c_im( ARR2D( a, j, 1, l ) ) + c_im( ARR2D( a, j, 4, l ) );
    x1 = c_re( ARR2D( a, j, 2, l ) ) + c_re( ARR2D( a, j, 3, l ) );
    y1 = c_im( ARR2D( a, j, 2, l ) ) + c_im( ARR2D( a, j, 3, l ) );
    x2 = c51 * (c_re( ARR2D( a, j, 1, l ) ) - c_re( ARR2D( a, j, 4, l ) ));
    y2 = c51 * (c_im( ARR2D( a, j, 1, l ) ) - c_im( ARR2D( a, j, 4, l ) ));
    x3 = c51 * (c_re( ARR2D( a, j, 2, l ) ) - c_re( ARR2D( a, j, 3, l ) ));
    y3 = c51 * (c_im( ARR2D( a, j, 2, l ) ) - c_im( ARR2D( a, j, 3, l ) ));
    x4 = x0 + x1;
    y4 = y0 + y1;
    x5 = c53 * (x0-x1);
    y5 = c53 * (y0-y1);
    x6 = c_re( ARR2D( a, j, 0, l ) ) - c54 * x4;
    y6 = c_im( ARR2D( a, j, 0, l ) ) - c54 * y4;
    x7 = x6 + x5;
    y7 = y6 + y5;
    x8 = x6 - x5;
    y8 = y6 - y5;
    x9 = y2 + c52*y3;
    y9 = -x2 - c52*x3;
    x10 = c52*y2 - y3;
    y10 = x3 - c52*x2;
    c_re( ARR2D( b, 0, j, 5 ) ) = c_re( ARR2D( a, j, 0, l ) ) + x4;
    c_im( ARR2D( b, 0, j, 5 ) ) = c_im( ARR2D( a, j, 0, l ) ) + y4;
    c_re( ARR2D( b, 1, j, 5 ) ) = wr1 * (x7+x9) - wi1 * (y7+y9);
    c_im( ARR2D( b, 1, j, 5 ) ) = wr1 * (y7+y9) + wi1 * (x7+x9);
    c_re( ARR2D( b, 2, j, 5 ) ) = wr2 * (x8+x10) - wi2 * (y8+y10);
    c_im( ARR2D( b, 2, j, 5 ) ) = wr2 * (y8+y10) + wi2 * (x8+x10);
    c_re( ARR2D( b, 3, j, 5 ) ) = wr3 * (x8-x10) - wi3 * (y8-y10);
    c_im( ARR2D( b, 3, j, 5 ) ) = wr3 * (y8-y10) + wi3 * (x8-x10);
    c_re( ARR2D( b, 4, j, 5 ) ) = wr4 * (x7-x9) - wi4 * (y7-y9);
    c_im( ARR2D( b, 4, j, 5 ) ) = wr4 * (y7-y9) + wi4 * (x7-x9);
  }
}

static void
fft5b(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  int i, j;
  double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
  double y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10;
  double wr1, wr2, wr3, wr4;
  double wi1, wi2, wi3, wi4;
  double c51 = 0.95105651629515357, c52 = 0.61803398874989485;
  double c53 = 0.55901699437494742, c54 = 0.25;

  for (i = 0; i < m; ++i) {
    x0 = c_re( ARR3D( a, i, 0, 1, m, l ) ) + c_re( ARR3D( a, i, 0, 4, m, l ) );
    y0 = c_im( ARR3D( a, i, 0, 1, m, l ) ) + c_im( ARR3D( a, i, 0, 4, m, l ) );
    x1 = c_re( ARR3D( a, i, 0, 2, m, l ) ) + c_re( ARR3D( a, i, 0, 3, m, l ) );
    y1 = c_im( ARR3D( a, i, 0, 2, m, l ) ) + c_im( ARR3D( a, i, 0, 3, m, l ) );
    x2 = c51 * (c_re( ARR3D( a, i, 0, 1, m, l ) ) - c_re( ARR3D( a, i, 0, 4, m, l ) ));
    y2 = c51 * (c_im( ARR3D( a, i, 0, 1, m, l ) ) - c_im( ARR3D( a, i, 0, 4, m, l ) ));
    x3 = c51 * (c_re( ARR3D( a, i, 0, 2, m, l ) ) - c_re( ARR3D( a, i, 0, 3, m, l ) ));
    y3 = c51 * (c_im( ARR3D( a, i, 0, 2, m, l ) ) - c_im( ARR3D( a, i, 0, 3, m, l ) ));
    x4 = x0 + x1;
    y4 = y0 + y1;
    x5 = c53 * (x0-x1);
    y5 = c53 * (y0-y1);
    x6 = c_re( ARR3D( a, i, 0, 0, m, l ) ) - c54 * x4;
    y6 = c_im( ARR3D( a, i, 0, 0, m, l ) ) - c54 * y4;
    x7 = x6 + x5;
    y7 = y6 + y5;
    x8 = x6 - x5;
    y8 = y6 - y5;
    x9 = y2 + c52 * y3;
    y9 = -x2 - c52 * x3;
    x10 = c52 * y2 - y3;
    y10 = x3 - c52 * x2;
    c_re( ARR3D( b, i, 0, 0, m, 5 ) ) = c_re( ARR3D( a, i, 0, 0, m, l ) ) + x4;
    c_im( ARR3D( b, i, 0, 0, m, 5 ) ) = c_im( ARR3D( a, i, 0, 0, m, l ) ) + y4;
    c_re( ARR3D( b, i, 1, 0, m, 5 ) ) = x7 + x9;
    c_im( ARR3D( b, i, 1, 0, m, 5 ) ) = y7 + y9;
    c_re( ARR3D( b, i, 2, 0, m, 5 ) ) = x8 + x10;
    c_im( ARR3D( b, i, 2, 0, m, 5 ) ) = y8 + y10;
    c_re( ARR3D( b, i, 3, 0, m, 5 ) ) = x8 - x10;
    c_im( ARR3D( b, i, 3, 0, m, 5 ) ) = y8 - y10;
    c_re( ARR3D( b, i, 4, 0, m, 5 ) ) = x7 - x9;
    c_im( ARR3D( b, i, 4, 0, m, 5 ) ) = y7 - y9;
  }

  for (j = 1; j < l; ++j) {
    wr1 = c_re( w[j] );
    wi1 = c_im( w[j] );
    wr2 = wr1 * wr1 - wi1*wi1;
    wi2 = wr1 * wi1 + wr1*wi1;
    wr3 = wr1 * wr2 - wi1*wi2;
    wi3 = wr1 * wi2 + wi1*wr2;
    wr4 = wr2 * wr2 - wi2*wi2;
    wi4 = wr2 * wi2 + wr2*wi2;
    for (i = 0; i < m; ++i) {
      x0 = c_re( ARR3D( a, i, j, 1, m, l ) ) + c_re( ARR3D( a, i, j, 4, m, l ) );
      y0 = c_im( ARR3D( a, i, j, 1, m, l ) ) + c_im( ARR3D( a, i, j, 4, m, l ) );
      x1 = c_re( ARR3D( a, i, j, 2, m, l ) ) + c_re( ARR3D( a, i, j, 3, m, l ) );
      y1 = c_im( ARR3D( a, i, j, 2, m, l ) ) + c_im( ARR3D( a, i, j, 3, m, l ) );
      x2 = c51 * (c_re( ARR3D( a, i, j, 1, m, l ) ) - c_re( ARR3D( a, i, j, 4, m, l ) ));
      y2 = c51 * (c_im( ARR3D( a, i, j, 1, m, l ) ) - c_im( ARR3D( a, i, j, 4, m, l ) ));
      x3 = c51 * (c_re( ARR3D( a, i, j, 2, m, l ) ) - c_re( ARR3D( a, i, j, 3, m, l ) ));
      y3 = c51 * (c_im( ARR3D( a, i, j, 2, m, l ) ) - c_im( ARR3D( a, i, j, 3, m, l ) ));
      x4 = x0 + x1;
      y4 = y0 + y1;
      x5 = c53 * (x0-x1);
      y5 = c53 * (y0-y1);
      x6 = c_re( ARR3D( a, i, j, 0, m, l ) ) - c54*x4;
      y6 = c_im( ARR3D( a, i, j, 0, m, l ) ) - c54*y4;
      x7 = x6 + x5;
      y7 = y6 + y5;
      x8 = x6 - x5;
      y8 = y6 - y5;
      x9 = y2 + c52 * y3;
      y9 = -x2 - c52 * x3;
      x10 = c52*y2 - y3;
      y10 = x3 - c52*x2;
      c_re( ARR3D( b, i, 0, j, m, 5 ) ) = c_re( ARR3D( a, i, j, 0, m, l ) ) + x4;
      c_im( ARR3D( b, i, 0, j, m, 5 ) ) = c_im( ARR3D( a, i, j, 0, m, l ) ) + y4;
      c_re( ARR3D( b, i, 1, j, m, 5 ) ) = wr1*(x7+x9) - wi1*(y7+y9);
      c_im( ARR3D( b, i, 1, j, m, 5 ) ) = wr1*(y7+y9) + wi1*(x7+x9);
      c_re( ARR3D( b, i, 2, j, m, 5 ) ) = wr2*(x8+x10) - wi2*(y8+y10);
      c_im( ARR3D( b, i, 2, j, m, 5 ) ) = wr2*(y8+y10) + wi2*(x8+x10);
      c_re( ARR3D( b, i, 3, j, m, 5 ) ) = wr3*(x8-x10) - wi3*(y8-y10);
      c_im( ARR3D( b, i, 3, j, m, 5 ) ) = wr3*(y8-y10) + wi3*(x8-x10);
      c_re( ARR3D( b, i, 4, j, m, 5 ) ) = wr4*(x7-x9) - wi4*(y7-y9);
      c_im( ARR3D( b, i, 4, j, m, 5 ) ) = wr4*(y7-y9) + wi4*(x7-x9);
    }
  }
}

static void
fft3(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  if (1 == m)
    fft3a( a, b, w, l );
  else
    fft3b( a, b, w, m, l );
}

static void
fft4(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  if (1 == m)
    fft4a( a, b, w, l );
  else
    fft4b( a, b, w, m, l );
}

static void
fft5(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  if (1 == m)
    fft5a( a, b, w, l );
  else
    fft5b( a, b, w, m, l );
}

static void
fft8(fftw_complex *a, fftw_complex *b, fftw_complex *w, int m, int l) {
  if (1 == m)
    fft8a( a, b, w, l );
  else
    fft8b( a, b, w, m, l );
}

int
HPCC_fft235(fftw_complex *a, fftw_complex *b, fftw_complex *w, int n, const int *ip) {
  int j, k, l, m, key, kp4, kp8;

  if (ip[0] != 1) {
    kp4 = 2 - (ip[0] + 2) % 3;
    kp8 = (ip[0]-kp4) / 3;
  } else {
    kp4 = 0;
    kp8 = 0;
  }

  key = 1;
  j = 0;
  l = n;
  m = 1;

  for (k = 0; k < kp8; ++k) {
    l >>= 3; /* divide by 8 */

    if (l >= 2) {
      if (key > 0)
        fft8( a, b, w + j, m, l );
      else
        fft8( b, a, w + j, m, l );

      key = -key;
    } else {
      if (key > 0)
        fft8( a, a, w + j, m, l );
      else
        fft8( b, a, w + j, m, l );
    }
    m <<= 3; /* multiply by 8 */
    j += l;
  }

  for (k = 0; k < ip[2]; ++k) {
    l /= 5;

    if (l >= 2) {
      if (key > 0)
        fft5( a, b, w+j, m, l );
      else
        fft5( b, a, w+j, m, l );

      key = -key;
    } else {
      if (key > 0)
        fft5( a, a, w+j, m, l );
      else
        fft5( b, a, w+j, m, l );
    }

    m *= 5;
    j += l;
  }

  for (k = 0; k < kp4; ++k) {
    l >>= 2; /* divide by 4 */

    if (l >= 2) {
      if (key > 0)
        fft4( a, b, w + j, m, l );
      else
        fft4( b, a, w + j, m, l );

      key = -key;
    } else {
      if (key > 0)
        fft4( a, a, w + j, m, l );
      else
        fft4( b, a, w + j, m, l );
    }
    m <<= 2; /* multiply by 4 */
    j += l;
  }

  for (k = 0; k < ip[1]; ++k) {
    l /= 3;

    if (l >= 2) {
      if (key > 0)
        fft3( a, b, w+j, m, l );
      else
        fft3( b, a, w+j, m, l );

      key = -key;
    } else {
      if (key > 0)
        fft3( a, a, w+j, m, l );
      else
        fft3( b, a, w+j, m, l );
    }

    m *= 3;
    j += l;
  }

  if (ip[0] == 1) {
    if (key > 0)
      fft2( a, a, m );
    else
      fft2( b, a, m );
  }

  return 0;
}

static int
settbl0(fftw_complex *w, int m, int l) {
  int i;
  double pi2, px;

  pi2 = 8.0 * atan(1.0);
  px = -pi2 / m / l;

  for (i = 0; i < l; ++i) {
    c_re(w[i]) = cos(px * i);
    c_im(w[i]) = sin(px * i);
  }

  return 0;
}

int
HPCC_settbl(fftw_complex *w, int n) {
  int j, k, l, kp4, kp8;
  int ip[3];

  HPCC_factor235( n, ip );

  if (1 != ip[0]) {
    kp4 = 2 - (ip[0] + 2) % 3;
    kp8 = (ip[0]-kp4) / 3;
  } else {
    kp4 = 0;
    kp8 = 0;
  }

  j = 0;
  l = n;

  for (k = 0; k < kp8; ++k) {
    l >>= 3; /* divide by 8 */
    settbl0( w + j, 8, l );
    j += l;
  }

  for (k = 0; k < ip[2]; ++k) {
    l /= 5;
    settbl0( w + j, 5, l );
    j += l;
  }

  for (k = 0; k < kp4; ++k) {
    l >>= 2; /* divide by 4 */
    settbl0( w + j, 4, l );
    j += l;
  }

  for (k = 0; k < ip[1]; ++k) {
    l /= 3;
    settbl0( w + j, 3, l );
    j += l;
  }

  return 0;
}	/* settbl */

int
HPCC_factor235(int n, int *ip) {
  ip[0] = ip[1] = ip[2] = 0;

  if (n % 2 != 0 && n % 3 != 0 && n % 5 != 0)
    return 1;

  if (n <= 1)
    return 1;

  /* count all 2 factors */
  for (; n > 1 && ! (n & 1); n >>= 1)
    ip[0]++;

  /* count all 3 factors */
  for (; n > 1 && ! (n % 3); n /= 3)
    ip[1]++;

  /* count all 5 factors */
  for (; n > 1 && ! (n % 5); n /= 5)
    ip[2]++;

  if (n != 1)
    return 1;

  return 0;
}

int
HPCC_factor235_8(s64Int_t n, int *ip) {
  ip[0] = ip[1] = ip[2] = 0;

  if (n % 2 != 0 && n % 3 != 0 && n % 5 != 0)
    return 1;

  if (n <= 1)
    return 1;

  /* count all 2 factors */
  for (; n > 1 && ! (n & 1); n >>= 1)
    ip[0]++;

  /* count all 3 factors */
  for (; n > 1 && ! (n % 3); n /= 3)
    ip[1]++;

  /* count all 5 factors */
  for (; n > 1 && ! (n % 5); n /= 5)
    ip[2]++;

  if (n != 1)
    return 1;

  return 0;
}

/*
 * End of fft235.c
 */


static void
MPIFFT0(HPCC_Params *params, int doIO, FILE *outFile, MPI_Comm comm, int locN,
        double *UGflops, s64Int_t *Un, double *UmaxErr, int *Ufailure) {
  int commRank, commSize, failure, flags;
  s64Int_t i, n;
  s64Int_t locn, loc0, alocn, aloc0, tls;
  double maxErr, tmp1, tmp2, tmp3, t0, t1, t2, t3, Gflops;
  double deps;
  fftw_complex *inout, *work;
  fftw_mpi_plan p;
  hpcc_fftw_mpi_plan ip;
  int sAbort, rAbort;
#ifdef USING_FFTW
  int ilocn, iloc0, ialocn, ialoc0, itls;
#endif

  failure = 1;
  Gflops = -1.0;
  deps = HPL_dlamch( HPL_MACH_EPS );
  maxErr = 1.0 / deps;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &commRank );

  n = locN;

  /* number of processes have been factored out - need to put it back in */
  n *= commSize;

  n *= commSize; /* global vector size */

#ifdef USING_FFTW
  /* FFTW ver. 2 only supports vector sizes that fit in 'int' */
  if (n > (1<<30)-1+(1<<30)) {
#ifdef HPCC_FFTW_CHECK32
    goto no_plan;
#else
  if (doIO) {
    fprintf( outFile, "Warning: problem size too large: %ld*%d*%d\n", (long)(n / commSize / commSize), commSize, commSize );
  }
#endif
  }
#endif

#ifdef HPCC_FFTW_ESTIMATE
  flags = FFTW_ESTIMATE;
#else
  flags = FFTW_MEASURE;
#endif

  t1 = -MPI_Wtime();
  p = fftw_mpi_create_plan( comm, n, FFTW_FORWARD, flags );
  t1 += MPI_Wtime();

  if (! p) goto no_plan;

#ifdef USING_FFTW
  fftw_mpi_local_sizes( p, &ilocn, &iloc0, &ialocn, &ialoc0, &itls );
  locn = ilocn;
  loc0 = iloc0;
  alocn = ialocn;
  aloc0 = ialoc0;
  tls = itls;
#else
  fftw_mpi_local_sizes( p, &locn, &loc0, &alocn, &aloc0, &tls );
#endif

  inout = (fftw_complex *)HPCC_fftw_malloc( tls * (sizeof *inout) );
  work  = (fftw_complex *)HPCC_fftw_malloc( tls * (sizeof *work) );

  sAbort = 0;
  if (! inout || ! work) sAbort = 1;
  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, comm );
  if (rAbort > 0) {
    fftw_mpi_destroy_plan( p );
    goto comp_end;
  }

  /* Make sure that `inout' and `work' are initialized in parallel if using
     Open MP: this will ensure better placement of pages if first-touch policy
     is used by a distrubuted shared memory machine. */
#ifdef _OPENMP
#pragma omp parallel for
  for (i = 0; i < tls; ++i) {
    c_re( inout[i] ) = c_re( work[i] ) = 0.0;
    c_re( inout[i] ) = c_im( work[i] ) = 0.0;
  }
#endif

  t0 = -MPI_Wtime();
  HPCC_bcnrand( 2 * tls, 53 * commRank * 2 * tls, inout );
  t0 += MPI_Wtime();

  t2 = -MPI_Wtime();
  fftw_mpi( p, 1, inout, work );
  t2 += MPI_Wtime();

  fftw_mpi_destroy_plan( p );

  ip = HPCC_fftw_mpi_create_plan( comm, n, FFTW_BACKWARD, FFTW_ESTIMATE );

  if (ip) {
    t3 = -MPI_Wtime();
    HPCC_fftw_mpi( ip, 1, inout, work );
    t3 += MPI_Wtime();

    HPCC_fftw_mpi_destroy_plan( ip );
  }

  HPCC_bcnrand( 2 * tls, 53 * commRank * 2 * tls, work ); /* regenerate data */

  maxErr = 0.0;
  for (i = 0; i < locn; ++i) {
    tmp1 = c_re( inout[i] ) - c_re( work[i] );
    tmp2 = c_im( inout[i] ) - c_im( work[i] );
    tmp3 = sqrt( tmp1*tmp1 + tmp2*tmp2 );
    maxErr = maxErr >= tmp3 ? maxErr : tmp3;
  }
  MPI_Allreduce( &maxErr, UmaxErr, 1, MPI_DOUBLE, MPI_MAX, comm );
  maxErr = *UmaxErr;
  if (maxErr / log(n) / deps < params->test.thrsh) failure = 0;

  if (t2 > 0.0) Gflops = 1e-9 * (5.0 * n * log(n) / log(2.0)) / t2;

  if (doIO) {
    fprintf( outFile, "Number of nodes: %d\n", commSize );
    fprintf( outFile, "Vector size: %20.0f\n", tmp1 = (double)n );
    fprintf( outFile, "Generation time: %9.3f\n", t0 );
    fprintf( outFile, "Tuning: %9.3f\n", t1 );
    fprintf( outFile, "Computing: %9.3f\n", t2 );
    fprintf( outFile, "Inverse FFT: %9.3f\n", t3 );
    fprintf( outFile, "max(|x-x0|): %9.3e\n", maxErr );
    fprintf( outFile, "Gflop/s: %9.3f\n", Gflops );
  }

  comp_end:

  if (work) HPCC_fftw_free( work );
  if (inout) HPCC_fftw_free( inout );

  no_plan:

  *UGflops = Gflops;
  *Un = n;
  *UmaxErr = maxErr;
  *Ufailure = failure;
}



int
HPCC_LocalVectorSize(HPCC_Params *params, int vecCnt, size_t size, int pow2) {
  int flg2, maxIntBits2;

  /* this is the maximum power of 2 that that can be held in a signed integer (for a 4-byte
     integer, 2**31-1 is the maximum integer, so the maximum power of 2 is 30) */
  maxIntBits2 = sizeof(int) * 8 - 2;

  /* flg2 = floor(log2(params->HPLMaxProcMem / size / vecCnt)) */
  for (flg2 = 1; params->HPLMaxProcMem / size / vecCnt >> flg2; ++flg2)
    ; /* EMPTY */
  --flg2;

  if (flg2 <= maxIntBits2) {
    if (pow2)
      return 1 << flg2;

    return params->HPLMaxProcMem / size / vecCnt;
  }

  return 1 << maxIntBits2;
}

int
HPCC_ipow(int x, int p) {
  int i, r;

  if (1 == x || 0 == x) return x;
  if (0 == p) return 1;
  if (-1 == x) return (p & 1) ? -1 : 1;
  if (p < 0) return 0;
  r = 1;
  for (i = 0; i < p; i++) r *= x;
  return r;
}



int
main (int argc, char * argv[]) {

  HPCC_Params p;
  HPCC_Params* params = &p;


  params->HPLMaxProcMem = 1;
  params->test.thrsh = 1.0;

  int commRank, commSize;
  int locN, procCnt, isComputing, doIO, failure = 0;
  s64Int_t n;
  double Gflops = -1.0, maxErr = -1.0;
  MPI_Comm comm = MPI_COMM_WORLD;
  FILE *outFile;
  MPI_Init( &argc, &argv );

  MPI_Comm_size( MPI_COMM_WORLD, &commSize );
  MPI_Comm_rank( MPI_COMM_WORLD, &commRank );

  doIO = commRank == 0 ? 1 : 0;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) outFile = stderr;
  }

  /*
  There are two vectors of size 'n'/'commSize': inout, work,
  and internal work: 2*'n'/'commSize'; it's 4 vectors then.

  FFTE requires that the global vector size 'n' has to be at least
  as big as square of number of processes. The square is calculated
  in each factor independently. In other words, 'n' has to have
  at least twice as many 2 factors as the process count, twice as many
  3 factors and twice as many 5 factors.
  */

#ifdef HPCC_FFT_235
  locN = 0; procCnt = commSize + 1;
  do {
    int f[3];

    procCnt--;

    for ( ; procCnt > 1 && HPCC_factor235( procCnt, f ); procCnt--)
      ; /* EMPTY */

    /* Make sure the local vector size is greater than 0 */
    locN = HPCC_LocalVectorSize( params, 4*procCnt, sizeof(fftw_complex), 0 );
    for ( ; locN >= 1 && HPCC_factor235( locN, f ); locN--)
      ; /* EMPTY */
  } while (locN < 1);
#else
  /* Find power of two that is smaller or equal to number of processes */
  for (procCnt = 1; procCnt <= (commSize >> 1); procCnt <<= 1)
    ; /* EMPTY */

  /* Make sure the local vector size is greater than 0 */
  while (1) {
    locN = HPCC_LocalVectorSize( params, 4*procCnt, sizeof(fftw_complex), 1 );
    if (locN) break;
    procCnt >>= 1;
  }
#endif

  isComputing = commRank < procCnt ? 1 : 0;

  HPCC_fft_timings_forward = params->MPIFFTtimingsForward;
  HPCC_fft_timings_backward = params->MPIFFTtimingsBackward;

  if (commSize == procCnt)
    comm = MPI_COMM_WORLD;
  else
    MPI_Comm_split( MPI_COMM_WORLD, isComputing ? 0 : MPI_UNDEFINED, commRank, &comm );

  if (isComputing)
    MPIFFT0( params, doIO, outFile, comm, locN, &Gflops, &n, &maxErr, &failure );

  if (commSize != procCnt && isComputing && comm != MPI_COMM_NULL)
    MPI_Comm_free( &comm );

  params->MPIFFT_N = n;
  params->MPIFFT_Procs = procCnt;
  params->MPIFFT_maxErr = maxErr;

  MPI_Bcast( &Gflops, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  params->MPIFFTGflops = Gflops;

  params->FFTEnblk = FFTE_NBLK;
  params->FFTEnp = FFTE_NP;
  params->FFTEl2size = FFTE_L2SIZE;

  if (failure)
    params->Failure = 1;

  if (doIO) if (outFile != stderr) fclose( outFile );

  MPI_Finalize();
  return 0;
}
