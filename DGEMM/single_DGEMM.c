#include <math.h>
#include <float.h>
#include <time.h>

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

int HPCC_external_init(int argc, char *argv[], void *extdata);
int HPCC_external_finalize(int argc, char *argv[], void *extdata);

int HPCC_Init(HPCC_Params *params);
int HPCC_Finalize(HPCC_Params *params);
int HPCC_LocalVectorSize(HPCC_Params *params, int vecCnt, size_t size, int pow2);
int
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

//extern int HPL_main(int ARGC, char **ARGV, HPL_RuntimeData *rdata, int *failure);
float HPL_slamch (const HPL_T_MACH);
double HPCC_dweps();
float HPCC_sweps();

/*
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
extern int HPCC_MPIFFT(HPCC_Params *params); */

//int HPCC_TestFFT(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure);
int HPCC_TestDGEMM(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure);
int MaxMem(int nprocs, int imrow, int imcol, int nmat, int *mval, int *nval, int nbmat,
  int *mbval, int *nbval, int ngrids, int *npval, int *nqval, long *maxMem);
/*
extern int HPCC_Stream(HPCC_Params *params, int doIO, MPI_Comm comm, int world_rank,
  double *copyGBs, double *scaleGBs, double *addGBs, double *triadGBs,
  int *failure);
extern void main_bench_lat_bw(HPCC_Params *params);
*/

/*
extern int pdtrans(char *trans, int *m, int *n, int * mb, int *nb, double *a, int *lda,
  double *beta, double *c__, int *ldc, int *imrow, int *imcol, double *work, int *iwork);
extern FILE* pdtransinfo(int *nmat, int *mval, int *nval, int *ldval,
  int *nbmat, int *mbval, int *nbval, int *ldnbval, int *ngrids, int *npval, int *nqval,
  int *ldpqval, int *iaseed, int *imrow, int *imcol, float *thresh, int *iam, int *nprocs,
  double *eps, char *infname, int *fcl, char *outfname);
int pdmatgen(int *ictxt, char *aform, char *diag, int *m, int *n, int *mb, int *nb, double*a,
  int *lda, int *iarow, int *iacol, int *iseed, int *iroff, int *irnum, int *icoff, int *icnum,
  int * myrow, int *mycol, int *nprow, int *npcol, double alpha);
extern int pdmatcmp(int *ictxt, int *m_, int *n_, double *a, int *lda_, double *aCopy, int *ldc_,
  double *error);
extern int pxerbla(int *ictxt, char *srname, int *info);
extern int slcombine_(int *ictxt, char *scope, char *op, char * timetype, int *n, int *ibeg,
  double *times);
extern int icopy_(int *n, int *sx, int *incx, int * sy, int *incy);
*/
int numroc_(int *, int *, int *, int *, int *);
//extern int slboot_(void);
//extern int sltimer_(int *i__);
int ilcm_(int *, int *);
int iceil_(int *, int *);
//extern double pdrand();
//extern int setran_(int *, int *, int *);
//extern int jumpit_(int *, int *, int *, int *);
//extern int xjumpm_(int *, int *, int *, int *, int *, int *, int *);
/* ---------------------------------------------------------------------- */

#define DPRN(i,v) do{printf(__FILE__ "(%d)@%d:" #v "=%g\n",__LINE__,i,(double)(v));fflush(stdout);}while(0)

#define BEGIN_IO(r,fn,f) if(0==r){f=fopen(fn,"a");if(!f)fprintf(f=stderr,"Problem with appending to file '%s'\n",fn)
#define END_IO(r,f) fflush(f); if (f!=stdout && f!=stderr) fclose(f);} f=(FILE*)(NULL)

#ifndef HPCCMEMA_H
#define HPCCMEMA_H 1

#ifdef HPCC_MEMALLCTR
int HPCC_alloc_init(size_t total_size);
int HPCC_alloc_finalize();
void *HPCC_malloc(size_t size);
void HPCC_free(void *ptr);
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

#endif

#define XMALLOC(t,s) ((t*)malloc(sizeof(t)*(s)))
#define XCALLOC(t,s) ((t*)calloc((s),sizeof(t)))

/* HPL dependant functions */
/** HPL_slamch.c **/

#ifndef FLT_DIGITS
#define FLT_DIGITS 24
#endif

#ifdef HPL_rone
#undef HPL_rone
#endif
#define HPL_rone 1.0f
#ifdef HPL_rtwo
#undef HPL_rtwo
#endif
#define HPL_rtwo 2.0f
#ifdef HPL_rzero
#undef HPL_rzero
#endif
#define HPL_rzero 0.0f
/*
 * ---------------------------------------------------------------------
 * Static function prototypes
 * ---------------------------------------------------------------------
 */
static void     HPL_slamc1
STDC_ARGS(
(  int *,           int *,           int *,           int * ) );
static void     HPL_slamc2
STDC_ARGS(
(  int *,           int *,           int *,           float *,
   int *,           float *,        int *,           float * ) );
static float   HPL_slamc3
STDC_ARGS(
(  const float,    const float ) );
static void     HPL_slamc4
STDC_ARGS(
(  int *,           const float,    const int ) );
static void     HPL_slamc5
STDC_ARGS(
(  const int,       const int,       const int,       const int,
   int *,           float * ) );
static float   HPL_sipow
STDC_ARGS(
(  const float,    const int ) );

#ifdef HPL_STDC_HEADERS
float HPL_slamch
(
   const HPL_T_MACH                 CMACH
)
#else
float HPL_slamch
( CMACH )
   const HPL_T_MACH                 CMACH;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_slamch determines  machine-specific  arithmetic constants such as
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
 * auxiliary function slamch.f  (version 2.0 -- 1992), that  was  itself
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
 *         Specifies the value to be returned by HPL_slamch
 *            = HPL_MACH_EPS,   HPL_slamch := eps (default)
 *            = HPL_MACH_SFMIN, HPL_slamch := sfmin
 *            = HPL_MACH_BASE,  HPL_slamch := base
 *            = HPL_MACH_PREC,  HPL_slamch := eps*base
 *            = HPL_MACH_MLEN,  HPL_slamch := t
 *            = HPL_MACH_RND,   HPL_slamch := rnd
 *            = HPL_MACH_EMIN,  HPL_slamch := emin
 *            = HPL_MACH_RMIN,  HPL_slamch := rmin
 *            = HPL_MACH_EMAX,  HPL_slamch := emax
 *            = HPL_MACH_RMAX,  HPL_slamch := rmax
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
   /*static*/ float              eps, sfmin, base, t, rnd, emin, rmin, emax,
                              rmax, prec;
   float                     small;
   /*static*/ int                 first=0/*1*/;
   int                        beta=0, imax=0, imin=0, it=0, lrnd=0;
/* ..
 * .. Executable Statements ..
 */
   eps = FLT_EPSILON / FLT_RADIX;
   base = FLT_RADIX;
   prec = FLT_EPSILON;
   t = FLT_DIGITS;
   rnd = FLT_ROUNDS < 2 ? HPL_rone : HPL_rzero;
   emin = FLT_MIN_EXP;
   rmin = FLT_MIN;
   emax = FLT_MAX_EXP;
   rmax = FLT_MAX;

   sfmin = rmin;
   small = HPL_rone / rmax;
   if (small >= sfmin)
     sfmin = small * ( HPL_rone + eps );

   if( first != 0 )
   {
      first = 0;
      HPL_slamc2( &beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax );
      base  = (float)(beta);  t     = (float)(it);
      if( lrnd != 0 )
      { rnd = HPL_rone;  eps = HPL_sipow( base, 1 - it ) / HPL_rtwo; }
      else
      { rnd = HPL_rzero; eps = HPL_sipow( base, 1 - it );            }
      prec  = eps * base;  emin  = (float)(imin); emax  = (float)(imax);
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
 * End of HPL_slamch
 */
}

#ifdef HPL_STDC_HEADERS
static void HPL_slamc1
(
   int                        * BETA,
   int                        * T,
   int                        * RND,
   int                        * IEEE1
)
#else
static void HPL_slamc1
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
 * HPL_slamc1  determines  the machine parameters given by BETA, T, RND,
 * and IEEE1.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function slamc1.f  (version 2.0 -- 1992), that  was  itself
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
   float                     a, b, c, f, one, qtr, savec, t1, t2;
   static int                 first=1, lbeta, lieee1, lrnd, lt;
/* ..
 * .. Executable Statements ..
 */
   if( first != 0 )
   {
      first = 0; one = HPL_rone;
/*
 * lbeta, lieee1, lt and lrnd are the local values of BETA, IEEE1, T and
 * RND. Throughout this routine we use the function HPL_slamc3 to ensure
 * that relevant values are stored and not held in registers, or are not
 * affected by optimizers.
 *
 * Compute  a = 2.0**m  with the  smallest  positive integer m such that
 * fl( a + 1.0 ) == a.
 */
      a = HPL_rone; c = HPL_rone;
      do
      { a *= HPL_rtwo; c = HPL_slamc3( a, one ); c = HPL_slamc3( c, -a ); }
      while( c == HPL_rone );
/*
 * Now compute b = 2.0**m with the smallest positive integer m such that
 * fl( a + b ) > a.
 */
      b = HPL_rone; c = HPL_slamc3( a, b );
      while( c == a ) { b *= HPL_rtwo; c = HPL_slamc3( a, b ); }
/*
 * Now compute the base.  a and c  are  neighbouring floating point num-
 * bers in the interval ( BETA**T, BETA**( T + 1 ) ) and so their diffe-
 * rence is BETA.  Adding 0.25 to c is to ensure that it is truncated to
 * BETA and not (BETA-1).
 */
      qtr = one / 4.0; savec = c;
      c   = HPL_slamc3( c, -a ); lbeta = (int)(c+qtr);
/*
 * Now  determine  whether  rounding or chopping occurs, by adding a bit
 * less than BETA/2 and a bit more than BETA/2 to a.
 */
      b = (float)(lbeta);
      f = HPL_slamc3( b / HPL_rtwo, -b / 100.0 ); c = HPL_slamc3( f, a );
      if( c == a ) { lrnd = 1; } else { lrnd = 0; }
      f = HPL_slamc3( b / HPL_rtwo,  b / 100.0 ); c = HPL_slamc3( f, a );
      if( ( lrnd != 0 ) && ( c == a ) ) lrnd = 0;
/*
 * Try  and decide whether rounding is done in the  IEEE  round to nea-
 * rest style.  b/2 is half a unit in the last place of the two numbers
 * a  and savec. Furthermore, a is even, i.e. has last bit zero, and sa-
 * vec is odd.  Thus adding b/2 to a should not change a, but adding b/2
 * to savec should change savec.
 */
      t1 = HPL_slamc3( b / HPL_rtwo, a );
      t2 = HPL_slamc3( b / HPL_rtwo, savec );
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
         lt++; a *= (float)(lbeta);
         c = HPL_slamc3( a, one ); c = HPL_slamc3( c,  -a );
      } while( c == HPL_rone );
   }

   *BETA  = lbeta; *T = lt; *RND = lrnd; *IEEE1 = lieee1;
}

#ifdef HPL_STDC_HEADERS
static void HPL_slamc2
(
   int                        * BETA,
   int                        * T,
   int                        * RND,
   float                     * EPS,
   int                        * EMIN,
   float                     * RMIN,
   int                        * EMAX,
   float                     * RMAX
)
#else
static void HPL_slamc2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
/*
 * .. Scalar Arguments ..
 */
   int                        * BETA, * EMAX, * EMIN, * RND, * T;
   float                     * EPS, * RMAX, * RMIN;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_slamc2  determines the machine  parameters specified in its argu-
 * ment list.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function  slamc2.f (version 2.0 -- 1992), that  was  itself
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
 * EPS     (local output)              float *
 *         The smallest positive number such that fl( 1.0 - EPS ) < 1.0,
 *         where fl denotes the computed value.
 *
 * EMIN    (local output)              int *
 *         The minimum exponent before (gradual) underflow occurs.
 *
 * RMIN    (local output)              float *
 *         The smallest  normalized  number  for  the  machine, given by
 *         BASE**( EMIN - 1 ), where  BASE  is the floating  point value
 *         of BETA.
 *
 * EMAX    (local output)              int *
 *         The maximum exponent before overflow occurs.
 *
 * RMAX    (local output)              float *
 *         The  largest  positive  number  for  the  machine,  given  by
 *         BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating  point
 *         value of BETA.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   static float              leps, lrmax, lrmin;
   float                     a, b, c, half, one, rbase, sixth, small,
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
 * Throughout this routine we use the function HPL_slamc3 to ensure that
 * relevant values are stored and not held in registers,  or are not af-
 * fected by optimizers.
 *
 * HPL_slamc1 returns the parameters  lbeta, lt, lrnd and lieee1.
 */
      HPL_slamc1( &lbeta, &lt, &lrnd, &lieee1 );
/*
 * Start to find eps.
 */
      b = (float)(lbeta); a = HPL_sipow( b, -lt ); leps = a;
/*
 * Try some tricks to see whether or not this is the correct  EPS.
 */
      b     = two / 3.0;
      half  = one / HPL_rtwo;
      sixth = HPL_slamc3( b, -half );
      third = HPL_slamc3( sixth, sixth );
      b     = HPL_slamc3( third, -half );
      b     = HPL_slamc3( b, sixth );
      b     = Mabs( b ); if( b < leps ) b = leps;

      leps = HPL_rone;

      while( ( leps > b ) && ( b > zero ) )
      {
         leps = b;
         c = HPL_slamc3( half * leps,
                         HPL_sipow( two, 5 ) * HPL_sipow( leps, 2 ) );
         c = HPL_slamc3( half, -c ); b = HPL_slamc3( half, c );
         c = HPL_slamc3( half, -b ); b = HPL_slamc3( half, c );
      }
      if( a < leps ) leps = a;
/*
 * Computation of EPS complete.
 *
 * Now find  EMIN.  Let a = + or - 1, and + or - (1 + BASE**(-3)).  Keep
 * dividing a by BETA until (gradual) underflow occurs. This is detected
 * when we cannot recover the previous a.
 */
      rbase = one / (float)(lbeta); small = one;
      for( i = 0; i < 3; i++ ) small = HPL_slamc3( small * rbase, zero );
      a = HPL_slamc3( one, small );
      HPL_slamc4( &ngpmin, one, lbeta ); HPL_slamc4( &ngnmin, -one, lbeta );
      HPL_slamc4( &gpmin,    a, lbeta ); HPL_slamc4( &gnmin,    -a, lbeta );

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
"out the  if  block  as marked within the code of routine  HPL_slamc2, ",
"otherwise supply EMIN explicitly." );
      }
/*
 * Assume IEEE arithmetic if we found denormalised  numbers above, or if
 * arithmetic seems to round in the  IEEE style,  determined  in routine
 * HPL_slamc1.  A true  IEEE  machine should have both things true; how-
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
         lrmin = HPL_slamc3( lrmin*rbase, zero );
/*
 * Finally, call HPL_slamc5 to compute emax and rmax.
 */
      HPL_slamc5( lbeta, lt, lemin, ieee, &lemax, &lrmax );
   }
   *BETA = lbeta; *T    = lt;    *RND  = lrnd;  *EPS  = leps;
   *EMIN = lemin; *RMIN = lrmin; *EMAX = lemax; *RMAX = lrmax;
}

#ifdef HPL_STDC_HEADERS
static float HPL_slamc3( const float A, const float B )
#else
static float HPL_slamc3( A, B )
/*
 * .. Scalar Arguments ..
 */
   const float               A, B;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_slamc3  is intended to force a and b  to be stored prior to doing
 * the addition of  a  and  b,  for  use  in situations where optimizers
 * might hold one of these in a register.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function slamc3.f (version 2.0 -- 1992).
 *
 * Arguments
 * =========
 *
 * A, B    (local input)               float
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
static void HPL_slamc4
(
   int                        * EMIN,
   const float               START,
   const int                  BASE
)
#else
static void HPL_slamc4( EMIN, START, BASE )
/*
 * .. Scalar Arguments ..
 */
   int                        * EMIN;
   const int                  BASE;
   const float               START;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_slamc4 is a service function for HPL_slamc2.
 *
 * Notes
 * =====
 *
 * This function has been manually translated from the Fortran 77 LAPACK
 * auxiliary function slamc4.f (version 2.0 -- 1992).
 *
 * Arguments
 * =========
 *
 * EMIN    (local output)              int *
 *         The minimum exponent before  (gradual) underflow, computed by
 *         setting A = START and dividing  by  BASE until the previous A
 *         can not be recovered.
 *
 * START   (local input)               float
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
   float                     a, b1, b2, c1, c2, d1, d2, one, rbase, zero;
   int                        i;
/* ..
 * .. Executable Statements ..
 */
   a     = START; one = HPL_rone; rbase = one / (float)(BASE);
   zero  = HPL_rzero;
   *EMIN = 1; b1 = HPL_slamc3( a * rbase, zero ); c1 = c2 = d1 = d2 = a;

   do
   {
      (*EMIN)--; a = b1;
      b1 = HPL_slamc3( a /  BASE,  zero );
      c1 = HPL_slamc3( b1 *  BASE, zero );
      d1 = zero; for( i = 0; i < BASE; i++ ) d1 = d1 + b1;
      b2 = HPL_slamc3( a * rbase,  zero );
      c2 = HPL_slamc3( b2 / rbase, zero );
      d2 = zero; for( i = 0; i < BASE; i++ ) d2 = d2 + b2;
   } while( ( c1 == a ) && ( c2 == a ) &&  ( d1 == a ) && ( d2 == a ) );
}

#ifdef HPL_STDC_HEADERS
static void HPL_slamc5
(
   const int                  BETA,
   const int                  P,
   const int                  EMIN,
   const int                  IEEE,
   int                        * EMAX,
   float                     * RMAX
)
#else
static void HPL_slamc5( BETA, P, EMIN, IEEE, EMAX, RMAX )
/*
 * .. Scalar Arguments ..
 */
   const int                  BETA, EMIN, IEEE, P;
   int                        * EMAX;
   float                     * RMAX;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_slamc5  attempts  to compute RMAX, the largest machine  floating-
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
 * auxiliary function slamc5.f (version 2.0 -- 1992).
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
 * RMAX    (local output)              float *
 *         The largest machine floating-point number.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   float                     oldy=HPL_rzero, recbas, y, z;
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
   recbas = HPL_rone / (float)(BETA);
   z      = (float)(BETA) - HPL_rone;
   y      = HPL_rzero;

   for( i = 0; i < P; i++ )
   { z *= recbas; if( y < HPL_rone ) oldy = y; y = HPL_slamc3( y, z ); }

   if( y >= HPL_rone ) y = oldy;
/*
 * Now multiply by BETA**EMAX to get RMAX.
 */
   for( i = 0; i < *EMAX; i++ ) y = HPL_slamc3( y * BETA, HPL_rzero );

   *RMAX = y;
/*
 * End of HPL_slamch
 */
}

#ifdef HPL_STDC_HEADERS
static float HPL_sipow
(
   const float               X,
   const int                  N
)
#else
static float HPL_sipow( X, N )
/*
 * .. Scalar Arguments ..
 */
   const int                  N;
   const float               X;
#endif
{
/*
 * Purpose
 * =======
 *
 * HPL_sipow computes the integer n-th power of a real scalar x.
 *
 * Arguments
 * =========
 *
 * X       (local input)               const float
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
   float                     r, y=HPL_rone;
   int                        k, n;
/* ..
 * .. Executable Statements ..
 */
   if( X == HPL_rzero ) return( HPL_rzero );
   if( N < 0 ) { n = -N; r = HPL_rone / X; } else { n = N; r = X; }
   for( k = 0; k < n; k++ ) y *= r;

   return( y );
}

/** HPL_pabort.c **/
#ifdef HPL_STDC_HEADERS
void HPL_pabort
(
   int                              LINE,
   const char *                     SRNAME,
   const char *                     FORM,
   ...                              
)
#else
void HPL_pabort( va_alist )
va_dcl
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pabort displays an error message on stderr and halts execution.
 * 
 *
 * Arguments
 * =========
 *
 * LINE    (local input)                 int
 *         On entry,  LINE  specifies the line  number in the file where
 *         the  error  has  occured.  When  LINE  is not a positive line
 *         number, it is ignored.
 *
 * SRNAME  (local input)                 const char *
 *         On entry, SRNAME  should  be the name of the routine  calling
 *         this error handler.
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
   int                        rank;
   char                       cline[128];
#ifndef HPL_STDC_HEADERS
   int                        LINE;
   char                       * FORM, * SRNAME;
#endif
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_STDC_HEADERS
   va_start( argptr, FORM );
#else
   va_start( argptr );
   LINE   = va_arg( argptr, int      );
   SRNAME = va_arg( argptr, char *   );
   FORM   = va_arg( argptr, char *   );
#endif
   (void) vsprintf( cline, FORM, argptr );
   va_end( argptr ); 

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
/*
 * Display an error message
 */
   if( LINE <= 0 )
      HPL_fprintf( stderr, "%s %s %d, %s %s:\n>>> %s <<< Abort ...\n\n",
                   "HPL ERROR", "from process #", rank, "in function",
                   SRNAME, cline );
   else
      HPL_fprintf( stderr,
                   "%s %s %d, %s %d %s %s:\n>>> %s <<< Abort ...\n\n",
                   "HPL ERROR", "from process #", rank, "on line", LINE,
                   "of function", SRNAME, cline );

   MPI_Abort( MPI_COMM_WORLD, -1 );
   exit( -1 );
/*
 * End of HPL_pabort
 */
}
/** HPL_reduce.c **/
#ifdef HPL_STDC_HEADERS
int HPL_reduce
(
   void *                           BUFFER,
   const int                        COUNT,
   const HPL_T_TYPE                 DTYPE,
   const HPL_T_OP                   OP,
   const int                        ROOT,
   MPI_Comm                         COMM
)
#else
int HPL_reduce
( BUFFER, COUNT, DTYPE, OP, ROOT, COMM )
   void *                           BUFFER;
   const int                        COUNT;
   const HPL_T_TYPE                 DTYPE;
   const HPL_T_OP                   OP;
   const int                        ROOT;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_reduce performs a global reduce operation across all processes of
 * a group.  Note that the input buffer is  used as workarray and in all
 * processes but the accumulating process corrupting the original data.
 *
 * Arguments
 * =========
 *
 * BUFFER  (local input/output)          void *
 *         On entry,  BUFFER  points to  the  buffer to be  reduced.  On
 *         exit,  and  in process of rank  ROOT  this array contains the
 *         reduced data.  This  buffer  is also used as workspace during
 *         the operation in the other processes of the group.
 *
 * COUNT   (global input)                const int
 *         On entry,  COUNT  indicates the number of entries in  BUFFER.
 *         COUNT must be at least zero.
 *
 * DTYPE   (global input)                const HPL_T_TYPE
 *         On entry,  DTYPE  specifies the type of the buffers operands.
 *
 * OP      (global input)                const HPL_T_OP 
 *         On entry, OP is a pointer to the local combine function.
 *
 * ROOT    (global input)                const int
 *         On entry, ROOT is the coordinate of the accumulating process.
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   MPI_Status                 status;
   void                       * buffer = NULL;
   int                        hplerr=MPI_SUCCESS, d=1, i, ip2=1, mask=0,
                              mpierr, mydist, partner, rank, size, 
                              tag = MSGID_BEGIN_COLL;
/* ..
 * .. Executable Statements ..
 */
   if( COUNT <= 0 ) return( MPI_SUCCESS );
   mpierr = MPI_Comm_size( COMM, &size );
   if( size  == 1 ) return( MPI_SUCCESS );
   mpierr = MPI_Comm_rank( COMM, &rank );
   i = size - 1; while( i > 1 ) { i >>= 1; d++; }

   if( DTYPE == HPL_INT )
      buffer = (void *)( (int *)   malloc( (size_t)(COUNT) * 
                                           sizeof( int    ) ) );
   else
      buffer = (void *)( (double *)malloc( (size_t)(COUNT) *
                                           sizeof( double ) ) );

   if( !( buffer ) )
   { HPL_pabort( __LINE__, "HPL_reduce", "Memory allocation failed" ); }

   if( ( mydist = MModSub( rank, ROOT, size ) ) == 0 )
   {
      do
      {
         mpierr = MPI_Recv( buffer, COUNT, HPL_2_MPI_TYPE( DTYPE ),
                            MModAdd( ROOT, ip2, size ), tag, COMM,
                            &status );
         if( mpierr != MPI_SUCCESS ) hplerr = mpierr;
         OP( COUNT, buffer, BUFFER, DTYPE );
         ip2 <<= 1; d--;
      } while( d );
   }
   else
   {
      do
      {
         if( ( mydist & mask ) == 0 )
         {
            partner = mydist ^ ip2;

            if( mydist & ip2 )
            {
               partner = MModAdd( ROOT, partner, size );
               mpierr = MPI_Send( BUFFER, COUNT, HPL_2_MPI_TYPE( DTYPE ),
                                  partner, tag, COMM );
            }
            else if( partner < size )
            {
               partner = MModAdd( ROOT, partner, size );
               mpierr  = MPI_Recv( buffer, COUNT, HPL_2_MPI_TYPE( DTYPE ),
                                   partner, tag, COMM, &status );
               OP( COUNT, buffer, BUFFER, DTYPE );
            }
            if( mpierr != MPI_SUCCESS ) hplerr = mpierr;
         }
         mask ^= ip2; ip2 <<= 1; d--;
      } while( d );
   }
   if( buffer ) free( buffer );

   return( hplerr );
/*
 * End of HPL_reduce
 */
}

/** HPL_min.c **/
#ifdef HPL_STDC_HEADERS
void HPL_min
(
   const int                        N,
   const void *                     IN,
   void *                           INOUT,
   const HPL_T_TYPE                 DTYPE
)
#else
void HPL_min
( N, IN, INOUT, DTYPE )
   const int                        N;
   const void *                     IN;
   void *                           INOUT;
   const HPL_T_TYPE                 DTYPE;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_min combines (min) two buffers.
 * 
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N  specifies  the  length  of  the  buffers  to  be
 *         combined. N must be at least zero.
 *
 * IN      (input)                       const void *
 *         On entry, IN points to the input-only buffer to be combined.
 *
 * INOUT   (input/output)                void *
 *         On entry, INOUT  points  to  the  input-output  buffer  to be
 *         combined.  On exit,  the  entries of this array contains  the
 *         combined results.
 *
 * DTYPE   (input)                       const HPL_T_TYPE
 *         On entry,  DTYPE  specifies the type of the buffers operands.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   register int               i;
/* ..
 * .. Executable Statements ..
 */
   if( DTYPE == HPL_INT )
   {
      const int       * a = (const int *)(IN);
      int             * b = (int *)(INOUT);
      for( i = 0; i < N; i++ ) b[i] = Mmin( a[i], b[i] );
   }
   else
   {
      const double    * a = (const double *)(IN);
      double          * b = (double *)(INOUT);
      for( i = 0; i < N; i++ ) b[i] = Mmin( a[i], b[i] );
   }
/*
 * End of HPL_min
 */
}

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

/** HPL_broadcast.c **/
#ifdef HPL_STDC_HEADERS
int HPL_broadcast
(
   void *                           BUFFER,
   const int                        COUNT,
   const HPL_T_TYPE                 DTYPE,
   const int                        ROOT,
   MPI_Comm                         COMM
)
#else
int HPL_broadcast
( BUFFER, COUNT, DTYPE, ROOT, COMM )
   void *                           BUFFER;
   const int                        COUNT;
   const HPL_T_TYPE                 DTYPE;
   const int                        ROOT;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_broadcast broadcasts  a message from the process with rank ROOT to
 * all processes in the group.
 *
 * Arguments
 * =========
 *
 * BUFFER  (local input/output)          void *
 *         On entry,  BUFFER  points to  the  buffer to be broadcast. On
 *         exit, this array contains the broadcast data and is identical
 *         on all processes in the group.
 *
 * COUNT   (global input)                const int
 *         On entry,  COUNT  indicates the number of entries in  BUFFER.
 *         COUNT must be at least zero.
 *
 * DTYPE   (global input)                const HPL_T_TYPE
 *         On entry,  DTYPE  specifies the type of the buffers operands.
 *
 * ROOT    (global input)                const int
 *         On entry, ROOT is the coordinate of the source process.
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   int                        hplerr=MPI_SUCCESS, ip2=1, kk, mask=1, 
                              mpierr, mydist, partner, rank, size, 
                              tag = MSGID_BEGIN_COLL;
   MPI_Status                 status;
/* ..
 * .. Executable Statements ..
 */
   if( COUNT <= 0 ) return( MPI_SUCCESS );
   mpierr = MPI_Comm_size( COMM, &size ); if( size <= 1 ) return( mpierr );
   mpierr = MPI_Comm_rank( COMM, &rank );

   kk = size - 1;
   while( kk > 1 ) { kk >>= 1; ip2 <<= 1; mask <<= 1; mask++; }
   mydist = MModSub( rank, ROOT, size );

   do
   {
      mask ^= ip2;
      if( ( mydist & mask ) == 0 )
      {
         partner = mydist ^ ip2;

         if( mydist & ip2 )
         {
            partner = MModAdd( ROOT, partner, size );
            mpierr  = MPI_Recv(  BUFFER, COUNT, HPL_2_MPI_TYPE( DTYPE ),
                                 partner, tag, COMM, &status );
         }
         else if( partner < size )
         {
            partner = MModAdd( ROOT, partner, size );
            mpierr  = MPI_Send( BUFFER, COUNT, HPL_2_MPI_TYPE( DTYPE ),
                                partner, tag, COMM );
         }
         if( mpierr != MPI_SUCCESS ) hplerr = mpierr;
      }
      ip2 >>= 1;
   } while( ip2 );

   return( hplerr );
/*
 * End of HPL_broadcast
 */
}

/** HPL_pdlamch.c **/
#ifdef HPL_STDC_HEADERS
double HPL_pdlamch
(
   MPI_Comm                         COMM,
   const HPL_T_MACH                 CMACH
)
#else
double HPL_pdlamch
( COMM, CMACH )
   MPI_Comm                         COMM;
   const HPL_T_MACH                 CMACH;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdlamch determines  machine-specific  arithmetic  constants  such  as
 * the relative machine precision (eps),  the safe minimum(sfmin) such that
 * 1/sfmin does not overflow, the base of the machine (base), the precision
 * (prec),  the  number  of  (base)  digits in the  mantissa  (t),  whether
 * rounding occurs in addition (rnd = 1.0 and 0.0 otherwise),  the  minimum
 * exponent before  (gradual)  underflow (emin),  the  underflow  threshold
 * (rmin)- base**(emin-1), the largest exponent before overflow (emax), the
 * overflow threshold (rmax)  - (base**emax)*(1-eps).
 *
 * Arguments
 * =========
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * CMACH   (global input)                const HPL_T_MACH
 *         Specifies the value to be returned by HPL_pdlamch            
 *            = HPL_MACH_EPS,   HPL_pdlamch := eps (default)            
 *            = HPL_MACH_SFMIN, HPL_pdlamch := sfmin                    
 *            = HPL_MACH_BASE,  HPL_pdlamch := base                     
 *            = HPL_MACH_PREC,  HPL_pdlamch := eps*base                 
 *            = HPL_MACH_MLEN,  HPL_pdlamch := t                        
 *            = HPL_MACH_RND,   HPL_pdlamch := rnd                      
 *            = HPL_MACH_EMIN,  HPL_pdlamch := emin                     
 *            = HPL_MACH_RMIN,  HPL_pdlamch := rmin                     
 *            = HPL_MACH_EMAX,  HPL_pdlamch := emax                     
 *            = HPL_MACH_RMAX,  HPL_pdlamch := rmax                     
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
   double                     param;
/* ..
 * .. Executable Statements ..
 */
   param = HPL_dlamch( CMACH );

   switch( CMACH )
   {
      case HPL_MACH_EPS   :
      case HPL_MACH_SFMIN :
      case HPL_MACH_EMIN  :
      case HPL_MACH_RMIN  :
         (void) HPL_all_reduce( (void *)(&param), 1, HPL_DOUBLE,
                                HPL_max, COMM );
         break;
      case HPL_MACH_EMAX  :
      case HPL_MACH_RMAX  :
         (void) HPL_all_reduce( (void *)(&param), 1, HPL_DOUBLE,
                                HPL_min, COMM );
         break;
      default             :
         break;
   } 

   return( param );
/*
 * End of HPL_pdlamch
 */
}

/** HPL_all_reduce.c **/

#ifdef HPL_STDC_HEADERS
int HPL_all_reduce
(
   void *                           BUFFER,
   const int                        COUNT,
   const HPL_T_TYPE                 DTYPE,
   const HPL_T_OP                   OP,
   MPI_Comm                         COMM
)
#else
int HPL_all_reduce
( BUFFER, COUNT, DTYPE, OP, COMM )
   void *                           BUFFER;
   const int                        COUNT;
   const HPL_T_TYPE                 DTYPE;
   const HPL_T_OP                   OP;
   MPI_Comm                         COMM;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_all_reduce performs   a   global   reduce  operation  across  all
 * processes of a group leaving the results on all processes.
 *
 * Arguments
 * =========
 *
 * BUFFER  (local input/global output)   void *
 *         On entry,  BUFFER  points to  the  buffer to be combined.  On
 *         exit, this array contains the combined data and  is identical
 *         on all processes in the group.
 *
 * COUNT   (global input)                const int
 *         On entry,  COUNT  indicates the number of entries in  BUFFER.
 *         COUNT must be at least zero.
 *
 * DTYPE   (global input)                const HPL_T_TYPE
 *         On entry,  DTYPE  specifies the type of the buffers operands.
 *
 * OP      (global input)                const HPL_T_OP 
 *         On entry, OP is a pointer to the local combine function.
 *
 * COMM    (global/local input)          MPI_Comm
 *         The MPI communicator identifying the process collection.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   int                        hplerr;
/* ..
 * .. Executable Statements ..
 */
   hplerr = HPL_reduce(   BUFFER, COUNT, DTYPE, OP, 0, COMM );
   if( hplerr != MPI_SUCCESS ) return( hplerr );
   return( HPL_broadcast( BUFFER, COUNT, DTYPE,     0, COMM ) );
/*
 * End of HPL_all_reduce
 */
}

/** HPL_max.c **/
#ifdef HPL_STDC_HEADERS
void HPL_max
(
   const int                        N,
   const void *                     IN,
   void *                           INOUT,
   const HPL_T_TYPE                 DTYPE
)
#else
void HPL_max
( N, IN, INOUT, DTYPE )
   const int                        N;
   const void *                     IN;
   void *                           INOUT;
   const HPL_T_TYPE                 DTYPE;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_max combines (max) two buffers.
 * 
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N  specifies  the  length  of  the  buffers  to  be
 *         combined. N must be at least zero.
 *
 * IN      (input)                       const void *
 *         On entry, IN points to the input-only buffer to be combined.
 *
 * INOUT   (input/output)                void *
 *         On entry, INOUT  points  to  the  input-output  buffer  to be
 *         combined.  On exit,  the  entries of this array contains  the
 *         combined results.
 *
 * DTYPE   (input)                       const HPL_T_TYPE
 *         On entry,  DTYPE  specifies the type of the buffers operands.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   register int               i;
/* ..
 * .. Executable Statements ..
 */
   if( DTYPE == HPL_INT )
   {
      const int       * a = (const int *)(IN);
      int             * b = (int *)(INOUT);
      for( i = 0; i < N; i++ ) b[i] = Mmax( a[i], b[i] );
   }
   else
   {
      const double    * a = (const double *)(IN);
      double          * b = (double *)(INOUT);
      for( i = 0; i < N; i++ ) b[i] = Mmax( a[i], b[i] );
   }
/*
 * End of HPL_max
 */
}

/** HPL_warn.c **/
#ifdef HPL_STDC_HEADERS
void HPL_pwarn
(
   FILE *                           STREAM,
   int                              LINE,
   const char *                     SRNAME,
   const char *                     FORM,
   ...                              
)
#else
void HPL_pwarn( va_alist )
va_dcl
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pwarn displays an error message.
 * 
 *
 * Arguments
 * =========
 *
 * STREAM  (local input)                 FILE *
 *         On entry, STREAM specifies the output stream.
 *
 * LINE    (local input)                 int
 *         On entry,  LINE  specifies the line  number in the file where
 *         the  error  has  occured.  When  LINE  is not a positive line
 *         number, it is ignored.
 *
 * SRNAME  (local input)                 const char *
 *         On entry, SRNAME  should  be the name of the routine  calling
 *         this error handler.
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
   int                        rank;
   char                       cline[128];
#ifndef HPL_STDC_HEADERS
   FILE                       * STREAM;
   int                        LINE;
   char                       * FORM, * SRNAME;
#endif
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_STDC_HEADERS
   va_start( argptr, FORM );
#else
   va_start( argptr );
   STREAM = va_arg( argptr, FILE * );
   LINE   = va_arg( argptr, int    );
   SRNAME = va_arg( argptr, char * );
   FORM   = va_arg( argptr, char * );
#endif
   (void) vsprintf( cline, FORM, argptr );
   va_end( argptr ); 

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
/*
 * Display an error message
 */
   if( LINE <= 0 )
      HPL_fprintf( STREAM, "%s %s %d, %s %s:\n>>> %s <<<\n\n",
                   "HPL ERROR", "from process #", rank, "in function",
                   SRNAME, cline );
   else if( LINE >  (1 << 30) )
      HPL_fprintf( STREAM, "%s %s %d, %s %d %s %s:\n>>> %s <<<\n\n",
                   "HPL WARNING", "from process #", rank, "on line", LINE - (1 << 30),
                   "of function", SRNAME, cline );
   else
      HPL_fprintf( STREAM, "%s %s %d, %s %d %s %s:\n>>> %s <<<\n\n",
                   "HPL ERROR", "from process #", rank, "on line", LINE,
                   "of function", SRNAME, cline );
/*
 * End of HPL_pwarn
 */
}


/** HPL_pdinfo.c **/
#ifdef HPL_STDC_HEADERS
void HPL_pdinfo
(
   HPL_T_test *                     TEST,
   int *                            NS,
   int *                            N,
   int *                            NBS,
   int *                            NB,
   HPL_T_ORDER *                    PMAPPIN,
   int *                            NPQS,
   int *                            P,
   int *                            Q,
   int *                            NPFS,
   HPL_T_FACT *                     PF,
   int *                            NBMS,
   int *                            NBM,
   int *                            NDVS,
   int *                            NDV,
   int *                            NRFS,
   HPL_T_FACT *                     RF,
   int *                            NTPS,
   HPL_T_TOP *                      TP,
   int *                            NDHS,
   int *                            DH,
   HPL_T_SWAP *                     FSWAP,
   int *                            TSWAP,
   int *                            L1NOTRAN,
   int *                            UNOTRAN,
   int *                            EQUIL,
   int *                            ALIGN
)
#else
void HPL_pdinfo
( TEST, NS, N, NBS, NB, PMAPPIN, NPQS, P, Q, NPFS, PF, NBMS, NBM, NDVS, NDV, NRFS, RF, NTPS, TP, NDHS, DH, FSWAP, TSWAP, L1NOTRAN, UNOTRAN, EQUIL, ALIGN )
   HPL_T_test *                     TEST;
   int *                            NS;
   int *                            N;
   int *                            NBS;
   int *                            NB;
   HPL_T_ORDER *                    PMAPPIN;
   int *                            NPQS;
   int *                            P;
   int *                            Q;
   int *                            NPFS;
   HPL_T_FACT *                     PF;
   int *                            NBMS;
   int *                            NBM;
   int *                            NDVS;
   int *                            NDV;
   int *                            NRFS;
   HPL_T_FACT *                     RF;
   int *                            NTPS;
   HPL_T_TOP *                      TP;
   int *                            NDHS;
   int *                            DH;
   HPL_T_SWAP *                     FSWAP;
   int *                            TSWAP;
   int *                            L1NOTRAN;
   int *                            UNOTRAN;
   int *                            EQUIL;
   int *                            ALIGN;
#endif
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdinfo reads  the  startup  information for the various tests and
 * transmits it to all processes.
 *
 * Arguments
 * =========
 *
 * TEST    (global output)               HPL_T_test *
 *         On entry, TEST  points to a testing data structure.  On exit,
 *         the fields of this data structure are initialized as follows:
 *         TEST->outfp  specifies the output file where the results will
 *         be printed.  It is only defined and used by  the process 0 of
 *         the grid.  TEST->thrsh specifies the threshhold value for the
 *         test ratio.  TEST->epsil is the relative machine precision of
 *         the distributed computer.  Finally  the test counters, kfail,
 *         kpass, kskip, ktest are initialized to zero.
 *
 * NS      (global output)               int *
 *         On exit,  NS  specifies the number of different problem sizes
 *         to be tested. NS is less than or equal to HPL_MAX_PARAM.
 *
 * N       (global output)               int *
 *         On entry, N is an array of dimension HPL_MAX_PARAM.  On exit,
 *         the first NS entries of this array contain the  problem sizes
 *         to run the code with.
 *
 * NBS     (global output)               int *
 *         On exit,  NBS  specifies the number of different distribution
 *         blocking factors to be tested. NBS must be less than or equal
 *         to HPL_MAX_PARAM.
 *
 * NB      (global output)               int *
 *         On exit,  PMAPPIN  specifies the process mapping onto the no-
 *         des of the  MPI machine configuration.  PMAPPIN  defaults  to
 *         row-major ordering.
 *
 * PMAPPIN (global output)               HPL_T_ORDER *
 *         On entry, NB is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first NBS entries of this array contain the values of the
 *         various distribution blocking factors, to run the code with.
 *
 * NPQS    (global output)               int *
 *         On exit, NPQS  specifies the  number of different values that
 *         can be used for P and Q, i.e., the number of process grids to
 *         run  the  code with.  NPQS must be  less  than  or  equal  to
 *         HPL_MAX_PARAM.
 *
 * P       (global output)               int *
 *         On entry, P  is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first NPQS entries of this array contain the values of P,
 *         the number of process rows of the  NPQS grids to run the code
 *         with.
 *
 * Q       (global output)               int *
 *         On entry, Q  is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first NPQS entries of this array contain the values of Q,
 *         the number of process columns of the  NPQS  grids to  run the
 *         code with.
 *
 * NPFS    (global output)               int *
 *         On exit, NPFS  specifies the  number of different values that
 *         can be used for PF : the panel factorization algorithm to run
 *         the code with. NPFS is less than or equal to HPL_MAX_PARAM.
 *
 * PF      (global output)               HPL_T_FACT *
 *         On entry, PF is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first  NPFS  entries  of this array  contain  the various
 *         panel factorization algorithms to run the code with.
 *
 * NBMS    (global output)               int *
 *         On exit,  NBMS  specifies  the  number  of  various recursive
 *         stopping criteria  to be tested.  NBMS  must be  less than or
 *         equal to HPL_MAX_PARAM.
 *
 * NBM     (global output)               int *
 *         On entry,  NBM  is an array of  dimension  HPL_MAX_PARAM.  On
 *         exit, the first NBMS entries of this array contain the values
 *         of the various recursive stopping criteria to be tested.
 *
 * NDVS    (global output)               int *
 *         On exit,  NDVS  specifies  the number  of various numbers  of
 *         panels in recursion to be tested.  NDVS is less than or equal
 *         to HPL_MAX_PARAM.
 *
 * NDV     (global output)               int *
 *         On entry,  NDV  is an array of  dimension  HPL_MAX_PARAM.  On
 *         exit, the first NDVS entries of this array contain the values
 *         of the various numbers of panels in recursion to be tested.
 *
 * NRFS    (global output)               int *
 *         On exit, NRFS  specifies the  number of different values that
 *         can be used for RF : the recursive factorization algorithm to
 *         be tested. NRFS is less than or equal to HPL_MAX_PARAM.
 *
 * RF      (global output)               HPL_T_FACT *
 *         On entry, RF is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first  NRFS  entries  of  this array contain  the various
 *         recursive factorization algorithms to run the code with.
 *
 * NTPS    (global output)               int *
 *         On exit, NTPS  specifies the  number of different values that
 *         can be used for the  broadcast topologies  to be tested. NTPS
 *         is less than or equal to HPL_MAX_PARAM.
 *
 * TP      (global output)               HPL_T_TOP *
 *         On entry, TP is an array of dimension HPL_MAX_PARAM. On exit,
 *         the  first NTPS  entries of this  array  contain  the various
 *         broadcast (along rows) topologies to run the code with.
 *
 * NDHS    (global output)               int *
 *         On exit, NDHS  specifies the  number of different values that
 *         can be used for the  lookahead depths to be  tested.  NDHS is
 *         less than or equal to HPL_MAX_PARAM.
 *
 * DH      (global output)               int *
 *         On entry,  DH  is  an array of  dimension  HPL_MAX_PARAM.  On
 *         exit, the first NDHS entries of this array contain the values
 *         of lookahead depths to run the code with.  Such a value is at
 *         least 0 (no-lookahead) or greater than zero.
 *
 * FSWAP   (global output)               HPL_T_SWAP *
 *         On exit, FSWAP specifies the swapping algorithm to be used in
 *         all tests.
 *
 * TSWAP   (global output)               int *
 *         On exit,  TSWAP  specifies the swapping threshold as a number
 *         of columns when the mixed swapping algorithm was chosen.
 *
 * L1NOTRA (global output)               int *
 *         On exit, L1NOTRAN specifies whether the upper triangle of the
 *         panels of columns  should  be stored  in  no-transposed  form
 *         (L1NOTRAN=1) or in transposed form (L1NOTRAN=0).
 *
 * UNOTRAN (global output)               int *
 *         On exit, UNOTRAN  specifies whether the panels of rows should
 *         be stored in  no-transposed form  (UNOTRAN=1)  or  transposed
 *         form (UNOTRAN=0) during their broadcast.
 *
 * EQUIL   (global output)               int *
 *         On exit,  EQUIL  specifies  whether  equilibration during the
 *         swap-broadcast  of  the  panel of rows  should  be  performed
 *         (EQUIL=1) or not (EQUIL=0).
 *
 * ALIGN   (global output)               int *
 *         On exit,  ALIGN  specifies the alignment  of  the dynamically
 *         allocated buffers in double precision words. ALIGN is greater
 *         than zero.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   char                       file[HPL_LINE_MAX], line[HPL_LINE_MAX],
                              auth[HPL_LINE_MAX], num [HPL_LINE_MAX];
   FILE                       * infp;
   int                        * iwork;
   char                       * lineptr;
   int                        error=0, fid, i, j, lwork, maxp, nprocs,
                              rank, size;
/* ..
 * .. Executable Statements ..
 */
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &size );
/*
 * Initialize the TEST data structure with default values
 */
   TEST->outfp = stderr; TEST->epsil = 2.0e-16; TEST->thrsh = 16.0;
   TEST->kfail = TEST->kpass = TEST->kskip = TEST->ktest = 0;
/*
 * Process 0 reads the input data, broadcasts to other processes and
 * writes needed information to TEST->outfp.
 */
   if( rank == 0 )
   {
/*
 * Open file and skip data file header
 */
#define INFILE "hpccinf.txt"
      if( ( infp = fopen( INFILE, "r" ) ) == NULL )
      { 
         HPL_pwarn( stderr, __LINE__ + (1 << 30), "HPL_pdinfo",
                    "cannot open file " INFILE );
         error = 1; /* goto label_error; */
      }

      if (infp) {
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) fgets( auth, HPL_LINE_MAX - 2, infp );
/*
 * Read name and unit number for summary output file
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", file );
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num  );
      fid = atoi( num );
      }

      fid = 8; /* always write to a file */
      strcpy( file, "hpccoutf.txt" );
      if     ( fid == 6 ) TEST->outfp = stdout;
      else if( fid == 7 ) TEST->outfp = stderr;
      else if( ( TEST->outfp = fopen( file, "a" ) ) == NULL )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "cannot open file %s.",
                    file );
         TEST->outfp = stderr;
         error = 1; goto label_error;
      }
      if (error == 1) goto label_error;
/*
 * Read and check the parameter values for the tests.
 *
 * Problem size (>=0) (N)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); 
      (void) sscanf( line, "%s", num ); *NS = atoi( num );
      if( ( *NS < 1 ) || ( *NS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %d",
                    "Number of values of N is less than 1 or greater than",
                    HPL_MAX_PARAM );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( N[ i ] = atoi( num ) ) < 0 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of N less than 0" );
            error = 1; goto label_error;
         }
      }
/*
 * Block size (>=1) (NB)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NBS = atoi( num );
      if( ( *NBS < 1 ) || ( *NBS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of NB is less than 1 or",
                    "greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NBS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( NB[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", 
                       "Value of NB less than 1" );
            error = 1; goto label_error;
         }
      }
/*
 * Process grids, mapping, (>=1) (P, Q)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num );
      *PMAPPIN = ( atoi( num ) == 1 ? HPL_COLUMN_MAJOR : HPL_ROW_MAJOR );

      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NPQS = atoi( num );
      if( ( *NPQS < 1 ) || ( *NPQS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of grids is less",
                    "than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NPQS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( P[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of P less than 1" );
            error = 1; goto label_error;
         }
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NPQS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( Q[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of Q less than 1" );
            error = 1; goto label_error;
         }
      }
/*
 * Check for enough processes in machine configuration
 */
      maxp = 0;
      for( i = 0; i < *NPQS; i++ )
      { nprocs   = P[i] * Q[i]; maxp = Mmax( maxp, nprocs ); }
      if( maxp > size )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                    "Need at least %d processes for these tests", maxp );
         error = 1; goto label_error;
      }
/*
 * Checking threshold value (TEST->thrsh)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); TEST->thrsh = atof( num );
/*
 * Panel factorization algorithm (PF)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NPFS = atoi( num );
      if( ( *NPFS < 1 ) || ( *NPFS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "number of values of PFACT",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NPFS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         j = atoi( num );
         if(      j == 0 ) PF[ i ] = HPL_LEFT_LOOKING;
         else if( j == 1 ) PF[ i ] = HPL_CROUT;
         else if( j == 2 ) PF[ i ] = HPL_RIGHT_LOOKING;
         else              PF[ i ] = HPL_RIGHT_LOOKING;
      }
/*
 * Recursive stopping criterium (>=1) (NBM)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NBMS = atoi( num );
      if( ( *NBMS < 1 ) || ( *NBMS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of NBMIN",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NBMS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( NBM[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of NBMIN less than 1" );
            error = 1; goto label_error;
         }
      }
/*
 * Number of panels in recursion (>=2) (NDV)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NDVS = atoi( num );
      if( ( *NDVS < 1 ) || ( *NDVS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of NDIV",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NDVS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( NDV[ i ] = atoi( num ) ) < 2 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of NDIV less than 2" );
            error = 1; goto label_error;
         }
      }
/*
 * Recursive panel factorization (RF)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NRFS = atoi( num );
      if( ( *NRFS < 1 ) || ( *NRFS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of RFACT",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NRFS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         j = atoi( num );
         if(      j == 0 ) RF[ i ] = HPL_LEFT_LOOKING;
         else if( j == 1 ) RF[ i ] = HPL_CROUT;
         else if( j == 2 ) RF[ i ] = HPL_RIGHT_LOOKING;
         else              RF[ i ] = HPL_RIGHT_LOOKING;
      }
/*
 * Broadcast topology (TP) (0=rg, 1=2rg, 2=rgM, 3=2rgM, 4=L)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NTPS = atoi( num );
      if( ( *NTPS < 1 ) || ( *NTPS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of BCAST",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NTPS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         j = atoi( num );
         if(      j == 0 ) TP[ i ] = HPL_1RING;
         else if( j == 1 ) TP[ i ] = HPL_1RING_M;
         else if( j == 2 ) TP[ i ] = HPL_2RING;
         else if( j == 3 ) TP[ i ] = HPL_2RING_M;
         else if( j == 4 ) TP[ i ] = HPL_BLONG;
         else if( j == 5 ) TP[ i ] = HPL_BLONG_M;
         else              TP[ i ] = HPL_1RING_M;
      }
/*
 * Lookahead depth (>=0) (NDH)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NDHS = atoi( num );
      if( ( *NDHS < 1 ) || ( *NDHS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of DEPTH",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NDHS; i++ )
      {
         (void) sscanf( lineptr, "%s", num );
         lineptr += strlen( num ) + 1;
         if( ( DH[ i ] = atoi( num ) ) < 0 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of DEPTH less than 0" );
            error = 1; goto label_error;
         }
      }
/*
 * Swapping algorithm (0,1 or 2) (FSWAP)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); j = atoi( num );
      if(      j == 0 ) *FSWAP = HPL_SWAP00;
      else if( j == 1 ) *FSWAP = HPL_SWAP01;
      else if( j == 2 ) *FSWAP = HPL_SW_MIX;
      else              *FSWAP = HPL_SWAP01;
/*
 * Swapping threshold (>=0) (TSWAP)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *TSWAP = atoi( num );
      if( *TSWAP <= 0 ) *TSWAP = 0;
/*
 * L1 in (no-)transposed form (0 or 1)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *L1NOTRAN = atoi( num );
      if( ( *L1NOTRAN != 0 ) && ( *L1NOTRAN != 1 ) ) *L1NOTRAN = 0; 
/*
 * U  in (no-)transposed form (0 or 1)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *UNOTRAN = atoi( num );
      if( ( *UNOTRAN != 0 ) && ( *UNOTRAN != 1 ) ) *UNOTRAN = 0;
/*
 * Equilibration (0=no, 1=yes)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *EQUIL = atoi( num );
      if( ( *EQUIL != 0 ) && ( *EQUIL != 1 ) ) *EQUIL = 1;
/*
 * Memory alignment in bytes (> 0) (ALIGN)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *ALIGN = atoi( num );
      if( *ALIGN <= 0 ) *ALIGN = 4;
/*
 * Close input file
 */
label_error:
      if (infp) fclose( infp );
   }
   else { TEST->outfp = NULL; }
/*
 * Check for error on reading input file
 */
   (void) HPL_all_reduce( (void *)(&error), 1, HPL_INT, HPL_max,
                          MPI_COMM_WORLD );
   if( error )
   {
     /*
      if( rank == 0 )
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                    "Illegal input in file " INFILE ". Exiting ..." );
      MPI_Finalize();
#ifdef HPL_CALL_VSIPL
      (void) vsip_finalize( NULL );
#endif
      exit( 1 );
      */
     HPCC_Defaults( TEST, /* use outfp, set threshold */
                    NS, N,
                    NBS, NB,
                    PMAPPIN,
                    NPQS, P, Q,
                    NPFS, PF,
                    NBMS, NBM,
                    NDVS, NDV,
                    NRFS, RF,
                    NTPS, TP,
                    NDHS, DH,
                    FSWAP,
                    TSWAP,
                    L1NOTRAN,
                    UNOTRAN,
                    EQUIL,
                    ALIGN,
                    MPI_COMM_WORLD );
   }
/*
 * Compute and broadcast machine epsilon
 */
   TEST->epsil = HPL_pdlamch( MPI_COMM_WORLD, HPL_MACH_EPS );
/*
 * Pack information arrays and broadcast
 */
   (void) HPL_broadcast( (void *)(&(TEST->thrsh)), 1, HPL_DOUBLE, 0,
                         MPI_COMM_WORLD );
/*
 * Broadcast array sizes
 */
   iwork = (int *)malloc( (size_t)(15) * sizeof( int ) );
   if( rank == 0 )
   {
      iwork[ 0] = *NS;      iwork[ 1] = *NBS;
      iwork[ 2] = ( *PMAPPIN == HPL_ROW_MAJOR ? 0 : 1 );
      iwork[ 3] = *NPQS;    iwork[ 4] = *NPFS;     iwork[ 5] = *NBMS;
      iwork[ 6] = *NDVS;    iwork[ 7] = *NRFS;     iwork[ 8] = *NTPS;
      iwork[ 9] = *NDHS;    iwork[10] = *TSWAP;    iwork[11] = *L1NOTRAN;
      iwork[12] = *UNOTRAN; iwork[13] = *EQUIL;    iwork[14] = *ALIGN;
   }
   (void) HPL_broadcast( (void *)iwork, 15, HPL_INT, 0, MPI_COMM_WORLD );
   if( rank != 0 )
   {
      *NS       = iwork[ 0]; *NBS   = iwork[ 1];
      *PMAPPIN  = ( iwork[ 2] == 0 ?  HPL_ROW_MAJOR : HPL_COLUMN_MAJOR );
      *NPQS     = iwork[ 3]; *NPFS  = iwork[ 4]; *NBMS     = iwork[ 5];
      *NDVS     = iwork[ 6]; *NRFS  = iwork[ 7]; *NTPS     = iwork[ 8];
      *NDHS     = iwork[ 9]; *TSWAP = iwork[10]; *L1NOTRAN = iwork[11];
      *UNOTRAN  = iwork[12]; *EQUIL = iwork[13]; *ALIGN    = iwork[14];
   }
   if( iwork ) free( iwork );
/*
 * Pack information arrays and broadcast
 */
   lwork = (*NS) + (*NBS) + 2 * (*NPQS) + (*NPFS) + (*NBMS) + 
           (*NDVS) + (*NRFS) + (*NTPS) + (*NDHS) + 1;
   iwork = (int *)malloc( (size_t)(lwork) * sizeof( int ) );
   if( rank == 0 )
   {
      j = 0;
      for( i = 0; i < *NS;   i++ ) { iwork[j] = N [i]; j++; }
      for( i = 0; i < *NBS;  i++ ) { iwork[j] = NB[i]; j++; }
      for( i = 0; i < *NPQS; i++ ) { iwork[j] = P [i]; j++; }
      for( i = 0; i < *NPQS; i++ ) { iwork[j] = Q [i]; j++; }
      for( i = 0; i < *NPFS; i++ )
      {
         if(      PF[i] == HPL_LEFT_LOOKING  ) iwork[j] = 0;
         else if( PF[i] == HPL_CROUT         ) iwork[j] = 1;
         else if( PF[i] == HPL_RIGHT_LOOKING ) iwork[j] = 2;
         j++;
      }
      for( i = 0; i < *NBMS; i++ ) { iwork[j] = NBM[i]; j++; }
      for( i = 0; i < *NDVS; i++ ) { iwork[j] = NDV[i]; j++; }
      for( i = 0; i < *NRFS; i++ )
      {
         if(      RF[i] == HPL_LEFT_LOOKING  ) iwork[j] = 0;
         else if( RF[i] == HPL_CROUT         ) iwork[j] = 1;
         else if( RF[i] == HPL_RIGHT_LOOKING ) iwork[j] = 2;
         j++;
      }
      for( i = 0; i < *NTPS; i++ )
      {
         if(      TP[i] == HPL_1RING   ) iwork[j] = 0;
         else if( TP[i] == HPL_1RING_M ) iwork[j] = 1;
         else if( TP[i] == HPL_2RING   ) iwork[j] = 2;
         else if( TP[i] == HPL_2RING_M ) iwork[j] = 3;
         else if( TP[i] == HPL_BLONG   ) iwork[j] = 4;
         else if( TP[i] == HPL_BLONG_M ) iwork[j] = 5;
         j++;
      }
      for( i = 0; i < *NDHS; i++ ) { iwork[j] = DH[i]; j++; }

      if(      *FSWAP == HPL_SWAP00 ) iwork[j] = 0;
      else if( *FSWAP == HPL_SWAP01 ) iwork[j] = 1;
      else if( *FSWAP == HPL_SW_MIX ) iwork[j] = 2;
      j++;
   }
   (void) HPL_broadcast( (void*)iwork, lwork, HPL_INT, 0,
                         MPI_COMM_WORLD );
   if( rank != 0 )
   {
      j = 0;
      for( i = 0; i < *NS;   i++ ) { N [i] = iwork[j]; j++; }
      for( i = 0; i < *NBS;  i++ ) { NB[i] = iwork[j]; j++; }
      for( i = 0; i < *NPQS; i++ ) { P [i] = iwork[j]; j++; }
      for( i = 0; i < *NPQS; i++ ) { Q [i] = iwork[j]; j++; }

      for( i = 0; i < *NPFS; i++ )
      {
         if(      iwork[j] == 0 ) PF[i] = HPL_LEFT_LOOKING;
         else if( iwork[j] == 1 ) PF[i] = HPL_CROUT;
         else if( iwork[j] == 2 ) PF[i] = HPL_RIGHT_LOOKING;
         j++;
      }
      for( i = 0; i < *NBMS; i++ ) { NBM[i] = iwork[j]; j++; }
      for( i = 0; i < *NDVS; i++ ) { NDV[i] = iwork[j]; j++; }
      for( i = 0; i < *NRFS; i++ )
      {
         if(      iwork[j] == 0 ) RF[i] = HPL_LEFT_LOOKING;
         else if( iwork[j] == 1 ) RF[i] = HPL_CROUT;
         else if( iwork[j] == 2 ) RF[i] = HPL_RIGHT_LOOKING;
         j++;
      }
      for( i = 0; i < *NTPS; i++ )
      {
         if(      iwork[j] == 0 ) TP[i] = HPL_1RING;
         else if( iwork[j] == 1 ) TP[i] = HPL_1RING_M;
         else if( iwork[j] == 2 ) TP[i] = HPL_2RING;
         else if( iwork[j] == 3 ) TP[i] = HPL_2RING_M;
         else if( iwork[j] == 4 ) TP[i] = HPL_BLONG;
         else if( iwork[j] == 5 ) TP[i] = HPL_BLONG_M;
         j++;
      }
      for( i = 0; i < *NDHS; i++ ) { DH[i] = iwork[j]; j++; }

      if(      iwork[j] == 0 ) *FSWAP = HPL_SWAP00;
      else if( iwork[j] == 1 ) *FSWAP = HPL_SWAP01;
      else if( iwork[j] == 2 ) *FSWAP = HPL_SW_MIX;
      j++;
   }
   if( iwork ) free( iwork );
/*
 * regurgitate input
 */
   if( rank == 0 )
   {
      HPL_fprintf( TEST->outfp, "%s%s\n",
                   "========================================",
                   "========================================" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "HPLinpack 2.0  --  High-Performance Linpack benchmark  --  ",
          " September 10, 2008" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Written by A. Petitet and R. Clint Whaley,  ",
          "Innovative Computing Laboratory, UTK" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Modified by Piotr Luszczek, ",
          "Innovative Computing Laboratory, UTK" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Modified by Julien Langou, ",
          "University of Colorado Denver");
      HPL_fprintf( TEST->outfp, "%s%s\n",
                   "========================================",
                   "========================================" );

      HPL_fprintf( TEST->outfp, "\n%s\n",
          "An explanation of the input/output parameters follows:" );
      HPL_fprintf( TEST->outfp, "%s\n",
          "T/V    : Wall time / encoded variant." );
      HPL_fprintf( TEST->outfp, "%s\n",
         "N      : The order of the coefficient matrix A." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "NB     : The partitioning blocking factor." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "P      : The number of process rows." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "Q      : The number of process columns." );
      HPL_fprintf( TEST->outfp, "%s\n",
         "Time   : Time in seconds to solve the linear system." );
      HPL_fprintf( TEST->outfp, "%s\n\n",
         "Gflops : Rate of execution for solving the linear system." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "The following parameter values will be used:" );
/*
 * Problem size
 */
      HPL_fprintf( TEST->outfp,       "\nN      :" );
      for( i = 0; i < Mmin( 8, *NS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", N[i]  );
      if( *NS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", N[i]  );
         if( *NS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", N[i]  );
         }
      }
/*
 * Distribution blocking factor
 */
      HPL_fprintf( TEST->outfp,       "\nNB     :" );
      for( i = 0; i < Mmin( 8, *NBS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", NB[i] );
      if( *NBS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NBS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", NB[i] );
         if( *NBS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NBS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", NB[i] );
         }
      }
/*
 * Process mapping
 */
      HPL_fprintf( TEST->outfp,       "\nPMAP   :" );
      if(      *PMAPPIN == HPL_ROW_MAJOR    )
         HPL_fprintf( TEST->outfp, " Row-major process mapping" );
      else if( *PMAPPIN == HPL_COLUMN_MAJOR )
         HPL_fprintf( TEST->outfp, " Column-major process mapping" );
/*
 * Process grid
 */
      HPL_fprintf( TEST->outfp,       "\nP      :" );
      for( i = 0; i < Mmin( 8, *NPQS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", P[i]  );
      if( *NPQS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NPQS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", P[i]  );
         if( *NPQS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NPQS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", P[i]  );
         }
      }
      HPL_fprintf( TEST->outfp,       "\nQ      :" );
      for( i = 0; i < Mmin( 8, *NPQS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", Q[i]  );
      if( *NPQS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NPQS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", Q[i]  );
         if( *NPQS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NPQS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", Q[i]  );
         }
      }
/*
 * Panel Factorization
 */
      HPL_fprintf( TEST->outfp,       "\nPFACT  :" );
      for( i = 0; i < Mmin( 8, *NPFS ); i++ )
      {
         if(      PF[i] == HPL_LEFT_LOOKING  )
            HPL_fprintf( TEST->outfp,       "    Left " );
         else if( PF[i] == HPL_CROUT         )
            HPL_fprintf( TEST->outfp,       "   Crout " );
         else if( PF[i] == HPL_RIGHT_LOOKING )
            HPL_fprintf( TEST->outfp,       "   Right " );
      }
      if( *NPFS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NPFS ); i++ )
         {
            if(      PF[i] == HPL_LEFT_LOOKING  )
               HPL_fprintf( TEST->outfp,       "    Left " );
            else if( PF[i] == HPL_CROUT         )
               HPL_fprintf( TEST->outfp,       "   Crout " );
            else if( PF[i] == HPL_RIGHT_LOOKING )
               HPL_fprintf( TEST->outfp,       "   Right " );
         }
         if( *NPFS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NPFS; i++ )
            {
               if(      PF[i] == HPL_LEFT_LOOKING  )
                  HPL_fprintf( TEST->outfp,       "    Left " );
               else if( PF[i] == HPL_CROUT         )
                  HPL_fprintf( TEST->outfp,       "   Crout " );
               else if( PF[i] == HPL_RIGHT_LOOKING )
                  HPL_fprintf( TEST->outfp,       "   Right " );
            }
         }
      }
/*
 * Recursive stopping criterium
 */
      HPL_fprintf( TEST->outfp,       "\nNBMIN  :" );
      for( i = 0; i < Mmin( 8, *NBMS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", NBM[i]  );
      if( *NBMS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NBMS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", NBM[i]  );
         if( *NBMS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NBMS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", NBM[i]  );
         }
      }
/*
 * Number of panels in recursion
 */
      HPL_fprintf( TEST->outfp,       "\nNDIV   :" );
      for( i = 0; i < Mmin( 8, *NDVS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", NDV[i]  );
      if( *NDVS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NDVS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", NDV[i]  );
         if( *NDVS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NDVS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", NDV[i]  );
         }
      }
/*
 * Recursive Factorization
 */
      HPL_fprintf( TEST->outfp,       "\nRFACT  :" );
      for( i = 0; i < Mmin( 8, *NRFS ); i++ )
      {
         if(      RF[i] == HPL_LEFT_LOOKING  )
            HPL_fprintf( TEST->outfp,       "    Left " );
         else if( RF[i] == HPL_CROUT         )
            HPL_fprintf( TEST->outfp,       "   Crout " );
         else if( RF[i] == HPL_RIGHT_LOOKING )
            HPL_fprintf( TEST->outfp,       "   Right " );
      }
      if( *NRFS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NRFS ); i++ )
         {
            if(      RF[i] == HPL_LEFT_LOOKING  )
               HPL_fprintf( TEST->outfp,       "    Left " );
            else if( RF[i] == HPL_CROUT         )
               HPL_fprintf( TEST->outfp,       "   Crout " );
            else if( RF[i] == HPL_RIGHT_LOOKING )
               HPL_fprintf( TEST->outfp,       "   Right " );
         }
         if( *NRFS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NRFS; i++ )
            {
               if(      RF[i] == HPL_LEFT_LOOKING  )
                  HPL_fprintf( TEST->outfp,       "    Left " );
               else if( RF[i] == HPL_CROUT         )
                  HPL_fprintf( TEST->outfp,       "   Crout " );
               else if( RF[i] == HPL_RIGHT_LOOKING )
                  HPL_fprintf( TEST->outfp,       "   Right " );
            }
         }
      }
/*
 * Broadcast topology
 */
      HPL_fprintf( TEST->outfp,       "\nBCAST  :" );
      for( i = 0; i < Mmin( 8, *NTPS ); i++ )
      {
         if(      TP[i] == HPL_1RING   )
            HPL_fprintf( TEST->outfp,       "   1ring " );
         else if( TP[i] == HPL_1RING_M )
            HPL_fprintf( TEST->outfp,       "  1ringM " );
         else if( TP[i] == HPL_2RING   )
            HPL_fprintf( TEST->outfp,       "   2ring " );
         else if( TP[i] == HPL_2RING_M )
            HPL_fprintf( TEST->outfp,       "  2ringM " );
         else if( TP[i] == HPL_BLONG   )
            HPL_fprintf( TEST->outfp,       "   Blong " );
         else if( TP[i] == HPL_BLONG_M )
            HPL_fprintf( TEST->outfp,       "  BlongM " );
      }
      if( *NTPS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NTPS ); i++ )
         {
            if(      TP[i] == HPL_1RING   )
               HPL_fprintf( TEST->outfp,       "   1ring " );
            else if( TP[i] == HPL_1RING_M )
               HPL_fprintf( TEST->outfp,       "  1ringM " );
            else if( TP[i] == HPL_2RING   )
               HPL_fprintf( TEST->outfp,       "   2ring " );
            else if( TP[i] == HPL_2RING_M )
               HPL_fprintf( TEST->outfp,       "  2ringM " );
            else if( TP[i] == HPL_BLONG   )
               HPL_fprintf( TEST->outfp,       "   Blong " );
            else if( TP[i] == HPL_BLONG_M )
               HPL_fprintf( TEST->outfp,       "  BlongM " );
         }
         if( *NTPS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NTPS; i++ )
            {
               if(      TP[i] == HPL_1RING   )
                  HPL_fprintf( TEST->outfp,       "   1ring " );
               else if( TP[i] == HPL_1RING_M )
                  HPL_fprintf( TEST->outfp,       "  1ringM " );
               else if( TP[i] == HPL_2RING   )
                  HPL_fprintf( TEST->outfp,       "   2ring " );
               else if( TP[i] == HPL_2RING_M )
                  HPL_fprintf( TEST->outfp,       "  2ringM " );
               else if( TP[i] == HPL_BLONG   )
                  HPL_fprintf( TEST->outfp,       "   Blong " );
               else if( TP[i] == HPL_BLONG_M )
                  HPL_fprintf( TEST->outfp,       "  BlongM " );
            }
         }
      }
/*
 * Lookahead depths
 */
      HPL_fprintf( TEST->outfp,       "\nDEPTH  :" );
      for( i = 0; i < Mmin( 8, *NDHS ); i++ )
         HPL_fprintf( TEST->outfp,       "%8d ", DH[i]  );
      if( *NDHS > 8 )
      {
         HPL_fprintf( TEST->outfp,    "\n        " );
         for( i = 8; i < Mmin( 16, *NDHS ); i++ )
            HPL_fprintf( TEST->outfp,    "%8d ", DH[i]  );
         if( *NDHS > 16 )
         {
            HPL_fprintf( TEST->outfp, "\n        " );
            for( i = 16; i < *NDHS; i++ )
               HPL_fprintf( TEST->outfp, "%8d ", DH[i]  );
         }
      }
/*
 * Swapping algorithm
 */
      HPL_fprintf( TEST->outfp,       "\nSWAP   :" );
      if(      *FSWAP == HPL_SWAP00 )
         HPL_fprintf( TEST->outfp, " Binary-exchange" );
      else if( *FSWAP == HPL_SWAP01 )
         HPL_fprintf( TEST->outfp, " Spread-roll (long)" );
      else if( *FSWAP == HPL_SW_MIX )
         HPL_fprintf( TEST->outfp, " Mix (threshold = %d)", *TSWAP );
/*
 * L1 storage form
 */
      HPL_fprintf( TEST->outfp,       "\nL1     :" );
      if(      *L1NOTRAN != 0 )
         HPL_fprintf( TEST->outfp, " no-transposed form" );
      else
         HPL_fprintf( TEST->outfp, " transposed form" );
/*
 * U  storage form
 */
      HPL_fprintf( TEST->outfp,       "\nU      :" );
      if(      *UNOTRAN != 0 )
         HPL_fprintf( TEST->outfp, " no-transposed form" );
      else
         HPL_fprintf( TEST->outfp, " transposed form" );
/*
 * Equilibration
 */
      HPL_fprintf( TEST->outfp,       "\nEQUIL  :" );
      if(      *EQUIL != 0 )
         HPL_fprintf( TEST->outfp, " yes" );
      else
         HPL_fprintf( TEST->outfp, " no" );
/*
 * Alignment
 */
      HPL_fprintf( TEST->outfp,       "\nALIGN  : %d double precision words",
                   *ALIGN );

      HPL_fprintf( TEST->outfp, "\n\n" );
/*
 * For testing only
 */
      if( TEST->thrsh > HPL_rzero )
      {
         HPL_fprintf( TEST->outfp, "%s%s\n\n",
                      "----------------------------------------",
                      "----------------------------------------" );
         HPL_fprintf( TEST->outfp, "%s\n",
            "- The matrix A is randomly generated for each test." );
         HPL_fprintf( TEST->outfp, "%s\n",
            "- The following scaled residual check will be computed:" );
         HPL_fprintf( TEST->outfp, "%s\n",
            "      ||Ax-b||_oo / ( eps * ( || x ||_oo * || A ||_oo + || b ||_oo ) * N )" );
         HPL_fprintf( TEST->outfp, "%s %21.6e\n",
            "- The relative machine precision (eps) is taken to be     ",
            TEST->epsil );
         HPL_fprintf( TEST->outfp, "%s   %11.1f\n\n",
            "- Computational tests pass if scaled residuals are less than      ",
            TEST->thrsh );
      }
   }
/*
 * End of HPL_pdinfo
 */
}


/* Dependant functions */
/** tstdgemm.c **/

/* Generates random matrix with entries between 0.0 and 1.0 */
static
void dmatgen(int m, int n, double *a, int lda, int seed) {
  int i, j;
  double *a0 = a, rcp = 1.0 / RAND_MAX;

  srand( seed );

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      a0[i] = rcp * rand();

    a0 += lda;
  }
}

static
double dnrm_inf(int m, int n, double *a, int lda) {
  int i, j, k, lnx;
  double mx, *a0;

  int nx = 10;
  double x[10];

  mx = 0.0;

  for (i = 0; i < m; i += nx) {
    lnx = Mmin( nx, m-i );
    for (k = 0; k < lnx; ++k) x[k] = 0.0;

    a0 = a + i;

    for (j = 0; j < n; ++j) {
      for (k = 0; k < lnx; ++k)
        x[k] += fabs( a0[k] );

      a0 += lda;
    }

    for (k = 0; k < lnx; ++k)
      if (mx < x[k]) mx = x[k];
  }

  return mx;
}
/** mem.c **/

static int
CheckNode(int imrow, int imcol, int nmat, int *mval, int *nval, int nbmat, int *mbval, int *nbval,
          int myrow, int mycol, int nprow, int npcol, long *maxMem) {
  int i__, ii, m, n, mb, nb, ierr[1];
  int lcm, np0, nq0, mp0, mq0, mg, ng, np, nq, mp, mq;
  long isw, ipw, ipiw, ipa, ipc;

  *maxMem = 0;
  for (i__ = 0; i__ < nmat; ++i__) {
    m = mval[i__];
    n = nval[i__];

/*           Make sure matrix information is correct */

      ierr[0] = 0;
      if (m < 1) {
        ierr[0] = 1;
      } else if (n < 1) {
        ierr[0] = 1;
      }

      if (ierr[0] > 0) {
        continue;
      }

      for (ii = 0; ii < nbmat; ++ii) { /* Loop over different block sizes */

        mb = mbval[ii];
        nb = nbval[ii];

/*              Make sure blocking sizes are legal */
        ierr[0] = 0;
        if (mb < 1) {
          ierr[0] = 1;
        } else if (nb < 1) {
          ierr[0] = 1;
        }

/*              Make sure no one had error */

        if (ierr[0] > 0) {
          continue;
        }

        mp = numroc_(&m, &mb, &myrow, &imrow, &nprow);
        mq = numroc_(&m, &mb, &mycol, &imcol, &npcol);
        np = numroc_(&n, &nb, &myrow, &imrow, &nprow);
        nq = numroc_(&n, &nb, &mycol, &imcol, &npcol);

        mg = iceil_(&m, &mb);
        ng = iceil_(&n, &nb);

        mp0 = iceil_(&mg, &nprow) * mb;
        mq0 = iceil_(&mg, &npcol) * mb;
        np0 = iceil_(&ng, &nprow) * nb;
        nq0 = iceil_(&ng, &npcol) * nb;

        lcm = ilcm_(&nprow, &npcol);
        ipc = 1;
        ipa = ipc + (long)np0 * (long)mq0;
        ipiw = (long)mp0 * (long)nq0 + ipa;
        ipw = ipiw;
        isw = ipw + (long)(iceil_(&mg, &lcm) << 1) * (long)mb * (long)iceil_(&ng, &lcm) * (long)nb;

        if (*maxMem < isw) *maxMem = isw;
      }
  }
  return 0;
}

int
MaxMem(int nprocs, int imrow, int imcol, int nmat, int *mval, int *nval, int nbmat, int *mbval,
       int *nbval, int ngrids, int *npval, int *nqval, long *maxMem) {
  int nprow, npcol, myrow, mycol;
  int j, ierr[1];
  long curMem;

  *maxMem = 0;
  for (j = 0; j < ngrids; ++j) {
    nprow = npval[j];
    npcol = nqval[j];

/*        Make sure grid information is correct */

    ierr[0] = 0;
    if (nprow < 1) {
     ierr[0] = 1;
    } else if (npcol < 1) {
      ierr[0] = 1;
    } else if (nprow * npcol > nprocs) {
      ierr[0] = 1;
    }

    if (ierr[0] > 0) {
      continue;
    }
    for (myrow = 0; myrow < nprow; myrow++)
      for (mycol = 0; mycol < npcol; mycol++) {
        CheckNode( imrow, imcol, nmat, mval, nval, nbmat, mbval, nbval, myrow, mycol, nprow,
                   npcol, &curMem );
        if (*maxMem < curMem) *maxMem = curMem;
      }
  }
  return 0;
}

// TODO Possible main interference
#ifdef HPCC_MEMMAIN
int iceil_(int *n,int *d) {return *n>0 ? (*n+*d-1)/ *d : *n/ *d;}
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs) {
  int ret_val, extrablks, mydist, nblocks;
  mydist = (*nprocs + *iproc - *isrcproc) % *nprocs;
  nblocks = *n / *nb;
  ret_val = nblocks / *nprocs * *nb;
  extrablks = nblocks % *nprocs;
  if (mydist < extrablks) {
    ret_val += *nb;
  } else if (mydist == extrablks) {
    ret_val += *n % *nb;
  }
  return ret_val;
}
int ilcm_(int *m, int *n) {
  int ret_val;
  int ia, iq, ir;
  if (*m >= *n) {
    ia = *m;
    ret_val = *n;
  } else {
    ia = *n;
    ret_val = *m;
  }
  for (;;) {
    iq = ia / ret_val;
    ir = ia - iq * ret_val;
    if (ir == 0) {
      ret_val = *m * *n / ret_val;
      return ret_val;
    }
    ia = ret_val;
    ret_val = ir;
  }
}
/* Commented because of CBLAS
int
main(int argc, char *argv[]) {
  int n, nb, nprow, npcol, ng, lcm;
  int nval[1], nbval[1];
  long maxMem;

  if (argc <= 1) {
    printf( "Usage:\n%s n nb nprow npcol\n", argv[0] );
  }

  if (argc <= 1 || sscanf( argv[1], "%d", &n  ) != 1 || n < 1)  n = 50000;
  if (argc <= 2 || sscanf( argv[2], "%d", &nb ) != 1 || nb < 1) nb = 80;
  if (argc <= 3 || sscanf( argv[3], "%d", &nprow ) != 1 || nprow < 1) nprow = 8;
  if (argc <= 4 || sscanf( argv[4], "%d", &npcol ) != 1 || npcol < 1) npcol = nprow;

  nval[0] = n;
  nbval[0] = nb;

  CheckNode( 0, 0, 1, nval, nval, 1, nbval, nbval, 0, 0, nprow, npcol, &maxMem );

  printf( "n=%d nb=%d nprow=%d npcol=%d lcm(nprow,npcol)=%d\n%ld\n", n, nb, nprow, npcol,
          ilcm_(&nprow, &npcol), maxMem );

  ng = iceil_(&n, &nb);
  lcm = ilcm_(&nprow, &npcol);
  printf( "%d %d %d\n", ng, lcm, (iceil_(&ng, &lcm) << 1) * nb * iceil_(&ng, &lcm) * nb );
  printf( "%d %d\n", (iceil_(&ng, &lcm) << 1), iceil_(&ng, &lcm) );

  return 0;
} */
#endif

/** noopt.c **/

double
HPCC_dweps() {
  double eps, one, half;

  one = 1.0;
  half = one / 2.0;
  eps = one;

  while (one + eps != one)
    eps *= half;

  return eps;
}

float
HPCC_sweps() {
  float eps, one, half;

  one = 1.0f;
  half = one / 2.0f;
  eps = one;

  while (one + eps != one)
    eps *= half;

  return eps;
}

/** io.c **/
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double HPCC_MemProc = -1.0, HPCC_MemVal = -1.0;
static int HPCC_MemSpec = -1;

static int
ReadInts(char *buf, int n, int *val) {
  int i, j;

  for (j = i = 0; i < n; i++) {
    if (sscanf( buf + j, "%d", val + i ) != 1) {
      i--;
      break;
    }
    for (; buf[j] && isdigit(buf[j]); j++)
      ; /* EMPTY */
    for (; buf[j] && ! isdigit(buf[j]); j++)
      ; /* EMPTY */
    if (! buf[j]) {
      i--;
      break;
    }
  }

  return i + 1;
}

static int
HPCC_InitHPL(HPCC_Params *p) {
  HPL_pdinfo( &p->test, &p->ns, p->nval, &p->nbs, p->nbval, &p->porder, &p->npqs, p->pval,
              p->qval, &p->npfs, p->pfaval, &p->nbms, p->nbmval, &p->ndvs, p->ndvval, &p->nrfs,
              p->rfaval, &p->ntps, p->topval, &p->ndhs, p->ndhval, &p->fswap, &p->tswap,
              &p->L1notran, &p->Unotran, &p->equil, &p->align );

  if (p->test.thrsh <= 0.0) p->Failure = 1;

  return 0;
}

static int
iiamax(int n, int *x, int incx) {
  int i, v, mx, idx = 0;

  idx = 0;
  mx = (x[0] < 0 ? -x[0] : x[0]);
  for (i = 0; i < n; i += incx) {
    v = (x[i] < 0 ? -x[i] : x[i]);
    if (mx < v) {mx = v; idx = i;}
  }

  return idx;
}

static void
icopy(int n, int *src, int sinc, int *dst, int dinc) {
  int i;

  for (i = n; i; i--) {
    *dst = *src;
    dst += dinc;
    src += sinc;
  }
}

int
HPCC_InputFileInit(HPCC_Params *params) {
  int myRank, commSize;
  int i, j, n, ioErr, lastConfigLine = 32, line, rv, maxHPLn;
  char buf[82]; int nbuf = 82;
  FILE *f, *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  if (0 == myRank) {
    f = fopen( params->inFname, "r" );
    if (! f) {
      ioErr = 1;
      goto ioEnd;
    }

    /* skip irrelevant lines in config file */
            for (line = 0; line < lastConfigLine; line++)
              if (! fgets( buf, nbuf, f )) break;

            if (line < lastConfigLine) { /* if didn't read all the required lines */
      ioErr = 1;
      goto ioEnd;
    }

    /* Get values of N for PTRANS */
    line++;
    fgets( buf, nbuf, f );
    rv = sscanf( buf, "%d", &n );
    if (rv != 1 || n < 0) { /* parse error or negative value*/
      n = 0;
      BEGIN_IO(myRank, params->outFname, outputFile);
      fprintf( outputFile, "Error in line %d of the input file.\n", line );
      END_IO( myRank, outputFile );
    }
    n = Mmin( n, HPL_MAX_PARAM );

    line++;
    fgets( buf, nbuf, f );
    ReadInts( buf, n, params->PTRANSnval );

    /* find the largest matrix for HPL */
    maxHPLn = params->nval[iiamax( params->ns, params->nval, 1 )];

    for (j = i = 0; i < n; i++) {
      /* if memory for PTRANS is at least 90% of what's used for HPL */
      if (params->PTRANSnval[i] >= 0.9486 * maxHPLn * 0.5) {
        params->PTRANSnval[j] = params->PTRANSnval[i];
        j++;
      }
    }
    n = j; /* only this many entries use enough memory */

    /* copy matrix sizes from HPL, divide by 2 so both PTRANS matrices (plus "work" arrays) occupy
       as much as HPL's one */
    for (i = 0; i < params->ns; i++)
      params->PTRANSnval[i + n] = params->nval[i] / 2;
    params->PTRANSns = n + params->ns;

    /* Get values of block sizes */
    line++;
    fgets( buf, nbuf, f );
    rv = sscanf( buf, "%d", &n );
    if (rv != 1 || n < 0) { /* parse error or negative value*/
      n = 0;
      BEGIN_IO( myRank, params->outFname, outputFile );
      fprintf( outputFile, "Error in line %d of the input file.\n", line );
      END_IO( myRank, outputFile );
    }
    n = Mmin( n, HPL_MAX_PARAM );

    line++;
    fgets( buf, nbuf, f );
    ReadInts( buf, n, params->PTRANSnbval );

    icopy( params->nbs, params->nbval, 1, params->PTRANSnbval + n, 1 );
    params->PTRANSnbs = n + params->nbs;

    ioErr = 0;
    ioEnd:
    if (f) fclose( f );
  }

  MPI_Bcast( &ioErr, 1, MPI_INT, 0, comm );
  if (ioErr) {
    /* copy matrix sizes from HPL, divide by 2 so both PTRANS matrices (plus "work" arrays) occupy
       as much as HPL's one */
    for (i = 0; i < params->ns; i++)
      params->PTRANSnval[i] = params->nval[i] / 2;
    params->PTRANSns = params->ns;

    icopy( params->nbs, params->nbval, 1, params->PTRANSnbval, 1 );
    params->PTRANSnbs = params->nbs;
  }

  /* broadcast what's been read on node 0 */
  MPI_Bcast( &params->PTRANSns, 1, MPI_INT, 0, comm );
  if (params->PTRANSns > 0)
    MPI_Bcast( &params->PTRANSnval, params->PTRANSns, MPI_INT, 0, comm );
  MPI_Bcast( &params->PTRANSnbs, 1, MPI_INT, 0, comm );
  if (params->PTRANSnbs > 0)
    MPI_Bcast( &params->PTRANSnbval, params->PTRANSnbs, MPI_INT, 0, comm );

  /* copy what HPL has */
  params->PTRANSnpqs = params->npqs;
  icopy( params->npqs, params->qval, 1, params->PTRANSqval, 1 );
  icopy( params->npqs, params->pval, 1, params->PTRANSpval, 1 );

  return ioErr;
}

static int
ErrorReduce(FILE *f, char *str, int eCode, MPI_Comm comm) {
  int rCode;

  if (eCode) eCode = 1; /* make sure error is indicated with 1 */

  MPI_Allreduce( &eCode, &rCode, 1, MPI_INT, MPI_SUM, comm );

  if (rCode) {
    if (f)
      fprintf( f, "%s", str );

    return -1;
  }

  return 0;
}

int
HPCC_Init(HPCC_Params *params) {
  int myRank, commSize;
  int i, nMax, nbMax, procCur, procMax, procMin, errCode;
  double totalMem;
  char inFname[12] = "hpccinf.txt", outFname[13] = "hpccoutf.txt";
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;
  time_t currentTime;
  char hostname[MPI_MAX_PROCESSOR_NAME + 1]; int hostnameLen;
#ifdef HPCC_MEMALLCTR
  size_t hpl_mem, ptrans_mem;
  long dMemSize;
#endif

  outputFile = NULL;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  strcpy( params->inFname, inFname );
  strcpy( params->outFname, outFname );

  if (0 == myRank)
    outputFile = fopen( params->outFname, "a" );

  errCode = 0;
  if (sizeof(u64Int) < 8 || sizeof(s64Int) < 8) errCode = 1;
  if (ErrorReduce( outputFile, "No 64-bit integer type available.", errCode, comm ))
    return -1;

  i = MPI_Get_processor_name( hostname, &hostnameLen );
  if (i) hostname[0] = 0;
  else hostname[Mmax(hostnameLen, MPI_MAX_PROCESSOR_NAME)] = 0;
  time( &currentTime );

  BEGIN_IO( myRank, params->outFname, outputFile );
  fprintf( outputFile,
            "########################################################################\n" );
  fprintf( outputFile,
            "This is the DARPA/DOE HPC Challenge Benchmark version %d.%d.%d October 2012\n",
            HPCC_VERSION_MAJOR, HPCC_VERSION_MINOR, HPCC_VERSION_MICRO );
  fprintf( outputFile, "Produced by Jack Dongarra and Piotr Luszczek\n" );
  fprintf( outputFile, "Innovative Computing Laboratory\n" );
  fprintf( outputFile, "University of Tennessee Knoxville and Oak Ridge National Laboratory\n\n" );
  fprintf( outputFile, "See the source files for authors of specific codes.\n" );
  fprintf( outputFile, "Compiled on %s at %s\n", __DATE__ , __TIME__ );
  fprintf( outputFile, "Current time (%ld) is %s\n",(long)currentTime,ctime(&currentTime));
  fprintf( outputFile, "Hostname: '%s'\n", hostname );
  fprintf( outputFile,
            "########################################################################\n" );
  END_IO( myRank, outputFile );

  params->Failure = 0;

  HPCC_InitHPL( params ); /* HPL calls exit() if there is a problem */
  HPCC_InputFileInit( params );

  params->RunHPL = 0;
  params->RunStarDGEMM = 0;
  params->RunSingleDGEMM = 0;
  params->RunPTRANS = 0;
  params->RunStarStream = 0;
  params->RunSingleStream = 0;
  params->RunMPIRandomAccess_LCG = 0;
  params->RunStarRandomAccess_LCG = 0;
  params->RunSingleRandomAccess_LCG = 0;
  params->RunMPIRandomAccess = 0;
  params->RunStarRandomAccess = 0;
  params->RunSingleRandomAccess = 0;
  params->RunLatencyBandwidth = 0;
  params->RunMPIFFT = 0;
  params->RunStarFFT = 0;
  params->RunSingleFFT = 0;
  params->RunHPL = params->RunStarDGEMM = params->RunSingleDGEMM =
  params->RunPTRANS = params->RunStarStream = params->RunSingleStream =
  params->RunMPIRandomAccess_LCG = params->RunStarRandomAccess_LCG = params->RunSingleRandomAccess_LCG =
  params->RunMPIRandomAccess = params->RunStarRandomAccess = params->RunSingleRandomAccess =
  params->RunMPIFFT = params->RunStarFFT = params->RunSingleFFT =
  params->RunLatencyBandwidth = 1;

  params->MPIRandomAccess_LCG_GUPs =
  params->MPIRandomAccess_GUPs = params->StarGUPs = params->SingleGUPs =
  params->StarDGEMMGflops = params->SingleDGEMMGflops = -1.0;
  params->StarStreamCopyGBs = params->StarStreamScaleGBs = params->StarStreamAddGBs =
  params->StarStreamTriadGBs = params->SingleStreamCopyGBs = params->SingleStreamScaleGBs =
  params->SingleStreamAddGBs = params->SingleStreamTriadGBs =
  params->SingleFFTGflops = params->StarFFTGflops = params->MPIFFTGflops = params->MPIFFT_maxErr =
  params->MaxPingPongLatency = params-> RandomlyOrderedRingLatency = params-> MinPingPongBandwidth =
  params->NaturallyOrderedRingBandwidth = params->RandomlyOrderedRingBandwidth =
  params->MinPingPongLatency = params->AvgPingPongLatency = params->MaxPingPongBandwidth =
  params->AvgPingPongBandwidth = params->NaturallyOrderedRingLatency = -1.0;

  params->HPLrdata.Gflops = -1000.0;
  params->HPLrdata.time = params->HPLrdata.eps = params->HPLrdata.RnormI = params->HPLrdata.Anorm1 = params->HPLrdata.AnormI = params->HPLrdata.Xnorm1 = params->HPLrdata.XnormI = -1.0;
  params->HPLrdata.N = params->HPLrdata.NB = params->HPLrdata.nprow = params->HPLrdata.npcol = params->HPLrdata.depth = params->HPLrdata.nbdiv = params->HPLrdata.nbmin = -1;
  params->HPLrdata.cpfact = params->HPLrdata.crfact = params->HPLrdata.ctop = params->HPLrdata.order = '-';

  params->PTRANSrdata.GBs = params->PTRANSrdata.time = params->PTRANSrdata.residual = -1.0;
  params->PTRANSrdata.n = params->PTRANSrdata.nb = params->PTRANSrdata.nprow =
  params->PTRANSrdata.npcol = -1;

  params->MPIRandomAccess_LCG_ErrorsFraction =
  params->MPIRandomAccess_ErrorsFraction =
  params->MPIRandomAccess_LCG_time = params->MPIRandomAccess_LCG_CheckTime =
  params->MPIRandomAccess_time = params->MPIRandomAccess_CheckTime =
  params->MPIRandomAccess_LCG_TimeBound =
  params->MPIRandomAccess_TimeBound = -1.0;

  params->DGEMM_N =
  params->FFT_N =
  params->StreamVectorSize =
  params->MPIRandomAccess_LCG_Algorithm =
  params->MPIRandomAccess_Algorithm =
  params->MPIFFT_Procs = -1;

  params->StreamThreads = 1;

  params->FFTEnblk = params->FFTEnp = params->FFTEl2size = -1;

  params->MPIFFT_N =
  params->RandomAccess_LCG_N =
  params->MPIRandomAccess_LCG_N =
  params->MPIRandomAccess_LCG_Errors =
  params->RandomAccess_N =
  params->MPIRandomAccess_N =
  params->MPIRandomAccess_Errors =
  params->MPIRandomAccess_LCG_ExeUpdates =
  params->MPIRandomAccess_ExeUpdates = (s64Int)(-1);

  procMax = procMin = params->pval[0] * params->qval[0];
  for (i = 1; i < params->npqs; ++i) {
    procCur = params->pval[i] * params->qval[i];
    if (procMax < procCur) procMax = procCur;
    if (procMin > procCur) procMin = procCur;
  }
  params->HPLMaxProc = procMax;
  params->HPLMinProc = procMin;

  nMax = params->nval[iiamax( params->ns, params->nval, 1 )];

  /* totalMem = (nMax*nMax) * sizeof(double) */
  totalMem = nMax;
  totalMem *= nMax;
  totalMem *= sizeof(double);
  params->HPLMaxProcMem = totalMem / procMin;

  for (i = 0; i < MPIFFT_TIMING_COUNT; i++)
    params->MPIFFTtimingsForward[i] = 0.0;

  i = iiamax( params->PTRANSnbs, params->PTRANSnbval, 1 );
  nbMax = params->PTRANSnbval[i];

#ifdef HPCC_MEMALLCTR
  MaxMem( commSize, 0, 0, params->PTRANSns, params->PTRANSnval, params->PTRANSnval, params->PTRANSnbs, params->PTRANSnbval, params->PTRANSnbval, params->PTRANSnpqs, params->PTRANSpval, params->PTRANSqval, &dMemSize );
  ptrans_mem = dMemSize * sizeof(double) + 3 * commSize * sizeof(int);
  hpl_mem = params->HPLMaxProcMem + (nMax + nbMax) * sizeof(double) * nbMax;
  HPCC_alloc_init( Mmax( ptrans_mem, hpl_mem ) );
#endif

  return 0;
}

int
HPCC_Finalize(HPCC_Params *params) {
  int myRank, commSize;
  int i;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;
  time_t currentTime;

#ifdef HPCC_MEMALLCTR
  HPCC_alloc_finalize();
#endif

  time( &currentTime );

  MPI_Comm_rank( comm, &myRank );
  MPI_Comm_size( comm, &commSize );

  BEGIN_IO(myRank, params->outFname, outputFile);
  fprintf( outputFile, "Begin of Summary section.\n" );
  fprintf( outputFile, "VersionMajor=%d\n", HPCC_VERSION_MAJOR );
  fprintf( outputFile, "VersionMinor=%d\n", HPCC_VERSION_MINOR );
  fprintf( outputFile, "VersionMicro=%d\n", HPCC_VERSION_MICRO );
  fprintf( outputFile, "VersionRelease=%c\n", HPCC_VERSION_RELEASE );
  fprintf( outputFile, "LANG=%s\n", "C" );
  fprintf( outputFile, "Success=%d\n", params->Failure ? 0 : 1 );
  fprintf( outputFile, "sizeof_char=%d\n", (int)sizeof(char) );
  fprintf( outputFile, "sizeof_short=%d\n", (int)sizeof(short) );
  fprintf( outputFile, "sizeof_int=%d\n", (int)sizeof(int) );
  fprintf( outputFile, "sizeof_long=%d\n", (int)sizeof(long) );
  fprintf( outputFile, "sizeof_void_ptr=%d\n", (int)sizeof(void*) );
  fprintf( outputFile, "sizeof_size_t=%d\n", (int)sizeof(size_t) );
  fprintf( outputFile, "sizeof_float=%d\n", (int)sizeof(float) );
  fprintf( outputFile, "sizeof_double=%d\n", (int)sizeof(double) );
  fprintf( outputFile, "sizeof_s64Int=%d\n", (int)sizeof(s64Int) );
  fprintf( outputFile, "sizeof_u64Int=%d\n", (int)sizeof(u64Int) );
  //fprintf( outputFile, "sizeof_struct_double_double=%d\n", (int)sizeof(struct{double HPCC_r,HPCC_i;}) );
  fprintf( outputFile, "CommWorldProcs=%d\n", commSize );
  fprintf( outputFile, "MPI_Wtick=%e\n", MPI_Wtick() );
  fprintf( outputFile, "HPL_Tflops=%g\n", params->HPLrdata.Gflops * 1e-3 );
  fprintf( outputFile, "HPL_time=%g\n", params->HPLrdata.time );
  fprintf( outputFile, "HPL_eps=%g\n", params->HPLrdata.eps );
  fprintf( outputFile, "HPL_RnormI=%g\n", params->HPLrdata.RnormI );
  fprintf( outputFile, "HPL_Anorm1=%g\n", params->HPLrdata.Anorm1 );
  fprintf( outputFile, "HPL_AnormI=%g\n", params->HPLrdata.AnormI );
  fprintf( outputFile, "HPL_Xnorm1=%g\n", params->HPLrdata.Xnorm1 );
  fprintf( outputFile, "HPL_XnormI=%g\n", params->HPLrdata.XnormI );
  fprintf( outputFile, "HPL_BnormI=%g\n", params->HPLrdata.BnormI );
  fprintf( outputFile, "HPL_N=%d\n", params->HPLrdata.N );
  fprintf( outputFile, "HPL_NB=%d\n", params->HPLrdata.NB );
  fprintf( outputFile, "HPL_nprow=%d\n", params->HPLrdata.nprow );
  fprintf( outputFile, "HPL_npcol=%d\n", params->HPLrdata.npcol );
  fprintf( outputFile, "HPL_depth=%d\n", params->HPLrdata.depth );
  fprintf( outputFile, "HPL_nbdiv=%d\n", params->HPLrdata.nbdiv );
  fprintf( outputFile, "HPL_nbmin=%d\n", params->HPLrdata.nbmin );
  fprintf( outputFile, "HPL_cpfact=%c\n", params->HPLrdata.cpfact );
  fprintf( outputFile, "HPL_crfact=%c\n", params->HPLrdata.crfact );
  fprintf( outputFile, "HPL_ctop=%c\n", params->HPLrdata.ctop );
  fprintf( outputFile, "HPL_order=%c\n", params->HPLrdata.order );
  fprintf( outputFile, "HPL_dMACH_EPS=%e\n",  HPL_dlamch( HPL_MACH_EPS ) );
  fprintf( outputFile, "HPL_dMACH_SFMIN=%e\n",HPL_dlamch( HPL_MACH_SFMIN ) );
  fprintf( outputFile, "HPL_dMACH_BASE=%e\n", HPL_dlamch( HPL_MACH_BASE ) );
  fprintf( outputFile, "HPL_dMACH_PREC=%e\n", HPL_dlamch( HPL_MACH_PREC ) );
  fprintf( outputFile, "HPL_dMACH_MLEN=%e\n", HPL_dlamch( HPL_MACH_MLEN ) );
  fprintf( outputFile, "HPL_dMACH_RND=%e\n",  HPL_dlamch( HPL_MACH_RND ) );
  fprintf( outputFile, "HPL_dMACH_EMIN=%e\n", HPL_dlamch( HPL_MACH_EMIN ) );
  fprintf( outputFile, "HPL_dMACH_RMIN=%e\n", HPL_dlamch( HPL_MACH_RMIN ) );
  fprintf( outputFile, "HPL_dMACH_EMAX=%e\n", HPL_dlamch( HPL_MACH_EMAX ) );
  fprintf( outputFile, "HPL_dMACH_RMAX=%e\n", HPL_dlamch( HPL_MACH_RMAX ) );
  fprintf( outputFile, "HPL_sMACH_EPS=%e\n",  (double)HPL_slamch( HPL_MACH_EPS ) );
  fprintf( outputFile, "HPL_sMACH_SFMIN=%e\n",(double)HPL_slamch( HPL_MACH_SFMIN ) );
  fprintf( outputFile, "HPL_sMACH_BASE=%e\n", (double)HPL_slamch( HPL_MACH_BASE ) );
  fprintf( outputFile, "HPL_sMACH_PREC=%e\n", (double)HPL_slamch( HPL_MACH_PREC ) );
  fprintf( outputFile, "HPL_sMACH_MLEN=%e\n", (double)HPL_slamch( HPL_MACH_MLEN ) );
  fprintf( outputFile, "HPL_sMACH_RND=%e\n",  (double)HPL_slamch( HPL_MACH_RND ) );
  fprintf( outputFile, "HPL_sMACH_EMIN=%e\n", (double)HPL_slamch( HPL_MACH_EMIN ) );
  fprintf( outputFile, "HPL_sMACH_RMIN=%e\n", (double)HPL_slamch( HPL_MACH_RMIN ) );
  fprintf( outputFile, "HPL_sMACH_EMAX=%e\n", (double)HPL_slamch( HPL_MACH_EMAX ) );
  fprintf( outputFile, "HPL_sMACH_RMAX=%e\n", (double)HPL_slamch( HPL_MACH_RMAX ) );
  fprintf( outputFile, "dweps=%e\n", HPCC_dweps() );
  fprintf( outputFile, "sweps=%e\n", (double)HPCC_sweps() );
  fprintf( outputFile, "HPLMaxProcs=%d\n", params->HPLMaxProc );
  fprintf( outputFile, "HPLMinProcs=%d\n", params->HPLMinProc );
  fprintf( outputFile, "DGEMM_N=%d\n", params->DGEMM_N );
  fprintf( outputFile, "StarDGEMM_Gflops=%g\n",   params->StarDGEMMGflops );
  fprintf( outputFile, "SingleDGEMM_Gflops=%g\n", params->SingleDGEMMGflops );
  fprintf( outputFile, "PTRANS_GBs=%g\n", params->PTRANSrdata.GBs );
  fprintf( outputFile, "PTRANS_time=%g\n", params->PTRANSrdata.time );
  fprintf( outputFile, "PTRANS_residual=%g\n", params->PTRANSrdata.residual );
  fprintf( outputFile, "PTRANS_n=%d\n", params->PTRANSrdata.n );
  fprintf( outputFile, "PTRANS_nb=%d\n", params->PTRANSrdata.nb );
  fprintf( outputFile, "PTRANS_nprow=%d\n", params->PTRANSrdata.nprow );
  fprintf( outputFile, "PTRANS_npcol=%d\n", params->PTRANSrdata.npcol );
  fprintf( outputFile, "MPIRandomAccess_LCG_N=" FSTR64 "\n", params->MPIRandomAccess_LCG_N );
  fprintf( outputFile, "MPIRandomAccess_LCG_time=%g\n", params->MPIRandomAccess_LCG_time );
  fprintf( outputFile, "MPIRandomAccess_LCG_CheckTime=%g\n", params->MPIRandomAccess_LCG_CheckTime );
  fprintf( outputFile, "MPIRandomAccess_LCG_Errors=" FSTR64 "\n", params->MPIRandomAccess_LCG_Errors );
  fprintf( outputFile, "MPIRandomAccess_LCG_ErrorsFraction=%g\n", params->MPIRandomAccess_LCG_ErrorsFraction );
  fprintf( outputFile, "MPIRandomAccess_LCG_ExeUpdates=" FSTR64 "\n", params->MPIRandomAccess_LCG_ExeUpdates );
  fprintf( outputFile, "MPIRandomAccess_LCG_GUPs=%g\n", params->MPIRandomAccess_LCG_GUPs );
  fprintf( outputFile, "MPIRandomAccess_LCG_TimeBound=%g\n", params->MPIRandomAccess_LCG_TimeBound );
  fprintf( outputFile, "MPIRandomAccess_LCG_Algorithm=%d\n", params->MPIRandomAccess_LCG_Algorithm );
  fprintf( outputFile, "MPIRandomAccess_N=" FSTR64 "\n", params->MPIRandomAccess_N );
  fprintf( outputFile, "MPIRandomAccess_time=%g\n", params->MPIRandomAccess_time );
  fprintf( outputFile, "MPIRandomAccess_CheckTime=%g\n", params->MPIRandomAccess_CheckTime );
  fprintf( outputFile, "MPIRandomAccess_Errors=" FSTR64 "\n", params->MPIRandomAccess_Errors );
  fprintf( outputFile, "MPIRandomAccess_ErrorsFraction=%g\n", params->MPIRandomAccess_ErrorsFraction );
  fprintf( outputFile, "MPIRandomAccess_ExeUpdates=" FSTR64 "\n", params->MPIRandomAccess_ExeUpdates );
  fprintf( outputFile, "MPIRandomAccess_GUPs=%g\n", params->MPIRandomAccess_GUPs );
  fprintf( outputFile, "MPIRandomAccess_TimeBound=%g\n", params->MPIRandomAccess_TimeBound );
  fprintf( outputFile, "MPIRandomAccess_Algorithm=%d\n", params->MPIRandomAccess_Algorithm );
  fprintf( outputFile, "RandomAccess_LCG_N=" FSTR64 "\n", params->RandomAccess_LCG_N );
  fprintf( outputFile, "StarRandomAccess_LCG_GUPs=%g\n", params->Star_LCG_GUPs );
  fprintf( outputFile, "SingleRandomAccess_LCG_GUPs=%g\n", params->Single_LCG_GUPs );
  fprintf( outputFile, "RandomAccess_N=" FSTR64 "\n", params->RandomAccess_N );
  fprintf( outputFile, "StarRandomAccess_GUPs=%g\n", params->StarGUPs );
  fprintf( outputFile, "SingleRandomAccess_GUPs=%g\n", params->SingleGUPs );
  fprintf( outputFile, "STREAM_VectorSize=%d\n", params->StreamVectorSize );
  fprintf( outputFile, "STREAM_Threads=%d\n", params->StreamThreads );
  fprintf( outputFile, "StarSTREAM_Copy=%g\n", params->StarStreamCopyGBs );
  fprintf( outputFile, "StarSTREAM_Scale=%g\n", params->StarStreamScaleGBs );
  fprintf( outputFile, "StarSTREAM_Add=%g\n", params->StarStreamAddGBs );
  fprintf( outputFile, "StarSTREAM_Triad=%g\n", params->StarStreamTriadGBs );
  fprintf( outputFile, "SingleSTREAM_Copy=%g\n", params->SingleStreamCopyGBs );
  fprintf( outputFile, "SingleSTREAM_Scale=%g\n", params->SingleStreamScaleGBs );
  fprintf( outputFile, "SingleSTREAM_Add=%g\n", params->SingleStreamAddGBs );
  fprintf( outputFile, "SingleSTREAM_Triad=%g\n", params->SingleStreamTriadGBs );
  fprintf( outputFile, "FFT_N=%d\n", params->FFT_N );
  fprintf( outputFile, "StarFFT_Gflops=%g\n",   params->StarFFTGflops );
  fprintf( outputFile, "SingleFFT_Gflops=%g\n", params->SingleFFTGflops );
  fprintf( outputFile, "MPIFFT_N=" FSTR64 "\n", params->MPIFFT_N );
  fprintf( outputFile, "MPIFFT_Gflops=%g\n", params->MPIFFTGflops );
  fprintf( outputFile, "MPIFFT_maxErr=%g\n", params->MPIFFT_maxErr );
  fprintf( outputFile, "MPIFFT_Procs=%d\n", params->MPIFFT_Procs );
  fprintf( outputFile, "MaxPingPongLatency_usec=%g\n", params->MaxPingPongLatency );
  fprintf( outputFile, "RandomlyOrderedRingLatency_usec=%g\n", params->RandomlyOrderedRingLatency );
  fprintf( outputFile, "MinPingPongBandwidth_GBytes=%g\n", params->MinPingPongBandwidth );
  fprintf( outputFile, "NaturallyOrderedRingBandwidth_GBytes=%g\n", params->NaturallyOrderedRingBandwidth );
  fprintf( outputFile, "RandomlyOrderedRingBandwidth_GBytes=%g\n", params->RandomlyOrderedRingBandwidth );
  fprintf( outputFile, "MinPingPongLatency_usec=%g\n", params->MinPingPongLatency );
  fprintf( outputFile, "AvgPingPongLatency_usec=%g\n", params->AvgPingPongLatency );
  fprintf( outputFile, "MaxPingPongBandwidth_GBytes=%g\n", params->MaxPingPongBandwidth );
  fprintf( outputFile, "AvgPingPongBandwidth_GBytes=%g\n", params->AvgPingPongBandwidth );
  fprintf( outputFile, "NaturallyOrderedRingLatency_usec=%g\n", params->NaturallyOrderedRingLatency );
  fprintf( outputFile, "FFTEnblk=%d\n", params->FFTEnblk );
  fprintf( outputFile, "FFTEnp=%d\n", params->FFTEnp );
  fprintf( outputFile, "FFTEl2size=%d\n", params->FFTEl2size );

#ifdef _OPENMP
  fprintf( outputFile, "M_OPENMP=%ld\n", (long)(_OPENMP) );
#pragma omp parallel
  {
#pragma omp single nowait
    {

      fprintf( outputFile, "omp_get_num_threads=%d\n", omp_get_num_threads() );
      fprintf( outputFile, "omp_get_max_threads=%d\n", omp_get_max_threads() );
      fprintf( outputFile, "omp_get_num_procs=%d\n",   omp_get_num_procs() );
    }
  }
#else
  fprintf( outputFile, "M_OPENMP=%ld\n", -1L );
  fprintf( outputFile, "omp_get_num_threads=%d\n", 0 );
  fprintf( outputFile, "omp_get_max_threads=%d\n", 0 );
  fprintf( outputFile, "omp_get_num_procs=%d\n",   0 );
#endif

  fprintf( outputFile, "MemProc=%g\n", HPCC_MemProc );
  fprintf( outputFile, "MemSpec=%d\n", HPCC_MemSpec );
  fprintf( outputFile, "MemVal=%g\n", HPCC_MemVal );

  for (i = 0; i < MPIFFT_TIMING_COUNT - 1; i++)
    fprintf( outputFile, "MPIFFT_time%d=%g\n", i, params->MPIFFTtimingsForward[i+1] - params->MPIFFTtimingsForward[i] );

  /* CPS: C Preprocessor Symbols */

  i = 0;
#ifdef HPCC_FFT_235
  i = 1;
#endif
  fprintf( outputFile, "CPS_HPCC_FFT_235=%d\n", i );

  i = 0;
#ifdef HPCC_FFTW_ESTIMATE
  i = 1;
#endif
  fprintf( outputFile, "CPS_HPCC_FFTW_ESTIMATE=%d\n", i );

  i = 0;
#ifdef HPCC_MEMALLCTR
  i = 1;
#endif
  fprintf( outputFile, "CPS_HPCC_MEMALLCTR=%d\n", i );

  i = 0;
#ifdef HPL_USE_GETPROCESSTIMES
  i = 1;
#endif
  fprintf( outputFile, "CPS_HPL_USE_GETPROCESSTIMES=%d\n", i );


  i = 0;
#ifdef RA_SANDIA_NOPT
  i = 1;
#endif
  fprintf( outputFile, "CPS_RA_SANDIA_NOPT=%d\n", i );

  i = 0;
#ifdef RA_SANDIA_OPT2
  i = 1;
#endif
  fprintf( outputFile, "CPS_RA_SANDIA_OPT2=%d\n", i );

  i = 0;
#ifdef USING_FFTW
  i = 1;
#endif
  fprintf( outputFile, "CPS_USING_FFTW=%d\n", i );

  fprintf( outputFile, "End of Summary section.%s\n", "" );
  fprintf( outputFile,
            "########################################################################\n" );
  fprintf( outputFile, "End of HPC Challenge tests.\n" );
  fprintf( outputFile, "Current time (%ld) is %s\n",(long)currentTime,ctime(&currentTime));
  fprintf( outputFile,
            "########################################################################\n" );
  END_IO( myRank, outputFile );

  return 0;
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
HPCC_ProcessGrid(int *P, int *Q, MPI_Comm comm) {
  int myRank, commSize;
  int i, p, q, nproc;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  for (nproc = commSize; ; --nproc) { /* this loop makes at most two iterations */

    for (i = (int)sqrt( nproc ); i > 1; --i) {
      q = nproc / i;
      p = nproc / q;
      if (p * q == nproc) {
        *P = p;
        *Q = q;
        return 0;
      }
    }

    /* if the code gets here `nproc' is small or is a prime */

    if (nproc < 20) { /* do 1D grid for small process counts */
      *P = 1;
      *Q = nproc;
      return 0;
    }
  }

  return 0;
}

size_t
HPCC_Memory(MPI_Comm comm) {
  int myRank, commSize;
  int num_threads;
  char memFile[13] = "hpccmemf.txt";
  char buf[HPL_LINE_MAX]; int nbuf = HPL_LINE_MAX;
  char *sVal;
  FILE *f;
  double mult, mval, procMem;
  size_t rv;

  mult = 1.0;
  num_threads = 1;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single nowait
    {
      num_threads = omp_get_num_threads();
    }
  }
#endif

  if (myRank == 0) {
    procMem = 64;

    f = fopen( memFile, "r" );
    if (f) {

      if (fgets( buf, nbuf, f )) {

        if (strncmp( "Total=", buf, 6 ) == 0) {
          mult = 1.0 / commSize;
          sVal = buf + 6;
          HPCC_MemSpec = 1;
        } else if (strncmp( "Thread=", buf, 7 ) == 0) {
          mult = num_threads;
          sVal = buf + 7;
          HPCC_MemSpec = 2;
        } else if (strncmp( "Process=", buf, 8 ) == 0) {
          mult = 1.0;
          sVal = buf + 8;
          HPCC_MemSpec = 3;
        } else
          sVal = NULL;

        if (sVal && 1 == sscanf( sVal, "%lf", &mval )) {
          procMem = mval * mult;
          HPCC_MemVal = mval;
        }
      }

      fclose( f );
    }
  }

  MPI_Bcast( &procMem, 1, MPI_DOUBLE, 0, comm );

  rv = procMem;
  rv *= 1024; rv *= 1024;

  HPCC_MemProc = procMem;

  return rv;
}

int
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
              HPL_T_SWAP *FSWAP, int *TSWAP, int *L1NOTRAN, int *UNOTRAN, int *EQUIL, int *ALIGN, MPI_Comm comm) {
  int nb = 80;
  double memFactor = 0.8;

  *NS = *NBS = *NPQS = *NPFS = *NBMS = *NDVS = *NRFS = *NTPS = *NDHS = 1;

  TEST->thrsh = 16.0;

  *NB = nb;

  *PMAPPIN = HPL_COLUMN_MAJOR;

  HPCC_ProcessGrid( P, Q, comm );

  *N = (int)sqrt( memFactor * (double)(*P * *Q) * (double)(HPCC_Memory( comm ) / sizeof(double)) ) / (2 * nb);
  *N *= 2*nb; /* make N multiple of 2*nb so both HPL and PTRANS see matrix
                 dimension divisible by nb */

  *PF = HPL_RIGHT_LOOKING;

  *NBM = 4;

  *NDV = 2;

  *RF = HPL_CROUT;

  *TP = HPL_1RING_M;

  *DH = 1;

  *FSWAP = HPL_SW_MIX;

  *TSWAP = 64;

  *L1NOTRAN = 0;

  *UNOTRAN = 0;

  *EQUIL = 1;

  *ALIGN = 8;

  return 0;
}

#ifdef XERBLA_MISSING

#ifdef Add_
#define F77xerbla xerbla_
#endif
#ifdef Add__
#define F77xerbla xerbla__
#endif
#ifdef NoChange
#define F77xerbla xerbla
#endif
#ifdef UpCase
#define F77xerbla XERBLA
#endif
#ifdef f77IsF2C
#define F77xerbla xerbla_
#endif

void
F77xerbla(char *srname, F77_INTEGER *info, long srname_len) {
  /*
  int i; char Cname[7];
  for (i = 0; i < 6; i++) Cname[i] = srname[i];
  Cname[6] = 0;
  printf("xerbla(%d)\n", *info);
  */
  printf("xerbla()\n");
  fflush(stdout);
}
#endif

#ifdef HPCC_MEMALLCTR
#define MEM_MAXCNT 7

typedef double Mem_t;

static Mem_t *Mem_base;
static size_t Mem_dsize;

/*
  Each entry can be in one of three states:
  1. Full (holds a block of allocated memory) if:
     ptr != NULL; size > 0; free == 0
  2. Free (holds block of unallocated memory) if:
     ptr != NULL; free = 1
  3  Empty (doesn't hold a block of memory) if:
     ptr == NULL; free = 1
 */
typedef struct {
  Mem_t *Mem_ptr;
  size_t Mem_size;
  int Mem_free;
} Mem_entry_t;

static Mem_entry_t Mem_blocks[MEM_MAXCNT];

static void
HPCC_alloc_set_empty(int idx) {
  int i, n0, n;

  if (MEM_MAXCNT == idx) {
    n0 = 0;
    n = idx;
  } else {
    n0 = idx;
    n = idx + 1;
  }

  /* initialize all blocks to empty */
  for (i = n0; i < n; ++i) {
    Mem_blocks[i].Mem_ptr = (Mem_t *)(NULL);
    Mem_blocks[i].Mem_size = 0;
    Mem_blocks[i].Mem_free = 1;
  }
}

static void
HPCC_alloc_set_free(int idx, Mem_t *dptr, size_t size) {
  Mem_blocks[idx].Mem_ptr = dptr;
  Mem_blocks[idx].Mem_size = size;
  Mem_blocks[idx].Mem_free = 1;
}

int
HPCC_alloc_init(size_t total_size) {
  size_t dsize;

  Mem_dsize = dsize = Mceil( total_size, sizeof(Mem_t) );

  Mem_base = (Mem_t *)malloc( dsize * sizeof(Mem_t) );

  HPCC_alloc_set_empty( MEM_MAXCNT );

  if (Mem_base) {
    HPCC_alloc_set_free( 0, Mem_base, dsize );
    return 0;
  }

  return -1;
}

int
HPCC_alloc_finalize() {
  free( Mem_base );
  HPCC_alloc_set_empty( MEM_MAXCNT );
  return 0;
}

void *
HPCC_malloc(size_t size) {
  size_t dsize, diff_size, cur_diff_size;
  int i, cur_best, cur_free;

  dsize = Mceil( size, sizeof(Mem_t) );

  cur_diff_size = Mem_dsize + 1;
  cur_free = cur_best = MEM_MAXCNT;

  for (i = 0; i < MEM_MAXCNT; ++i) {
    /* skip full spots */
    if (! Mem_blocks[i].Mem_free)
      continue;

    /* find empty spot */
    if (! Mem_blocks[i].Mem_ptr) {
      cur_free = i;
      continue;
    }

    diff_size = Mem_blocks[i].Mem_size - dsize;

    if (Mem_blocks[i].Mem_size >= dsize && diff_size < cur_diff_size) {
      /* a match that's the best (so far) was found */
      cur_diff_size = diff_size;
      cur_best = i;
    }
  }

  /* found a match */
  if (cur_best < MEM_MAXCNT) {
    if (cur_free < MEM_MAXCNT && cur_diff_size > 0) {
      /* create a new free block */
      HPCC_alloc_set_free( cur_free, Mem_blocks[cur_best].Mem_ptr + dsize,
                           cur_diff_size );

      Mem_blocks[cur_best].Mem_size = dsize; /* shrink the best match */
    }

    Mem_blocks[cur_best].Mem_free = 0;

    return (void *)(Mem_blocks[cur_best].Mem_ptr);
  }

  return NULL;
}

void
HPCC_free(void *ptr) {
  Mem_t *dptr = (Mem_t *)ptr;
  int cur_blk = MEM_MAXCNT, made_changes, i, j;

  /* look for the block being freed */
  for (i = 0; i < MEM_MAXCNT; ++i) {
    if (Mem_blocks[i].Mem_free)
      continue;

    if (Mem_blocks[i].Mem_ptr == dptr) {
      cur_blk = i;
      break;
    }
  }

  /* not finding the pointer (including NULL) causes abort */
  if (MEM_MAXCNT == cur_blk) {
    HPL_pabort( __LINE__, "HPCC_free", "Unknown pointer in HPCC_free()." );
  }

  /* double-free causes abort */
  if (1 == Mem_blocks[cur_blk].Mem_free) {
    HPL_pabort( __LINE__, "HPCC_free", "Second call to HPCC_free() with the same pointer." );
  }

  Mem_blocks[cur_blk].Mem_free = 1;

  /* merge as many blocks as possible */
  for (made_changes = 1; made_changes;) {
    made_changes = 0;

    for (i = 0; i < MEM_MAXCNT; ++i) {

      /* empty or full blocks can't be merged */
      if (! Mem_blocks[i].Mem_free || ! Mem_blocks[i].Mem_ptr)
        continue;

      for (j = 0; j < MEM_MAXCNT; ++j) {

        /* empty or occupied blocks can't be merged */
        if (! Mem_blocks[j].Mem_free || ! Mem_blocks[j].Mem_ptr)
          continue;

        if (Mem_blocks[i].Mem_ptr + Mem_blocks[i].Mem_size ==
            Mem_blocks[j].Mem_ptr) {
          Mem_blocks[i].Mem_size += Mem_blocks[j].Mem_size;

          HPCC_alloc_set_empty( j );

          made_changes = 1;
        }
      }
    }
  }
}
#endif


/* Functions */
int
HPCC_TestDGEMM(HPCC_Params *params,
     int doIO, double *UGflops, int *Un, int *Ufailure) {
  int n, lda, ldb, ldc, failure = 1;
  double *a, *b, *c, *x, *y, *z, alpha, beta, sres, cnrm, xnrm;
  double Gflops = 0.0, dn, t0, t1;
  long l_n;
  FILE *outFile;
  int seed_a, seed_b, seed_c, seed_x;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
  }

  n = (int)(sqrt( params->HPLMaxProcMem / sizeof(double) / 3 + 0.25 ) - 0.5);
  if (n < 0) n = -n; /* if 'n' has overflown an integer */
  l_n = n;
  lda = ldb = ldc = n;

  a = HPCC_XMALLOC( double, l_n * l_n );
  b = HPCC_XMALLOC( double, l_n * l_n );
  c = HPCC_XMALLOC( double, l_n * l_n );

  x = HPCC_XMALLOC( double, l_n );
  y = HPCC_XMALLOC( double, l_n );
  z = HPCC_XMALLOC( double, l_n );

  if (! a || ! b || ! c || ! x || ! y || ! z) {
    goto comp_end;
  }

  seed_a = (int)time( NULL );
  dmatgen( n, n, a, n, seed_a );

  seed_b = (int)time( NULL );
  dmatgen( n, n, b, n, seed_b );

  seed_c = (int)time( NULL );
  dmatgen( n, n, c, n, seed_c );

  seed_x = (int)time( NULL );
  dmatgen( n, 1, x, n, seed_x );

  alpha = a[n / 2];
  beta  = b[n / 2];

  t0 = MPI_Wtime();
  HPL_dgemm( HplColumnMajor, HplNoTrans, HplNoTrans, n, n, n, alpha, a, n, b, n, beta, c, n );
  t1 = MPI_Wtime();

  t1 -= t0;
  dn = (double)n;
  if (t1 != 0.0 && t1 != -0.0)
    Gflops = 2.0e-9 * dn * dn * dn / t1;
  else
    Gflops = 0.0;

  cnrm = dnrm_inf( n, n, c, n );
  xnrm = dnrm_inf( n, 1, x, n );

  /* y <- c*x */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, 1.0, c, ldc, x, 1, 0.0, y, 1 );

  /* z <- b*x */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, 1.0, b, ldb, x, 1, 0.0, z, 1 );

  /* y <- alpha * a * z - y */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, alpha, a, lda, z, 1, -1.0, y, 1 );

  dmatgen( n, n, c, n, seed_c );

  /* y <- beta * c_orig * x + y */
  HPL_dgemv( HplColumnMajor, HplNoTrans, n, n, beta, c, ldc, x, 1, 1.0, y, 1 );

  sres = dnrm_inf( n, 1, y, n ) / cnrm / xnrm / n / HPL_dlamch( HPL_MACH_EPS );

  if (doIO) fprintf( outFile, "Scaled residual: %g\n", sres );

  if (sres < params->test.thrsh)
    failure = 0;

  comp_end:

  if (z) HPCC_free( z );
  if (y) HPCC_free( y );
  if (x) HPCC_free( x );
  if (c) HPCC_free( c );
  if (b) HPCC_free( b );
  if (a) HPCC_free( a );

  if (doIO) {
    fflush( outFile );
    fclose( outFile );
  }

  if (UGflops) *UGflops = Gflops;
  if (Un) *Un = n;
  if (Ufailure) *Ufailure = failure;

  return 0;
}

int
HPCC_external_init(int argc, char *argv[], void *extdata) {
  return 0;
}

int
HPCC_external_finalize(int argc, char *argv[], void *extdata) {
  return 0;
}

int
main (int argc, char * argv[]) {
  HPCC_Params p;
  HPCC_Params* params = &p;

  params->HPLMaxProcMem = 1;
  params->test.thrsh = 1.0;

  int myRank, commSize;
  char *outFname;
  void *extdata;
  int rv, errCount, rank, failure = 0;
  double localGflops;
  int n;
  double scl = 1.0 / RAND_MAX;
  FILE *outputFile;
  MPI_Comm comm = MPI_COMM_WORLD;


  MPI_Init( &argc, &argv );

  if (HPCC_external_init( argc, argv, &extdata ))
    goto hpcc_end;

  if (HPCC_Init( params ))
    goto hpcc_end;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &myRank );

  localGflops = 0.0;


  srand(time(NULL));
  scl *= (double)commSize;

  /* select a node at random, but not node 0 (unless there is just one node) */
  if (1 == commSize)
    rank = 0;
  else
    for (rank = 0; ; rank = (int)(scl * rand())) {
      if (rank > 0 && rank < commSize) break;
    }

  MPI_Bcast( &rank, 1, MPI_INT, 0, comm ); /* broadcast the rank selected on node 0 */

  if (myRank == rank) /* if this node has been selected */
    rv = HPCC_TestDGEMM( params, 0 == myRank ? 1 : 0, &localGflops, &n, &failure );

  MPI_Bcast( &rv, 1, MPI_INT, rank, comm ); /* broadcast error code */
  MPI_Bcast( &failure, 1, MPI_INT, rank, comm ); /* broadcast failure indication */
  errCount = rv;
  if (failure) params->Failure = 1;

  /* broadcast result */
  MPI_Bcast( &localGflops, 1, MPI_DOUBLE, rank, comm );
  params->SingleDGEMMGflops = localGflops;

  BEGIN_IO( myRank, params->outFname, outputFile);
  fprintf( outputFile, "Node(s) with error %d\n", errCount );
  fprintf( outputFile, "Node selected %d\n", rank );
  fprintf( outputFile, "Single DGEMM Gflop/s %.6f\n", localGflops );
  END_IO( myRank, outputFile );

  hpcc_end:
  HPCC_Finalize( params );

  HPCC_external_finalize( argc, argv, extdata );

  MPI_Finalize();

  return 0;
}
