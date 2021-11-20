using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZCHUD
{
    public static void zchud(ref Complex[] r, int ldr, int p, Complex[] x,
            ref Complex[] z, int ldz, int nz, ref Complex[] y, ref double[] rho,
            ref double[] c, ref Complex[] s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZCHUD updates an augmented Cholesky decomposition.
        //
        //  Discussion:
        //
        //    ZCHUD updates an augmented Cholesky decomposition of the
        //    triangular part of an augmented QR decomposition.  Specifically,
        //    given an upper triangular matrix R of order P, a row vector
        //    X, a column vector Z, and a scalar Y, ZCHUD determines a
        //    unitary matrix U and a scalar ZETA such that
        //
        //           ( R  Z )     ( RR   ZZ  )
        //      U  * (      )  =  (          ),
        //           ( X  Y )     (  0  ZETA )
        //
        //    where RR is upper triangular.  If R and Z have been
        //    obtained from the factorization of a least squares
        //    problem, then RR and ZZ are the factors corresponding to
        //    the problem with the observation (X,Y) appended.  In this
        //    case, if RHO is the norm of the residual vector, then the
        //    norm of the residual vector of the updated problem is
        //    sqrt ( RHO**2 + ZETA**2 ).  ZCHUD will simultaneously update
        //    several triplets (Z,Y,RHO).
        //
        //    For a less terse description of what ZCHUD does and how
        //    it may be applied see the LINPACK guide.
        //
        //    The matrix U is determined as the product U(P)*...*U(1),
        //    where U(I) is a rotation in the (I,P+1) plane of the
        //    form
        //
        //      (          C(I)    S(I) )
        //      (                       ).
        //      ( -Complex.Conjugateg ( S(I) )  C(I) )
        //
        //    The rotations are chosen so that C(I) is real.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //
        //  Parameters:
        //
        //    Input/output, Complex R[LDR*P], the upper triangular matrix
        //    that is to be updated.  The part of R below the diagonal is
        //    not referenced.
        //
        //    Input, int LDR, the leading dimension of R.
        //    P <= LDR.
        //
        //    Input, int P, the order of the matrix.
        //
        //    Input, Complex X[P], the row to be added to R.
        //
        //    Input/output, Complex Z[LDZ*NZ], NZ P-vectors to
        //    be updated with R.
        //
        //    Input, int LDZ, the leading dimension of Z.
        //    P <= LDZ.
        //
        //    Input, int NZ, the number of vectors to be updated.
        //    NZ may be zero, in which case Z, Y, and RHO are not referenced.
        //
        //    Input, Complex Y[NZ], the scalars for updating the vectors Z.
        //
        //    Input/output, double RHO[NZ]; on input, the norms of the residual
        //    vectors that are to be updated.  If RHO(J) is negative, it is
        //    left unaltered.  On output, the updated values.
        //
        //    Output, double C[P]. the cosines of the transforming rotations.
        //
        //    Output, Complex S[P], the sines of the transforming rotations.
        //
    {
        int i;
        int j;
        Complex t;
        //
        //  Update R.
        //
        for (j = 1; j <= p; j++)
        {
            Complex xj = x[j - 1];
            //
            //  Apply the previous rotations.
            //
            for (i = 1; i <= j - 1; i++)
            {
                t = c[i - 1] * r[i - 1 + (j - 1) * ldr] + s[i - 1] * xj;
                xj = c[i - 1] * xj - Complex.Conjugate(s[i - 1]) * r[i - 1 + (j - 1) * ldr];
                r[i - 1 + (j - 1) * ldr] = t;
            }

            //
            //  Compute the next rotation.
            //
            BLAS1Z.zrotg(ref r[+j - 1 + (j - 1) * ldr], xj, ref c[+j - 1], ref s[+j - 1]);
        }

        //
        //  If required, update Z and RHO.
        //
        for (j = 1; j <= nz; j++)
        {
            Complex zeta = y[j - 1];

            for (i = 1; i <= p; i++)
            {
                t = c[i - 1] * z[i - 1 + (j - 1) * ldz] + s[i - 1] * zeta;
                zeta = c[i - 1] * zeta - Complex.Conjugate(s[i - 1]) * z[i - 1 + (j - 1) * ldz];
                z[i - 1 + (j - 1) * ldz] = t;
            }

            double azeta = Complex.Abs(zeta);

            if (azeta == 0.0 || !(0.0 <= rho[j - 1]))
            {
                continue;
            }

            double scale = azeta + rho[j - 1];
            rho[j - 1] = scale * Math.Sqrt(Math.Pow(azeta / scale, 2)
                                           + Math.Pow(rho[j - 1] / scale, 2));
        }
    }

}