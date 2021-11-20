using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DCHUD
{
    public static void dchud(ref double[] r, int ldr, int p, double[] x, ref double[] z, int ldz,
            int nz, double[] y, ref double[] rho, ref double[] c, ref double[] s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DCHUD updates an augmented Cholesky decomposition.
        //
        //  Discussion:
        //
        //    DCHUD can also update the triangular part of an augmented QR
        //    decomposition.
        //
        //    Specifically, given an upper triangular matrix R of order P, a row vector
        //    X, a column vector Z, and a scalar Y, DCHUD determines a unitary matrix
        //    U and a scalar ZETA such that
        //
        //           (R  Z)     (RR   ZZ )
        //      U  * (    )  =  (        ),
        //           (X  Y)     ( 0  ZETA)
        //
        //    where RR is upper triangular.
        //
        //    If R and Z have been obtained from the factorization of a least squares
        //    problem, then RR and ZZ are the factors corresponding to the problem
        //    with the observation (X,Y) appended.  In this case, if RHO is the
        //    norm of the residual vector, then the norm of the residual vector of
        //    the updated problem is sqrt ( RHO * RHO + ZETA * ZETA ).  DCHUD will
        //    simultaneously update several triplets (Z, Y, RHO).
        //
        //    For a less terse description of what DCHUD does and how
        //    it may be applied, see the LINPACK guide.
        //
        //    The matrix U is determined as the product U(P)*...*U(1),
        //    where U(I) is a rotation in the (I,P+1) plane of the form
        //
        //      (     C(I)      S(I) )
        //      (                    ).
        //      (    -S(I)      C(I) )
        //
        //    The rotations are chosen so that C(I) is real.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 June 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
        //    Pete Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input/output, double R[LDR*P], the upper triangular matrix to be
        //    updated.  The part of R below the diagonal is not referenced.
        //    On output, the matrix has been updated.
        //
        //    Input, int LDR, the leading dimension of the array R.
        //    LDR must be at least equal to P.
        //
        //    Input, int P, the order of the matrix R.
        //
        //    Input, double X[P], the row to be added to R.
        //
        //    Input/output, double Z[LDZ*NZ], contains NZ P-vectors
        //    to be updated with R.
        //
        //    Input, int LDZ, the leading dimension of the array Z.
        //    LDZ must be at least P.
        //
        //    Input, int NZ, the number of vectors to be updated.  NZ may be
        //    zero, in which case Z, Y, and RHO are not referenced.
        //
        //    Input, double Y[NZ], the scalars for updating the vectors Z.
        //
        //    Input/output, double RHO[NZ].  On input, the norms of the
        //    residual vectors to be updated.  If RHO(J) is negative, it is left
        //    unaltered.
        //
        //    Output, double C[P], S[P], the cosines and sines of the
        //    transforming rotations.
        //
    {
        int i;
        int j;
        double t;
        //
        //  Update R.
        //
        for (j = 1; j <= p; j++)
        {
            double xj = x[j - 1];
            //
            //  Apply the previous rotations.
            //
            for (i = 1; i <= j - 1; i++)
            {
                t = c[i - 1] * r[i - 1 + (j - 1) * ldz] + s[i - 1] * xj;
                xj = c[i - 1] * xj - s[i - 1] * r[i - 1 + (j - 1) * ldz];
                r[i - 1 + (j - 1) * ldz] = t;
            }

            //
            //  Compute the next rotation.
            //
            BLAS1D.drotg(ref r[j - 1 + (j - 1) * ldr], ref xj, ref c[j - 1], ref s[j - 1]);
        }

        //
        //  If required, update Z and RHO.
        //
        for (j = 1; j <= nz; j++)
        {
            double zeta = y[j - 1];
            for (i = 1; i <= p; i++)
            {
                t = c[i - 1] * z[i - 1 + (j - 1) * ldz] + s[i - 1] * zeta;
                zeta = c[i - 1] * zeta - s[i - 1] * z[i - 1 + (j - 1) * ldz];
                z[i - 1 + (j - 1) * ldz] = t;
            }

            double azeta = Math.Abs(zeta);

            if (azeta == 0.0 || !(0.0 <= rho[j - 1]))
            {
                continue;
            }

            double scale = azeta + rho[j - 1];
            rho[j - 1] = scale * Math.Sqrt(
                Math.Pow(azeta / scale, 2) + Math.Pow(rho[j - 1] / scale, 2));

        }
    }
}