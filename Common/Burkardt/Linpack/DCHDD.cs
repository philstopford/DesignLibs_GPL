using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DCHDD
    {
        public static int dchdd(ref double[] r, int ldr, int p, double[] x, ref double[] z, int ldz,
                int nz, double[] y, ref double[] rho, ref double[] c, ref double[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DCHDD downdates an augmented Cholesky decomposition.
            //
            //  Discussion:
            //
            //    DCHDD can also downdate the triangular factor of an augmented QR
            //    decomposition.
            //
            //    Specifically, given an upper triangular matrix R of order P, a
            //    row vector X, a column vector Z, and a scalar Y, DCHDD
            //    determines an orthogonal matrix U and a scalar ZETA such that
            //
            //          (R   Z )     (RR  ZZ)
            //      U * (      )  =  (      ),
            //          (0 ZETA)     ( X   Y)
            //
            //    where RR is upper triangular.
            //
            //    If R and Z have been obtained from the factorization of a least squares
            //    problem, then RR and ZZ are the factors corresponding to the problem
            //    with the observation (X,Y) removed.  In this case, if RHO
            //    is the norm of the residual vector, then the norm of
            //    the residual vector of the downdated problem is
            //    sqrt ( RHO * RHO - ZETA * ZETA ). DCHDD will simultaneously downdate
            //    several triplets (Z, Y, RHO) along with R.
            //
            //    For a less terse description of what DCHDD does and how
            //    it may be applied, see the LINPACK guide.
            //
            //    The matrix U is determined as the product U(1)*...*U(P)
            //    where U(I) is a rotation in the (P+1,I)-plane of the form
            //
            //      ( C(I)      -S(I)    )
            //      (                    ).
            //      ( S(I)       C(I)    )
            //
            //    The rotations are chosen so that C(I) is real.
            //
            //    The user is warned that a given downdating problem may be impossible
            //    to accomplish or may produce inaccurate results.  For example, this
            //    can happen if X is near a vector whose removal will reduce the
            //    rank of R.  Beware.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 June 2009
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
            //    Input/output, double R[LDR*P], the upper triangular matrix that
            //    is to be  downdated.  The part of R below the diagonal is not referenced.
            //
            //    Input, int LDR, the leading dimension of the array R.
            //    LDR must be at least P.
            //
            //    Input, int P, the order of the matrix R.
            //
            //    Input, double X[P], the row vector that is to be removed from R.
            //
            //    Input/output, double Z[LDZ*NZ], an array of NZ P-vectors
            //    which are to be downdated along with R.
            //
            //    Input, int LDZ, the leading dimension of the array Z.
            //    LDZ must be at least P.
            //
            //    Input, int NZ, the number of vectors to be downdated.
            //    NZ may be zero, in which case Z, Y, and RHO are not referenced.
            //
            //    Input, double Y[NZ], the scalars for the downdating of
            //    the vectors Z.
            //
            //    Input/output, double RHO[NZ], the norms of the residual vectors.
            //    On output these have been changed along with R and Z.
            //
            //    Output, double C[P], S[P], the cosines and sines of the
            //    transforming rotations.
            //
            //    Output, int DCHDD, return flag.
            //     0, the entire downdating was successful.
            //    -1, if R could not be downdated.  In this case, all quantities
            //        are left unaltered.
            //     1, if some RHO could not be downdated.  The offending RHO's are
            //        set to -1.
            //
        {
            double a;
            double alpha;
            double azeta;
            double b;
            int i;
            int ii;
            int info;
            int j;
            double norm;
            double scale;
            double t;
            double xx;
            double zeta;
            //
            //  Solve R' * A = X, placing the result in the array S.
            //
            info = 0;
            s[0] = x[0] / r[0 + 0 * ldr];

            for (j = 2; j <= p; j++)
            {
                s[j - 1] = x[j - 1] - BLAS1D.ddot(j - 1, r, 1, s, 1, xIndex: +0 + (j - 1) * ldr);
                s[j - 1] = s[j - 1] / r[j - 1 + (j - 1) * ldr];
            }

            norm = BLAS1D.dnrm2(p, s, 1);

            if (1.0 <= norm)
            {
                info = -1;
                return info;
            }

            alpha = Math.Sqrt(1.0 - norm * norm);
            //
            //  Determine the transformations.
            //
            for (ii = 1; ii <= p; ii++)
            {
                i = p - ii + 1;
                scale = alpha + Math.Abs(s[i - 1]);
                a = alpha / scale;
                b = s[i - 1] / scale;
                norm = Math.Sqrt(a * a + b * b);
                c[i - 1] = a / norm;
                s[i - 1] = b / norm;
                alpha = scale * norm;
            }

            //
            //  Apply the transformations to R.
            //
            for (j = 1; j <= p; j++)
            {
                xx = 0.0;
                for (ii = 1; ii <= j; ii++)
                {
                    i = j - ii + 1;
                    t = c[i - 1] * xx + s[i - 1] * r[i - 1 + (j - 1) * ldr];
                    r[i - 1 + (j - 1) * ldr] = c[i - 1] * r[i - 1 + (j - 1) * ldr] - s[i - 1] * xx;
                    xx = t;
                }
            }

            //
            //  If required, downdate Z and RHO.
            //
            for (j = 1; j <= nz; j++)
            {
                zeta = y[j - 1];
                for (i = 1; i <= p; i++)
                {
                    z[i - 1 + (j - 1) * ldz] = (z[i - 1 + (j - 1) * ldz] - s[i - 1] * zeta) / c[i - 1];
                    zeta = c[i - 1] * zeta - s[i - 1] * z[i - 1 + (j - 1) * ldz];
                }

                azeta = Math.Abs(zeta);

                if (rho[j - 1] < azeta)
                {
                    info = 1;
                    rho[j - 1] = -1.0;
                }
                else
                {
                    rho[j - 1] = rho[j - 1] * Math.Sqrt(1.0 - Math.Pow(azeta / rho[j - 1], 2));
                }
            }

            return info;
        }
    }
}