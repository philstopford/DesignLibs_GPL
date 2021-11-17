using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZCDHD
{
    public static int zcdhd(ref Complex[] r, int ldr, int p, Complex[] x,
            ref Complex[] z, int ldz, int nz, Complex[] y, ref double[] rho,
            ref double[] c, ref Complex[] s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZCHDD downdates an augmented Cholesky decomposition.
        //
        //  Discussion:
        //
        //    zcdhD downdates an augmented Cholesky decomposition or the
        //    triangular factor of an augmented QR decomposition.
        //    Specifically, given an upper triangular matrix R of order P,  a
        //    row vector X, a column vector Z, and a scalar Y, ZCHDD
        //    determines a unitary matrix U and a scalar ZETA such that
        //
        //          ( R   Z  )     ( RR  ZZ )
        //      U * (        )  =  (        ),
        //          ( 0 ZETA )     (  X   Y )
        //
        //    where RR is upper triangular.  If R and Z have been obtained
        //    from the factorization of a least squares problem, then
        //    RR and ZZ are the factors corresponding to the problem
        //    with the observation (X,Y) removed.  In this case, if RHO
        //    is the norm of the residual vector, then the norm of
        //    the residual vector of the downdated problem is
        //      Math.Sqrt ( RHO**2 - ZETA**2 ).
        //    zcdhD will simultaneously downdate several triplets (Z,Y,RHO)
        //    along with R.
        //
        //    For a less terse description of what ZCHDD does and how
        //    it may be applied, see the LINPACK guide.
        //
        //    The matrix U is determined as the product U(1)*...*U(P)
        //    where U(I) is a rotation in the (P+1,I)-plane of the
        //    form
        //
        //      ( C(I)  -Complex.Conjugate ( S(I) ) )
        //      (                       ).
        //      ( S(I)           C(I)   )
        //
        //    The rotations are chosen so that C(I) is real.
        //
        //    The user is warned that a given downdating problem may
        //    be impossible to accomplish or may produce
        //    inaccurate results.  For example, this can happen
        //    if X is near a vector whose removal will reduce the
        //    rank of R.  Beware.
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
        //    Input/output, Complex R[LDR*P]; on input, the upper triangular matrix
        //    that is to be downdated.  On output, the downdated matrix.  The
        //    part of R below the diagonal is not referenced.
        //
        //    Input, int LDR, the leading dimension of R.  P <= LDR.
        //
        //    Input, int P, the order of the matrix.
        //
        //    Input, Complex X(P), the row vector that is to
        //    be removed from R.
        //
        //    Input/output, Complex Z[LDZ*NZ]; on input, an array of NZ
        //    P-vectors which are to be downdated along with R.  On output,
        //    the downdated vectors.
        //
        //    Input, int LDZ, the leading dimension of Z.  P <= LDZ.
        //
        //    Input, int NZ, the number of vectors to be downdated.
        //    NZ may be zero, in which case Z, Y, and R are not referenced.
        //
        //    Input, Complex Y[NZ], the scalars for the downdating
        //    of the vectors Z.
        //
        //    Input/output, double RHO[NZ].  On input, the norms of the residual
        //    vectors that are to be downdated.  On output, the downdated norms.
        //
        //    Output, double C[P], the cosines of the transforming rotations.
        //
        //    Output, Complex S[P], the sines of the transforming rotations.
        //
        //    Output, int ZCHDD:
        //     0, if the entire downdating was successful.
        //    -1, if R could not be downdated.  In this case, all quantities
        //        are left unaltered.
        //     1, if some RHO could not be downdated.  The offending RHO's are
        //        set to -1.
        //
    {
        double a;
        double alpha;
        double azeta;
        Complex b;
        int i;
        int ii;
        int info;
        int j;
        double norm;
        double scale;
        Complex t;
        Complex xx;
        Complex zeta;
        //
        //  Solve the system hermitian(R) * A = X, placing the result in S.
        //
        info = 0;
        s[0] = Complex.Conjugate(x[0]) / Complex.Conjugate(r[0 + 0 * ldr]);

        for (j = 2; j <= p; j++)
        {
            s[j - 1] = Complex.Conjugate(x[j - 1]) - BLAS1Z.zdotc(j - 1, r, 1, s, 1, xIndex: +0 + (j - 1) * ldr);
            s[j - 1] /= Complex.Conjugate(r[j - 1 + (j - 1) * ldr]);
        }

        norm = BLAS1Z.dznrm2(p, s, 1);

        switch (norm)
        {
            case >= 1.0:
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
            scale = alpha + Complex.Abs(s[i - 1]);
            a = alpha / scale;
            b = s[i - 1] / scale;
            norm = Math.Sqrt(a * a + b.Real * b.Real + b.Imaginary * b.Imaginary);
            c[i - 1] = a / norm;
            s[i - 1] = Complex.Conjugate(b) / norm;
            alpha = scale * norm;
        }

        //
        //  Apply the transformations to R.
        //
        for (j = 1; j <= p; j++)
        {
            xx = new Complex(0.0, 0.0);
            for (ii = 1; ii <= j; ii++)
            {
                i = j - ii + 1;
                t = c[i - 1] * xx + s[i - 1] * r[i - 1 + (j - 1) * ldr];
                r[i - 1 + (j - 1) * ldr] = c[i - 1] * r[i - 1 + (j - 1) * ldr] - Complex.Conjugate(s[i - 1]) * xx;
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
                z[i - 1 + (j - 1) * ldz] = (z[i - 1 + (j - 1) * ldz]
                                            - Complex.Conjugate(s[i - 1]) * zeta) / c[i - 1];
                zeta = c[i - 1] * zeta - s[i - 1] * z[i - 1 + (j - 1) * ldz];
            }

            azeta = Complex.Abs(zeta);

            if (rho[j - 1] < azeta)
            {
                info = 1;
                rho[j - 1] = -1.0;
            }
            else
            {
                rho[j - 1] *= Math.Sqrt(1.0 - azeta / rho[j - 1] * (azeta / rho[j - 1]));
            }
        }

        return info;
    }

}