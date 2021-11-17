using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DCHEX
{
    public static void dchex(ref double[] r, int ldr, int p, int k, int l, ref double[] z, int ldz,
            int nz, ref double[] c, ref double[] s, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DCHEX updates the Cholesky factorization of a positive definite matrix.
        //
        //  Discussion:
        //
        //    The factorization has the form
        //
        //      A = R' * R
        //
        //    where A is a positive definite matrix of order P.
        //
        //    The updating involves diagonal permutations of the form
        //
        //      E' * A * E
        //
        //    where E is a permutation matrix.  Specifically, given
        //    an upper triangular matrix R and a permutation matrix
        //    E (which is specified by K, L, and JOB), DCHEX determines
        //    an orthogonal matrix U such that
        //
        //      U * R * E = RR,
        //
        //    where RR is upper triangular.  At the user's option, the
        //    transformation U will be multiplied into the array Z.
        //    If A = X'*X, so that R is the triangular part of the
        //    QR factorization of X, then RR is the triangular part of the
        //    QR factorization of X*E, that is, X with its columns permuted.
        //
        //    For a less terse description of what DCHEX does and how
        //    it may be applied, see the LINPACK guide.
        //
        //    The matrix Q is determined as the product U(L-K)*...*U(1)
        //    of plane rotations of the form
        //
        //      (    C(I)       S(I) )
        //      (                    ),
        //      (   -S(I)       C(I) )
        //
        //    where C(I) is real, the rows these rotations operate on
        //    are described below.
        //
        //    There are two types of permutations, which are determined
        //    by the value of JOB.
        //
        //    1, right circular shift.  The columns are rearranged in the order:
        //
        //         1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
        //
        //       U is the product of L-K rotations U(I), where U(I)
        //       acts in the (L-I,L-I+1)-plane.
        //
        //    2, left circular shift: the columns are rearranged in the order
        //
        //         1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
        //
        //       U is the product of L-K rotations U(I), where U(I)
        //       acts in the (K+I-1,K+I)-plane.
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
        //    Input/output, double R[LDR*P].  On input, the upper
        //    triangular factor that is to be updated.  Elements of R below the
        //    diagonal are not referenced.  On output, R has been updated.
        //
        //    Input, int LDR, the leading dimension of the array R.
        //    LDR must be at least P.
        //
        //    Input, int P, the order of the matrix R.
        //
        //    Input, int K, the first column to be permuted.
        //
        //    Input, int L, the last column to be permuted.
        //    L must be strictly greater than K.
        //
        //    Input/output double Z[LDZ*NZ], an array of NZ P-vectors into
        //    which the transformation U is multiplied.  Z is not referenced if NZ = 0.
        //    On output, Z has been updated.
        //
        //    Input, int LDZ, the leading dimension of the array Z.
        //    LDZ must be at least P.
        //
        //    Input, int NZ, the number of columns of the matrix Z.
        //
        //    Output, double C[P], S[P], the cosines and sines of the
        //    transforming rotations.
        //
        //    Input, int JOB, determines the type of permutation.
        //    1, right circular shift.
        //    2, left circular shift.
        //
    {
        int i;
        int ii;
        int il;
        int iu;
        int j;
        int jj;
        int lm1;
        int lmk;
        double t;
        //
        //  Initialize
        //
        lmk = l - k;
        lm1 = l - 1;
        switch (job)
        {
            //
            //  Right circular shift.
            //
            case 1:
            {
                //
                //  Reorder the columns.
                //
                for (i = 1; i <= l; i++)
                {
                    ii = l - i + 1;
                    s[i - 1] = r[ii - 1 + (l - 1) * ldr];
                }

                for (jj = k; jj <= lm1; jj++)
                {
                    j = lm1 - jj + k;
                    for (i = 1; i <= j; i++)
                    {
                        r[i - 1 + j * ldr] = r[i - 1 + (j - 1) * ldr];
                    }

                    r[j + j * ldr] = 0.0;
                }

                for (i = 1; i <= k - 1; i++)
                {
                    ii = l - i + 1;
                    r[i - 1 + (k - 1) * ldr] = s[ii - 1];
                }

                //
                //  Calculate the rotations.
                //
                t = s[0];
                for (i = 1; i <= lmk; i++)
                {
                    BLAS1D.drotg(ref s[i], ref t, ref c[i - 1], ref s[i - 1]);
                    t = s[i];
                }

                r[k - 1 + (k - 1) * ldr] = t;

                for (j = k + 1; j <= p; j++)
                {
                    il = Math.Max(1, l - j + 1);
                    for (ii = il; ii <= lmk; ii++)
                    {
                        i = l - ii;
                        t = c[ii - 1] * r[i - 1 + (j - 1) * ldr] + s[ii - 1] * r[i + (j - 1) * ldr];
                        r[i + (j - 1) * ldr] = c[ii - 1] * r[i + (j - 1) * ldr] - s[ii - 1] * r[i - 1 + (j - 1) * ldr];
                        r[i - 1 + (j - 1) * ldr] = t;
                    }
                }

                //
                //  If required, apply the transformations to Z.
                //
                for (j = 1; j <= nz; j++)
                {
                    for (ii = 1; ii <= lmk; ii++)
                    {
                        i = l - ii;
                        t = c[ii - 1] * z[i - 1 + (j - 1) * ldr] + s[ii - 1] * z[i + (j - 1) * ldr];
                        z[i + (j - 1) * ldr] = c[ii - 1] * z[i + (j - 1) * ldr] - s[ii - 1] * z[i - 1 + (j - 1) * ldr];
                        z[i - 1 + (j - 1) * ldr] = t;
                    }
                }

                break;
            }
            //
            default:
            {
                //
                //  Reorder the columns.
                //
                for (i = 1; i <= k; i++)
                {
                    ii = lmk + i;
                    s[ii - 1] = r[i - 1 + (k - 1) * ldr];
                }

                for (j = k; j <= lm1; j++)
                {
                    for (i = 1; i <= j; i++)
                    {
                        r[i - 1 + (j - 1) * ldr] = r[i - 1 + j * ldr];
                    }

                    jj = j - k + 1;
                    s[jj - 1] = r[j + j * ldr];
                }

                for (i = 1; i <= k; i++)
                {
                    ii = lmk + i;
                    r[i - 1 + (l - 1) * ldr] = s[ii - 1];
                }

                for (i = k + 1; i <= l; i++)
                {
                    r[i - 1 + (l - 1) * ldr] = 0.0;
                }

                //
                //  Reduction loop.
                //
                for (j = k; j <= p; j++)
                {
                    //
                    //  Apply the rotations.
                    //
                    if (j != k)
                    {
                        iu = Math.Min(j - 1, l - 1);

                        for (i = k; i <= iu; i++)
                        {
                            ii = i - k + 1;
                            t = c[ii - 1] * r[i - 1 + (j - 1) * ldr] + s[ii - 1] * r[i + (j - 1) * ldr];
                            r[i + (j - 1) * ldr] = c[ii - 1] * r[i + (j - 1) * ldr]
                                                   - s[ii - 1] * r[i - 1 + (j - 1) * ldr];
                            r[i - 1 + (j - 1) * ldr] = t;
                        }
                    }

                    if (j < l)
                    {
                        jj = j - k + 1;
                        t = s[jj - 1];
                        BLAS1D.drotg(ref r[j - 1 + (j - 1) * ldr], ref t, ref c[jj - 1], ref s[jj - 1]);
                    }
                }

                //
                //  Apply the rotations to Z.
                //
                for (j = 1; j <= nz; j++)
                {
                    for (i = k; i <= lm1; i++)
                    {
                        ii = i - k + 1;
                        t = c[ii - 1] * z[i - 1 + (j - 1) * ldr] + s[ii - 1] * z[i + (j - 1) * ldr];
                        z[i + (j - 1) * ldr] = c[ii - 1] * z[i + (j - 1) * ldr] - s[ii - 1] * z[i - 1 + (j - 1) * ldr];
                        z[i - 1 + (j - 1) * ldr] = t;
                    }
                }

                break;
            }
        }
    }
}