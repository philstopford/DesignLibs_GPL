using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DGBFA
    {
        public static int dgbfa(ref double[] abd, int lda, int n, int ml, int mu, ref int[] ipvt )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGBFA factors a real band matrix by elimination.
        //
        //  Discussion:
        //
        //    DGBFA is usually called by DGBCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2005
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
        //    Input/output, double ABD[LDA*N].  On input, the matrix in band
        //    storage.  The columns of the matrix are stored in the columns of ABD
        //    and the diagonals of the matrix are stored in rows ML+1 through
        //    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
        //    and the multipliers which were used to obtain it.  The factorization
        //    can be written A = L*U where L is a product of permutation and unit lower
        //    triangular matrices and U is upper triangular.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //    2*ML + MU + 1 <= LDA is required.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, MU, the number of diagonals below and above the
        //    main diagonal.  0 <= ML < N, 0 <= MU < N.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, integer DGBFA, error flag.
        //    0, normal value.
        //    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
        //      subroutine, but it does indicate that DGBSL will divide by zero if
        //      called.  Use RCOND in DGBCO for a reliable indication of singularity.
        //
        {
            int i;
            int i0;
            int info;
            int j;
            int j0;
            int j1;
            int ju;
            int jz;
            int k;
            int l;
            int lm;
            int m;
            int mm;
            double t;

            m = ml + mu + 1;
            info = 0;
            //
            //  Zero initial fill-in columns.
            //
            j0 = mu + 2;
            j1 = Math.Min(n, m) - 1;

            for (jz = j0; jz <= j1; jz++)
            {
                i0 = m + 1 - jz;
                for (i = i0; i <= ml; i++)
                {
                    abd[i - 1 + (jz - 1) * lda] = 0.0;
                }
            }

            jz = j1;
            ju = 0;
            //
            //  Gaussian elimination with partial pivoting.
            //
            for (k = 1; k <= n - 1; k++)
            {
                //
                //  Zero out the next fill-in column.
                //
                jz = jz + 1;
                if (jz <= n)
                {
                    for (i = 1; i <= ml; i++)
                    {
                        abd[i - 1 + (jz - 1) * lda] = 0.0;
                    }
                }

                //
                //  Find L = pivot index.
                //
                lm = Math.Min(ml, n - k);
                l = BLAS1D.idamax(lm + 1, abd, 1, index:  + m - 1 + (k - 1) * lda) + m - 1;
                ipvt[k - 1] = l + k - m;
                //
                //  Zero pivot implies this column already triangularized.
                //
                if (abd[l - 1 + (k - 1) * lda] == 0.0)
                {
                    info = k;
                }
                //
                //  Interchange if necessary.
                //
                else
                {
                    if (l != m)
                    {
                        t = abd[l - 1 + (k - 1) * lda];
                        abd[l - 1 + (k - 1) * lda] = abd[m - 1 + (k - 1) * lda];
                        abd[m - 1 + (k - 1) * lda] = t;
                    }

                    //
                    //  Compute multipliers.
                    //
                    t = -1.0 / abd[m - 1 + (k - 1) * lda];
                    BLAS1D.dscal(lm, t, ref abd, 1, index:  + m + (k - 1) * lda);
                    //
                    //  Row elimination with column indexing.
                    //
                    ju = Math.Min(Math.Max(ju, mu + ipvt[k - 1]), n);
                    mm = m;

                    for (j = k + 1; j <= ju; j++)
                    {
                        l = l - 1;
                        mm = mm - 1;
                        t = abd[l - 1 + (j - 1) * lda];
                        if (l != mm)
                        {
                            abd[l - 1 + (j - 1) * lda] = abd[mm - 1 + (j - 1) * lda];
                            abd[mm - 1 + (j - 1) * lda] = t;
                        }

                        BLAS1D.daxpy(lm, t, abd, 1, ref abd, 1, xIndex:  + m + (k - 1) * lda, yIndex:  + mm + (j - 1) * lda);
                    }

                }

            }

            ipvt[n - 1] = n;

            if (abd[m - 1 + (n - 1) * lda] == 0.0)
            {
                info = n;
            }

            return info;
        }
    }
}