using System;
using Burkardt.BLAS;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static int dgb_fa ( int n, int ml, int mu, ref double[] a, ref int[] pivot )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGB_FA performs a LINPACK-style PLU factorization of a DGB matrix.
        //
        //  Discussion:
        //
        //    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
        //    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
        //    which may be required to store nonzero entries generated during Gaussian
        //    elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2003
        //
        //  Author:
        //
        //    FORTRAN77 original version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N-1.
        //
        //    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.
        //    On output, A has been overwritten by the LU factors.
        //
        //    Output, int PIVOT[N], the pivot vector.
        //
        //    Output, int DGB_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;
        int jz;
        int k;

        int m = ml + mu + 1;
        //
        //  Zero out the initial fill-in columns.
        //
        int j0 = mu + 2;
        int j1 = Math.Min ( n, m ) - 1;

        for ( jz = j0; jz <= j1; jz++ )
        {
            int i0 = m + 1 - jz;
            for ( i = i0; i <= ml; i++ )
            {
                a[i-1+(jz-1)*col] = 0.0;
            }
        }

        jz = j1;
        int ju = 0;

        for ( k = 1; k <= n-1; k++ )
        {
            //
            //  Zero out the next fill-in column.
            //
            jz += 1;
            if ( jz <= n )
            {
                for ( i = 1; i <= ml; i++ )
                {
                    a[i-1+(jz-1)*col] = 0.0;
                }
            }
            //
            //  Find L = pivot index.
            //
            int lm = Math.Min ( ml, n-k );
            int l = m;

            int j;
            for ( j = m+1; j <= m + lm; j++ )
            {
                if ( Math.Abs ( a[l-1+(k-1)*col] ) < Math.Abs ( a[j-1+(k-1)*col] ) )
                {
                    l = j;
                }
            }

            pivot[k-1] = l + k - m;
            switch (a[l-1+(k-1)*col])
            {
                //
                //  Zero pivot implies this column already triangularized.
                //
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("DGB_FA - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + k + "");
                    return k;
            }
            //
            //  Interchange if necessary.
            //
            double t = a[l-1+(k-1)*col];
            a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
            a[m-1+(k-1)*col] = t;
            //
            //  Compute multipliers.
            //
            for ( i = m+1; i <= m+lm; i++ )
            {
                a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
            }
            //
            //  Row elimination with column indexing.
            //
            ju = Math.Max ( ju, mu + pivot[k-1] );
            ju = Math.Min ( ju, n );
            int mm = m;

            for ( j = k+1; j <= ju; j++ )
            {
                l -= 1;
                mm -= 1;

                if ( l != mm )
                {
                    t                 = a[l-1+(j-1)*col];
                    a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
                    a[mm-1+(j-1)*col] = t;
                }
                for ( i = 1; i <= lm; i++ )
                {
                    a[mm+i-1+(j-1)*col] += a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
                }
            }
        }

        pivot[n-1] = n;

        switch (a[m-1+(n-1)*col])
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("DGB_FA - Fatal error!");
                Console.WriteLine("  Zero pivot on step " + n + "");
                return n;
            default:
                return 0;
        }
    }
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
        int jz;
        int k;

        int m = ml + mu + 1;
        int info = 0;
        //
        //  Zero initial fill-in columns.
        //
        int j0 = mu + 2;
        int j1 = Math.Min(n, m) - 1;

        for (jz = j0; jz <= j1; jz++)
        {
            int i0 = m + 1 - jz;
            for (i = i0; i <= ml; i++)
            {
                abd[i - 1 + (jz - 1) * lda] = 0.0;
            }
        }

        jz = j1;
        int ju = 0;
        //
        //  Gaussian elimination with partial pivoting.
        //
        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Zero out the next fill-in column.
            //
            jz += 1;
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
            int lm = Math.Min(ml, n - k);
            int l = BLAS1D.idamax(lm + 1, abd, 1,  + m - 1 + (k - 1) * lda) + m - 1;
            ipvt[k - 1] = l + k - m;
            switch (abd[l - 1 + (k - 1) * lda])
            {
                //
                //  Zero pivot implies this column already triangularized.
                //
                case 0.0:
                    info = k;
                    break;
                //
                default:
                {
                    double t;
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
                    BLAS1D.dscal(lm, t, ref abd, 1,  + m + (k - 1) * lda);
                    //
                    //  Row elimination with column indexing.
                    //
                    ju = Math.Min(Math.Max(ju, mu + ipvt[k - 1]), n);
                    int mm = m;

                    int j;
                    for (j = k + 1; j <= ju; j++)
                    {
                        l -= 1;
                        mm -= 1;
                        t = abd[l - 1 + (j - 1) * lda];
                        if (l != mm)
                        {
                            abd[l - 1 + (j - 1) * lda] = abd[mm - 1 + (j - 1) * lda];
                            abd[mm - 1 + (j - 1) * lda] = t;
                        }

                        BLAS1D.daxpy(lm, t, abd, 1, ref abd, 1,  + m + (k - 1) * lda,  + mm + (j - 1) * lda);
                    }

                    break;
                }
            }

        }

        ipvt[n - 1] = n;

        info = abd[m - 1 + (n - 1) * lda] switch
        {
            0.0 => n,
            _ => info
        };

        return info;
    }
}