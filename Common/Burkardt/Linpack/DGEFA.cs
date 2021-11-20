using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DGEFA
{
    public static int dgefa(ref double[] a, int lda, int n, ref int[] ipvt )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGEFA factors a real general matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 May 2005
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
        //    Input/output, double A[LDA*N].
        //    On intput, the matrix to be factored.
        //    On output, an upper triangular matrix and the multipliers used to obtain
        //    it.  The factorization can be written A=L*U, where L is a product of
        //    permutation and unit lower triangular matrices, and U is upper triangular.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, int DGEFA, singularity indicator.
        //    0, normal value.
        //    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
        //    but it does indicate that DGESL or DGEDI will divide by zero if called.
        //    Use RCOND in DGECO for a reliable indication of singularity.
        //
    {
        int k;
        //
        //  Gaussian elimination with partial pivoting.
        //
        int info = 0;

        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Find L = pivot index.
            //
            int l = BLAS1D.idamax(n - k + 1, a, 1, index:  + (k - 1) + (k - 1) * lda) + k - 1;
            ipvt[k - 1] = l;
            switch (a[l - 1 + (k - 1) * lda])
            {
                //
                //  Zero pivot implies this column already triangularized.
                //
                case 0.0:
                    info = k;
                    continue;
            }

            //
            //  Interchange if necessary.
            //
            double t;
            if (l != k)
            {
                t = a[l - 1 + (k - 1) * lda];
                a[l - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda];
                a[k - 1 + (k - 1) * lda] = t;
            }

            //
            //  Compute multipliers.
            //
            t = -1.0 / a[k - 1 + (k - 1) * lda];

            BLAS1D.dscal(n - k, t, ref a, 1, index: + k + (k - 1) * lda);
            //
            //  Row elimination with column indexing.
            //
            int j;
            for (j = k + 1; j <= n; j++)
            {
                t = a[l - 1 + (j - 1) * lda];
                if (l != k)
                {
                    a[l - 1 + (j - 1) * lda] = a[k - 1 + (j - 1) * lda];
                    a[k - 1 + (j - 1) * lda] = t;
                }

                BLAS1D.daxpy(n - k, t, a, 1, ref a, 1, xIndex:  + k + (k - 1) * lda, yIndex:  + k + (j - 1) * lda);
            }

        }

        ipvt[n - 1] = n;

        info = a[n - 1 + (n - 1) * lda] switch
        {
            0.0 => n,
            _ => info
        };

        return info;
    }
}