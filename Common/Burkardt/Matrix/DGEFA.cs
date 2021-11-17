using System;
using Burkardt.BLAS;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double dge_det(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGE_DET computes the determinant of a square matrix in DGE storage.
        //
        //  Discussion:
        //
        //    The DGE storage format is used for a general M by N matrix.  A storage
        //    space is made for each logical entry.  The two dimensional logical
        //    array is mapped to a vector, in which storage is by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dongarra, Bunch, Moler, Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input/output, double A[N*N], the matrix to be analyzed.
        //    On output, the matrix has been overwritten by factorization information.
        //
        //    Output, double DGE_DET, the determinant of the matrix.
        //
    {
        double det;
        int i;
        int j;
        int k;
        int l;
        double t;

        det = 1.0;

        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Find L, the index of the pivot row.
            //
            l = k;
            for (i = k + 1; i <= n; i++)
            {
                if (Math.Abs(a[l - 1 + (k - 1) * n]) < Math.Abs(a[i - 1 + (k - 1) * n]))
                {
                    l = i;
                }
            }

            det *= a[l - 1 + (k - 1) * n];

            switch (a[l - 1 + (k - 1) * n])
            {
                case 0.0:
                    return det;
            }

            //
            //  Interchange rows L and K if necessary.
            //
            if (l != k)
            {
                t = a[l - 1 + (k - 1) * n];
                a[l - 1 + (k - 1) * n] = a[k - 1 + (k - 1) * n];
                a[k - 1 + (k - 1) * n] = t;
            }

            //
            //  Normalize the values that lie below the pivot entry A(K,K).
            //
            for (i = k + 1; i <= n; i++)
            {
                a[i - 1 + (k - 1) * n] = -a[i - 1 + (k - 1) * n] / a[k - 1 + (k - 1) * n];
            }

            //
            //  Row elimination with column indexing.
            //
            for (j = k + 1; j <= n; j++)
            {
                if (l != k)
                {
                    t = a[l - 1 + (j - 1) * n];
                    a[l - 1 + (j - 1) * n] = a[k - 1 + (j - 1) * n];
                    a[k - 1 + (j - 1) * n] = t;
                }

                for (i = k + 1; i <= n; i++)
                {
                    a[i - 1 + (j - 1) * n] += a[i - 1 + (k - 1) * n] * a[k - 1 + (j - 1) * n];
                }
            }
        }

        det *= a[n - 1 + (n - 1) * n];

        return det;
    }

    public static double dge_det(int n, double[] a, int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGE_DET computes the determinant of a matrix factored by SGE_FA.
        //
        //  Discussion:
        //
        //    The doubly dimensioned array A is treated as a one dimensional vector,
        //    stored by COLUMNS:
        //
        //      A(0,0), A(1,0), A(2,0), ..., A(N-1,0) // A(1,0), A(1,1), ... A(N-1,1)
        //
        //    Entry A(I,J) is stored as A[I+J*N]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, James Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[N*N], the LU factors computed by DGE_FA.
        //
        //    Input, int PIVOT[N], as computed by DGE_FA.
        //
        //    Output, double DGE_DET, the determinant of the matrix.
        //
    {
        double det;
        int i;

        det = 1.0;

        for (i = 0; i < n; i++)
        {
            det *= a[i + i * n];
            if (pivot[i] != i + 1)
            {
                det = -det;
            }
        }

        return det;
    }

    public static int dge_fa(int n, ref double[] a, ref int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGE_FA factors a general matrix.
        //
        //  Discussion:
        //
        //    DGE_FA is a simplified version of the LINPACK routine SGEFA.
        //
        //    The doubly dimensioned array A is treated as a one dimensional vector,
        //    stored by COLUMNS:
        //
        //      A(0,0), A(1,0), A(2,0), ..., A(N-1,0) // A(1,0), A(1,1), ... A(N-1,1)
        //
        //    Entry A(I,J) is stored as A[I+J*N]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, James Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input/output, double A[N*N], the matrix to be factored.
        //    On output, A contains an upper triangular matrix and the multipliers
        //    which were used to obtain it.  The factorization can be written
        //    A = L * U, where L is a product of permutation and unit lower
        //    triangular matrices and U is upper triangular.
        //
        //    Output, int PIVOT[N], a vector of pivot indices.
        //
        //    Output, int DGE_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the DGE_FA-th step.
        //
    {
        int i;
        int ii;
        int info;
        int j;
        int k;
        int l;
        double t;

        info = 0;

        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Find L, the index of the pivot row.
            //
            l = k;
            for (i = k + 1; i <= n; i++)
            {
                if (Math.Abs(a[l - 1 + (k - 1) * n]) < Math.Abs(a[i - 1 + (k - 1) * n]))
                {
                    l = i;
                }
            }

            pivot[k - 1] = l;
            switch (a[l - 1 + (k - 1) * n])
            {
                //
                //  If the pivot index is zero, the algorithm has failed.
                //
                case 0.0:
                    info = k;
                    Console.WriteLine("");
                    Console.WriteLine("DGE_FA - Warning!");
                    Console.WriteLine("  Zero pivot on step " + info + "");
                    return info;
            }

            //
            //  Interchange rows L and K if necessary.
            //
            if (l != k)
            {
                t = a[l - 1 + (k - 1) * n];
                a[l - 1 + (k - 1) * n] = a[k - 1 + (k - 1) * n];
                a[k - 1 + (k - 1) * n] = t;
            }

            //
            //  Normalize the values that lie below the pivot entry A(K,K).
            //
            for (j = k + 1; j <= n; j++)
            {
                a[j - 1 + (k - 1) * n] = -a[j - 1 + (k - 1) * n] / a[k - 1 + (k - 1) * n];
            }

            //
            //  Row elimination with column indexing.
            //
            for (j = k + 1; j <= n; j++)
            {
                if (l != k)
                {
                    t = a[l - 1 + (j - 1) * n];
                    a[l - 1 + (j - 1) * n] = a[k - 1 + (j - 1) * n];
                    a[k - 1 + (j - 1) * n] = t;
                }

                for (ii = k; ii < n; ii++)
                {
                    a[ii + (j - 1) * n] += a[ii + (k - 1) * n] * a[k - 1 + (j - 1) * n];
                }
            }
        }

        pivot[n - 1] = n;

        switch (a[n - 1 + (n - 1) * n])
        {
            case 0.0:
                info = n;
                Console.WriteLine("");
                Console.WriteLine("DGE_FA - Warning!");
                Console.WriteLine("  Zero pivot on step " + info + "");
                break;
        }

        return info;
    }

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
        int info;
        int j;
        int k;
        int l;
        double t;
        //
        //  Gaussian elimination with partial pivoting.
        //
        info = 0;

        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Find L = pivot index.
            //
            l = BLAS1D.idamax(n - k + 1, a, 1,  + (k - 1) + (k - 1) * lda) + k - 1;
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

            BLAS1D.dscal(n - k, t, ref a, 1,  + k + (k - 1) * lda);
            //
            //  Row elimination with column indexing.
            //
            for (j = k + 1; j <= n; j++)
            {
                t = a[l - 1 + (j - 1) * lda];
                if (l != k)
                {
                    a[l - 1 + (j - 1) * lda] = a[k - 1 + (j - 1) * lda];
                    a[k - 1 + (j - 1) * lda] = t;
                }

                BLAS1D.daxpy(n - k, t, a, 1, ref a, 1,  + k + (k - 1) * lda,  + k + (j - 1) * lda);
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