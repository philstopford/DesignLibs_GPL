using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[][] r8mat_to_r8cmat_new(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TO_R8CMAT_NEW copies data from an R8MAT to an R8CMAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a column-major array stored as a vector, so
        //    that element (I,J) of the M by N array is stored in location
        //    I+J*M.
        //
        //    An R8CMAT is a column-major array, storing element (I,J)
        //    as A[J][I], and can be created by a command like:
        //      double **a;
        //      a = r8cmat_new ( m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the data, stored as an R8MAT.
        //
        //    Output, double R8MAT_TO_R8CMAT_NEW[M][N], the data, stored as an R8CMAT.
        //
    {
        int j;

        double[][] b = r8cmat_new(m, n);

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                b[j][i] = a[i + j * m];
            }
        }

        return b;
    }

    public static int r8mat_to_r8plu(int n, double[] a, ref int[] pivot, ref double[] lu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TO_R8PLU factors a general matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    This routine is a simplified version of the LINPACK routine DGEFA.
        //    Fortran conventions are used to index doubly-dimensioned arrays.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[N*N], the matrix to be factored.
        //
        //    Output, int PIVOT[N], a vector of pivot indices.
        //
        //    Output, double LU[N*N], an upper triangular matrix U and the multipliers
        //    L which were used to obtain it.  The factorization can be written
        //    A = L * U, where L is a product of permutation and unit lower
        //    triangular matrices and U is upper triangular.
        //
        //    Output, int R8MAT_TO_R8PLU, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the R8MAT_TO_R8PLU-th step.
        //
    {
        int i;
        int j;
        int k;

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                lu[i + j * n] = a[i + j * n];
            }
        }

        int info = 0;

        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Find L, the index of the pivot row.
            //
            int l = k;
            for (i = k + 1; i <= n; i++)
            {
                if (Math.Abs(lu[l - 1 + (k - 1) * n]) < Math.Abs(lu[i - 1 + (k - 1) * n]))
                {
                    l = i;
                }
            }

            pivot[k - 1] = l;
            switch (lu[l - 1 + (k - 1) * n])
            {
                //
                //  If the pivot index is zero, the algorithm has failed.
                //
                case 0.0:
                    info = k;
                    return info;
            }

            //
            //  Interchange rows L and K if necessary.
            //
            double temp;
            if (l != k)
            {
                temp = lu[l - 1 + (k - 1) * n];
                lu[l - 1 + (k - 1) * n] = lu[k - 1 + (k - 1) * n];
                lu[k - 1 + (k - 1) * n] = temp;
            }

            //
            //  Normalize the values that lie below the pivot entry A(K,K).
            //
            for (i = k + 1; i <= n; i++)
            {
                lu[i - 1 + (k - 1) * n] = -lu[i - 1 + (k - 1) * n] / lu[k - 1 + (k - 1) * n];
            }

            //
            //  Row elimination with column indexing.
            //
            for (j = k + 1; j <= n; j++)
            {
                if (l != k)
                {
                    temp = lu[l - 1 + (j - 1) * n];
                    lu[l - 1 + (j - 1) * n] = lu[k - 1 + (j - 1) * n];
                    lu[k - 1 + (j - 1) * n] = temp;
                }

                for (i = k + 1; i <= n; i++)
                {
                    lu[i - 1 + (j - 1) * n] += lu[i - 1 + (k - 1) * n] * lu[k - 1 + (j - 1) * n];
                }
            }
        }

        pivot[n - 1] = n;

        info = lu[n - 1 + (n - 1) * n] switch
        {
            0.0 => n,
            _ => info
        };

        return info;
    }

    public static double[][] r8mat_to_r8rmat(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TO_R8RMAT copies data from an R8MAT to an R8RMAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a column-major array stored as a vector, so
        //    that element (I,J) of the M by N array is stored in location
        //    I+J*M.
        //
        //    An R8RMAT is a row-major array that was created by a 
        //    command like:
        //
        //      double **a;
        //      a = r8rmat_new ( m, n );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the data, stored as an R8MAT.
        //
        //    Output, double R8RMAT_TO_R8MAT[M][N], the data, stored as an R8RMAT.
        //
    {
        int j;

        double[][] b = r8rmat_new(m, n);

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                b[i][j] = a[i + j * m];
            }
        }

        return b;
    }
        
}