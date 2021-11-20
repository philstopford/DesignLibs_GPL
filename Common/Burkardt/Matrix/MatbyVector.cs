﻿using System;

namespace Burkardt.MatrixNS;

public static class MatbyVector
{
    public static double[] mv_gb(int m, int n, int ml, int mu, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MV_GB multiplies a banded matrix by an R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //
        //    Input, double A[(2*ML+MU+1)*N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double MV_GB[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
        }

        for (i = 0; i < m; i++)
        {
            int jlo = Math.Max(0, i - ml);
            int jhi = Math.Min(n - 1, i + mu);
            int j;
            for (j = jlo; j <= jhi; j++)
            {
                b[i] += a[i - j + ml + mu + j * (2 * ml + mu + 1)] * x[j];
            }
        }

        return b;
    }

    public static double[] mv_ge(int m, int n, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MV_GE multiplies a GE matrix by an R8VEC.
        //
        //  Discussion:
        //
        //    The GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    GE storage is used by LINPACK and LAPACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, double A(M,N), the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double MV_GE[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                b[i] += a[i + j * m] * x[j];
            }
        }

        return b;
    }

    public static double[] mv_st(int m, int n, int nz_num, int[] row, int[] col, double[] a,
            double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MV_ST multiplies a sparse triple matrix times a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int NZ_NUM, the number of nonzero values.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
        //    column indices.
        //
        //    Input, double A[NZ_NUM], the nonzero values in the matrix.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Output, double MV_ST[M], the product A*X.
        //
    {
        int i;
        int k;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
        }

        for (k = 0; k < nz_num; k++)
        {
            b[row[k]] += a[k] * x[col[k]];
        }

        return b;
    }
}