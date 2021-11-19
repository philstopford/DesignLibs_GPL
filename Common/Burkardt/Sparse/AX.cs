namespace Burkardt.Sparse;

public static class AX
{
    public static void ax ( double[] a, int[] ia, int[] ja, double[] x, ref double[] w, int n, int nz_num, int xIndex = 0, int wIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AX computes A * X for a sparse matrix.
        //
        //  Discussion:
        //
        //    The matrix A is assumed to be sparse.  To save on storage, only
        //    the nonzero entries of A are stored.  For instance, the K-th nonzero
        //    entry in the matrix is stored by:
        //
        //      A(K) = value of entry,
        //      IA(K) = row of entry,
        //      JA(K) = column of entry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2006
        //
        //  Author:
        //
        //    Original C version by Lili Ju.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
        //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
        //    Charles Romine, Henk van der Vorst,
        //    Templates for the Solution of Linear Systems:
        //    Building Blocks for Iterative Methods,
        //    SIAM, 1994,
        //    ISBN: 0898714710,
        //    LC: QA297.8.T45.
        //
        //    Tim Kelley,
        //    Iterative Methods for Linear and Nonlinear Equations,
        //    SIAM, 2004,
        //    ISBN: 0898713528,
        //    LC: QA297.8.K45.
        //
        //    Yousef Saad,
        //    Iterative Methods for Sparse Linear Systems,
        //    Second Edition,
        //    SIAM, 2003,
        //    ISBN: 0898715342,
        //    LC: QA188.S17.
        //
        //  Parameters:
        //
        //    Input, double A[NZ_NUM], the matrix values.
        //
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
        //    of the matrix values.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double W[N], the value of A*X.
        //
        //    Input, int N, the order of the system.
        //
        //    Input, int NZ_NUM, the number of nonzeros.
        //
    {
        int i;
        int k;

        for ( i = 0; i < n; i++ )
        {
            w[(i + wIndex) % w.Length] = 0.0;
        }

        for ( k = 0; k < nz_num; k++ )
        {
            i = ia[k] - 1;
            int j = ja[k] - 1;
            w[(i + wIndex) % w.Length] += a[k] * x[(j + xIndex) % x.Length];
        }
    }
}